"""LAYER 2 - how accurate is a finite-difference Jacobian?  Error vs step h and solver tolerance tol.

A is computed by perturbing each variable by a relative step h and re-solving the engine:
    forward  dF/dv_j ~ [F(v + s e_j) - F(v)] / s                  truncation error ~ (M2/2) s
    central  dF/dv_j ~ [F(v + s e_j) - F(v - s e_j)] / (2 s)      truncation error ~ (M3/6) s^2
Every F evaluation runs the inner flow-matching solver, which stops at max|r| < tol, so F carries an
error that grows with tol and is divided by s: total error ~ truncation + noise/s -> a U-shaped curve
in log-log axes, with an optimum step
    forward  h* ~ 2 sqrt(eps/M2),    central  h* ~ (3 eps/M3)^(1/3)       (eps = noise level of F).
Two inner-solver starts are studied (see linearize.start_vector):
  'converged' - warm start from the operating point's own solution, as a linearization program does;
  'offset'    - start 0.1 % away from it, so the solver always iterates (classical noisy-F case).

Reference Jacobian: central differences + Richardson extrapolation at tol = 1e-13, checked with two
step sizes. (A complex-step reference is NOT possible for this engine: the T-MATS component code uses
table look-ups, clamps and branches that are not complex-analytic. Complex step is verified on an
analytic test function in tests/test_mathcore.py only.)

Tests hypothesis H3: "the best forward step moves to about 1e-3 when tol is about 1e-6".

Run: python -m eigenmap.layer2_fd [--fast]
Course topic: Jacobians / partial derivatives; Taylor series (truncation error); Eigenvalues
"""
from __future__ import annotations

import numpy as np

from . import config
from . import mathcore as mc
from .common import Checks, finish, write_json
from .linearize import fd_AB, make_F, reference_AB, v0_of
from .opoints import named_points, standard_points

EPS_MACH = float(np.finfo(float).eps)


def eig_error(A: np.ndarray, lam_ref: np.ndarray) -> float:
    """Largest relative eigenvalue error max_i |lam_FD,i - lam_ref,i| / |lam_ref,i| (nearest-neighbour matching)."""
    if not np.all(np.isfinite(A)):
        return float("nan")
    lam, _ = mc.eig2x2(A)
    return float(max(min(abs(l - m) for l in lam) / abs(m) for m in lam_ref))


def noise_level(engine, p, tol: float, start: str, n: int = 6) -> float:
    """eps(tol): relative error of F at this tolerance, measured against tol = 1e-14 at n nearby states.

    Scale for 'relative' = |A| x0 (the size of F for a 100 % change), so eps is directly comparable with h.
    Course topic: Taylor series (error analysis)
    """
    f_tol = make_F(engine, p.sp, tol, p.jac, start)
    f_ref = make_F(engine, p.sp, 1e-14, p.jac)
    v0 = v0_of(p.sp)
    worst = 0.0
    for k in range(n):
        dv = np.zeros(3)
        dv[k % 3] = (1e-7, 1e-5, 1e-3)[k // 3 % 3] * v0[k % 3]
        worst = max(worst, float(np.max(np.abs(f_tol(v0 + dv) - f_ref(v0 + dv)))))
    return worst


def derivative_scales(engine, p, A_ref: np.ndarray) -> tuple:
    """M2, M3: sizes of the 2nd and 3rd derivatives of F along each relative coordinate (tight tol, h = 1e-2).

    Course topic: Taylor series
    """
    f = make_F(engine, p.sp, 1e-13, p.jac)
    v0 = v0_of(p.sp)
    M2 = M3 = 0.0
    for j in range(3):
        s = 1e-2 * v0[j]
        e = np.zeros(3)
        e[j] = s
        fm2, fm1, f0, f1, f2 = (f(v0 - 2 * e), f(v0 - e), f(v0), f(v0 + e), f(v0 + 2 * e))
        d2 = (f1 - 2 * f0 + fm1)            # = s^2 F''  -> in relative units (per unit h^2): divide by 1e-4
        d3 = (f2 - 2 * f1 + 2 * fm1 - fm2) / 2.0
        M2 = max(M2, float(np.max(np.abs(d2))) / 1e-4)
        M3 = max(M3, float(np.max(np.abs(d3))) / 1e-6)
    return M2, M3


def study_point(engine, p, cfg, start: str) -> dict:
    """The full h x tol sweep at one operating point, for one inner-solver start."""
    ref = reference_AB(engine, p.sp, tol=1e-13, jac=p.jac)
    lam_ref, _ = mc.eig2x2(ref["A"])
    # independent check of the reference: plain central difference, h = 1e-6, tol = 1e-14
    A_c, _B = fd_AB(engine, p.sp, 1e-6, "central", 1e-14, p.jac)
    ref_check = eig_error(A_c, lam_ref)
    scale = float(np.max(np.abs(np.hstack([ref["A"] * p.sp.x0, ref["B"] * p.sp.ctrl.Wf]))))   # |F| per 100 % change
    M2, M3 = derivative_scales(engine, p, ref["A"])
    hs = list(cfg.h_grid)
    out = {"id": p.id, "label": p.label, "lam_ref": lam_ref, "A_ref": ref["A"], "ref_rel_diff_A": ref["rel_diff_A"],
           "ref_check_central_1e-6": ref_check, "h": hs, "F_scale": scale, "M2_rel": M2 / scale, "M3_rel": M3 / scale,
           "start": start, "tols": {}}
    for tol in cfg.fd_tols:
        eps_abs = noise_level(engine, p, tol, start)
        eps = max(eps_abs / scale, EPS_MACH)
        row = {"eps_F_rel": eps, "forward": [], "central": []}
        for scheme in ("forward", "central"):
            for h in hs:
                try:
                    A, _ = fd_AB(engine, p.sp, h, scheme, tol, p.jac, start)
                    err = eig_error(A, lam_ref)
                except Exception:          # the engine cannot be solved that far from the point
                    err = float("nan")
                row[scheme].append(err)
        for scheme in ("forward", "central"):
            e = np.array(row[scheme], dtype=float)
            k = int(np.nanargmin(e))
            row[scheme + "_h_opt"] = hs[k]
            row[scheme + "_err_min"] = float(e[k])
            row[scheme + "_U_ratio_left"] = float(np.nanmax(e[:k + 1]) / e[k])
            row[scheme + "_U_ratio_right"] = float(np.nanmax(e[k:]) / e[k])
            # left-branch shape: classical noise/h rise (gradual) or solver dead-zone plateau with a cliff (>10x jump)
            left = e[:k + 1]
            jumps = left[:-1] / np.maximum(left[1:], 1e-300) if len(left) > 1 else np.array([1.0])
            # plateau = a >10x cliff AND, left of the cliff, at least 3 points that stay within a factor 2 of each other
            j = int(np.nanargmax(jumps)) if len(left) > 1 else 0
            before = left[:j + 1]
            flat = len(before) >= 3 and np.nanmax(before[-3:]) / max(np.nanmin(before[-3:]), 1e-300) < 2.0
            row[scheme + "_left_shape"] = ("dead-zone plateau + cliff" if (np.nanmax(jumps) > 10.0 and flat)
                                           else "gradual rise (noise/h-like)")
            # is the minimum well defined? error must rise >= 3x within one decade on BOTH sides
            per_dec = 4 if len(hs) > 20 else 1
            lo_i, hi_i = max(0, k - per_dec), min(len(e) - 1, k + per_dec)
            row[scheme + "_min_sharp"] = bool(e[lo_i] >= 3 * e[k] and e[hi_i] >= 3 * e[k])
        row["forward_h_theory"] = mc.optimal_h_forward(eps, M2 / scale)
        row["central_h_theory"] = mc.optimal_h_central(eps, M3 / scale)
        out["tols"]["%g" % tol] = row
    return out


def h3_verdict(pts: list) -> dict:
    """H3 'best forward step moves to about 1e-3 when tol is about 1e-6', decided by an explicit rule."""
    nums = {}
    per_start = {}
    for start in ("converged", "offset"):
        shift_ok, loc_ok = [], []
        for s in [q for q in pts if q["start"] == start]:
            h6 = s["tols"]["1e-06"]["forward_h_opt"]
            h13 = s["tols"]["1e-13"]["forward_h_opt"]
            nums["%s/%s" % (start, s["id"])] = {"h_opt_fwd_tol1e-13": h13, "h_opt_fwd_tol1e-6": h6, "shift_factor": h6 / h13,
                                                "h_theory_fwd_tol1e-6": s["tols"]["1e-06"]["forward_h_theory"],
                                                "err_min_fwd_tol1e-6": s["tols"]["1e-06"]["forward_err_min"]}
            shift_ok.append(h6 / h13 >= 10.0)
            loc_ok.append(1e-3 / 3.5 <= h6 <= 1e-3 * 3.5)
        per_start[start] = {"shift_everywhere": all(shift_ok), "shift_somewhere": any(shift_ok),
                            "near_1e-3_everywhere": all(loc_ok), "near_1e-3_somewhere": any(loc_ok)}
    loc_nowhere = all(not per_start[k]["near_1e-3_somewhere"] for k in per_start)
    full = [per_start[k]["shift_everywhere"] and per_start[k]["near_1e-3_everywhere"] for k in per_start]
    none = [not per_start[k]["shift_somewhere"] for k in per_start]
    verdict = "supported" if all(full) else ("not supported" if (loc_nowhere or all(none)) else "mixed")
    return {"verdict": verdict, "numbers": nums, "by_start": per_start,
            "rule": "the hypothesis has two parts, (a) the best forward step moves with tol (>= 10x from tol=1e-13 to tol=1e-6 at "
                    "every Layer-2 point) and (b) at tol=1e-6 it lies near 1e-3 (within a factor 3.5). supported = (a) and (b) hold for "
                    "both inner-solver starts; not supported = (b) fails at EVERY point for both starts, or (a) fails everywhere; "
                    "otherwise mixed"}


def run(cfg: config.StudyConfig) -> Checks:
    ck = Checks("layer2_fd")
    engine, _ = standard_points(cfg)
    named = named_points(cfg)
    order = ["idle", "cruise", "max"] if cfg.fast else ["idle", "cruise_idle", "cruise", "max"]
    pts = []
    for start in ("converged", "offset"):
        for key in order:
            p = named[key]
            s = study_point(engine, p, cfg, start)
            s["name"] = key
            pts.append(s)
            if start == "converged":
                ck.add("%s reference eigenvalues vs plain central h=1e-6 at tol 1e-14 (rel)" % p.id,
                       s["ref_check_central_1e-6"] < 1e-6, s["ref_check_central_1e-6"], 1e-6)
            for tol in ("1e-06", "0.0001"):
                r = s["tols"][tol]
                u = min(r["forward_U_ratio_left"], r["forward_U_ratio_right"])
                if start == "converged":     # the configuration a linearization program uses: acceptance check
                    ck.add("%s tol %s: forward-difference error has an interior minimum, rising >=10x on both sides "
                           "(left side shape: %s)" % (p.id, tol, r["forward_left_shape"]), u >= 10.0, u, 10.0)
                else:                        # extra experiment: recorded, reported in RESULTS_FACTS, not a pass/fail item
                    print("INFO layer2_fd: %s [offset start] tol %s U-ratio %.3g" % (p.id, tol, u))
    h3 = h3_verdict(pts)
    # how the best forward step scales with tol (log-log slope over tol = 1e-9 .. 1e-4): 0.5 = classical sqrt(noise) law
    scaling = {}
    for s in pts:
        tl = [t for t in sorted(s["tols"], key=float) if float(t) >= 1e-9]
        x = np.log10([float(t) for t in tl])
        y = np.log10([s["tols"][t]["forward_h_opt"] for t in tl])
        scaling["%s/%s" % (s["start"], s["id"])] = float(np.polyfit(x, y, 1)[0])
    h3["h_opt_vs_tol_loglog_slope"] = scaling
    print("H3 best h shifts with tol: %s" % h3["verdict"])
    for k, v in h3["numbers"].items():
        print("   %-22s best forward h: tol 1e-13 -> %.1e, tol 1e-6 -> %.1e (theory %.1e)" % (
            k, v["h_opt_fwd_tol1e-13"], v["h_opt_fwd_tol1e-6"], v["h_theory_fwd_tol1e-6"]))
    write_json("layer2_fd", {"points": pts, "verdicts": {"H3_best_h_shifts_with_tol": h3}, "fast": cfg.fast,
                             "h_definition": "relative step: perturbation of variable j = h * |v0_j|, v = [N_LP, N_HP, W_f]",
                             "error_definition": "max over the two eigenvalues of |lam_FD - lam_ref| / |lam_ref|"}, ck)
    return ck


if __name__ == "__main__":
    c = config.from_argv()
    c.ensure_dirs()
    finish(run(c))
