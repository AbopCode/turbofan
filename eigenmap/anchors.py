"""Published-example anchor: turboshaft JetCat SPT5 (illustration).

This module checks the OWN math core (`eigenmap.mathcore`) on three small 2x2
state-space models that are printed in the open literature, so that every
number can be recomputed by hand:

    M. Pakmehr and T. Yucelen, "Model Reference Adaptive Control of Systems with
    Gain Scheduled Reference Models", arXiv:1403.3738, Sec. 3.3, eqs. (57)-(59).

IMPORTANT - what this example is and is not:
  * It is a small TURBOSHAFT engine (JetCat SPT5), used as an ILLUSTRATION only.
    It is NOT the AGTF30 geared turbofan studied in the rest of the project and
    no number from here may be quoted as a turbofan result.
  * The states are unlabeled / normalized in the source, the time unit is not
    stated ("model time units" below) and the matrices are printed to 1-2
    significant figures only.
  * The matrices below were transcribed from the project brief; this module does
    not re-open the paper (no network at run time).

Run:  python -m eigenmap.anchors [--fast]     ->  outputs/data/anchors.json

Course topic: Eigenvalues / Diagonalization / Stability / Linear ODE systems /
Gauss-Seidel / Matrix norms and conditioning
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import scipy.linalg

from eigenmap import mathcore as mc
from eigenmap.config import DATA, StudyConfig, from_argv

MODULE = "anchors"
LABEL = "turboshaft JetCat SPT5 (illustration)"
DISCLAIMER = (LABEL + ": published small-turboshaft example used only to illustrate and check the "
              "mathematics; it is NOT the AGTF30 turbofan. States unlabeled/normalized, time unit not "
              "stated (model time units), matrices printed to 1-2 significant figures.")
SOURCE = ("M. Pakmehr and T. Yucelen, 'Model Reference Adaptive Control of Systems with Gain Scheduled "
          "Reference Models', arXiv:1403.3738, Sec. 3.3, eqs. (57)-(59)")

# ---- the three printed operating points (alpha = scheduling parameter of the source) ----
ORDER = ("idle", "mid", "cruise")
POINTS = {
    "idle": {"alpha": 0.3361, "A": [[-0.38, -0.0008], [0.26, -0.34]], "B": [[0.7, 0.0], [0.1, -0.0024]]},
    "mid": {"alpha": 0.6473, "A": [[-0.85, 0.032], [0.32, -0.64]], "B": [[1.0, 0.0], [0.17, -0.011]]},
    "cruise": {"alpha": 0.8818, "A": [[-1.7, 0.1], [0.6, -1.1]], "B": [[1.2, 0.0], [0.3, -0.023]]},
}
INPUTS = ("fuel flow", "propeller pitch")           # columns of B

# ---- numbers quoted in the project brief (rounded as printed there) ----
BRIEF = {
    "idle": {"lam": (-0.346, -0.374), "tau": (2.89, 2.67), "kappaV": 18.9},
    "mid": {"lam": (-0.599, -0.891), "tau": (1.67, 1.12), "kappaV": 2.39},
    "cruise": {"lam": (-1.013, -1.787), "tau": (0.99, 0.56), "kappaV": 1.84},
}
BRIEF_GS = ((0.70588, 0.65775), (0.74457, 0.67886), (0.74582, 0.67954), (0.74586, 0.67956))
BRIEF_GS_EXACT = (0.74586, 0.67956)
BRIEF_GS_RHO = 0.032
BRIEF_SS_CRUISE = (0.746, 0.680)
BRIEF_FLIP = {"a12_printed": -0.0008, "a12_changed": -0.0018, "re": -0.36, "im": 0.008}
BRIEF_MC = {"E_norm": 1e-3, "idle": 0.016, "cruise": 0.0012, "ratio": 13.0}

# ---- tolerances of CONTRACT.md section 9 (published turboshaft example) ----
TOL_ABS = 2e-3            # eigenvalues / tau: absolute
TOL_TAU_PRINTED = 5e-3    # half a unit of the brief's last printed digit of tau (2 decimals)
TOL_KAPPA_REL = 0.01      # kappa(V): 1 %
TOL_GS = 0.5e-5           # Gauss-Seidel table: 5 decimal places
TOL_MC_REL = 0.15         # Monte Carlo worst shift: 15 % of the brief's rounded figure
MC_RATIO_RANGE = (10.0, 16.0)
TOL_LIB = 1e-8            # own routine vs numpy / scipy

# ---- plot grid for the step responses (model time units; not a study tunable) ----
T_END = 20.0
T_PER_UNIT = 10           # samples per model time unit -> 201 points


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _plain(obj):
    """Convert numpy / complex objects into plain JSON types (complex -> [re, im]).

    Course topic: (none - project plumbing)
    """
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _plain(obj.tolist())
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (complex, np.complexfloating)):
        return [float(obj.real), float(obj.imag)]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        return float(obj)
    return obj


class Checks:
    """Collects the PASS/FAIL items (printed AND stored in the JSON `checks` list).

    Course topic: (none - project plumbing)
    """

    def __init__(self) -> None:
        self.items: list[dict] = []

    def add(self, claim: str, ok: bool, value, limit) -> bool:
        """Print `PASS anchors: <claim> (<number> vs <limit>)` and store the item.

        Course topic: (none - project plumbing)
        """
        ok = bool(ok)
        name = f"{LABEL}: {claim}"
        value, limit = _plain(value), _plain(limit)
        fmt = lambda v: f"{v:.6g}" if isinstance(v, float) else str(v)   # noqa: E731
        print(f"{'PASS' if ok else 'FAIL'} {MODULE}: {name} ({fmt(value)} vs {fmt(limit)})")
        self.items.append({"name": name, "pass": ok, "value": value, "limit": limit})
        return ok


def _np_eig_sorted(A: np.ndarray) -> np.ndarray:
    """numpy eigenvalues in the mathcore ordering (slow first / +imag first).

    Course topic: Eigenvalues
    """
    w = np.linalg.eigvals(A).astype(complex)
    if np.any(w.imag != 0.0):
        return np.array(sorted(w, key=lambda z: -z.imag))
    return np.array(sorted(w, key=lambda z: (abs(z.real), -z.real)))


# ---------------------------------------------------------------------------
# 1. eigenvalues, time constants, eigenvector conditioning, step response
# ---------------------------------------------------------------------------
def analyse_point(key: str, chk: Checks, rounding_only: list) -> dict:
    """All eigen-quantities of one printed operating point, with brief and numpy checks.

    Course topic: Eigenvalues
    """
    A = np.array(POINTS[key]["A"], dtype=float)
    B = np.array(POINTS[key]["B"], dtype=float)
    brief = BRIEF[key]
    lam, info = mc.eig2x2(A)
    V = mc.eigvec2x2(A, lam)
    kap = mc.kappa2(V)
    tau_slow, tau_fast = mc.time_constants(lam)
    P, D, Pinv = mc.diagonalize(A)

    # --- against the brief ---
    d_lam = float(np.max(np.abs(lam.real - np.array(brief["lam"]))))
    chk.add(f"{key}: eigenvalues {lam[0].real:.5f} / {lam[1].real:.5f} vs brief "
            f"{brief['lam'][0]} / {brief['lam'][1]} (max abs diff)", d_lam <= TOL_ABS and np.all(lam.imag == 0.0),
            d_lam, TOL_ABS)
    for name, val, ref in (("tau_slow", tau_slow, brief["tau"][0]), ("tau_fast", tau_fast, brief["tau"][1])):
        diff = abs(val - ref)
        if diff <= TOL_ABS:
            chk.add(f"{key}: {name} {val:.5f} vs brief {ref} (abs diff, model time units)", True, diff, TOL_ABS)
        else:
            # not within 2e-3; does it at least round to the printed value?
            ok = diff <= TOL_TAU_PRINTED
            if ok:
                rounding_only.append({"point": key, "quantity": name, "value": val, "brief": ref,
                                      "abs_diff": diff, "contract_tol": TOL_ABS})
            chk.add(f"{key}: {name} {val:.5f} matches brief {ref} ONLY to its printed rounding "
                    f"(abs diff exceeds the 2e-3 tolerance; brief gives 2 decimals)", ok, diff,
                    f"{TOL_TAU_PRINTED} (printed rounding; contract 2e-3 NOT met)")
    rel_k = abs(kap - brief["kappaV"]) / brief["kappaV"]
    chk.add(f"{key}: kappa(V) {kap:.4f} vs brief {brief['kappaV']} (relative diff; unit-norm eigenvector "
            "columns, matrices as printed)", rel_k <= TOL_KAPPA_REL, rel_k, TOL_KAPPA_REL)

    # --- own routines against numpy / scipy ---
    lam_np = _np_eig_sorted(A)
    w_np, V_np = np.linalg.eig(A)
    kap_np = float(np.linalg.cond(V_np, 2))
    err_lam = float(np.max(np.abs(lam - lam_np)))
    err_kap = max(abs(kap - kap_np), abs(kap - float(np.linalg.cond(V, 2)))) / kap_np
    err_rec = float(np.max(np.abs(P @ D @ Pinv - A)))
    chk.add(f"{key}: own eig2x2 vs numpy.linalg.eigvals", err_lam <= TOL_LIB, err_lam, TOL_LIB)
    chk.add(f"{key}: own kappa2 vs numpy.linalg.cond (own V and numpy V)", err_kap <= TOL_LIB, err_kap, TOL_LIB)
    chk.add(f"{key}: own diagonalization P D P^-1 reproduces A", err_rec <= TOL_LIB, err_rec, TOL_LIB)
    stable = mc.is_stable2x2(A)
    chk.add(f"{key}: stable (tr < 0, det > 0, Re(lambda) < 0, diag(A) < 0)",
            stable and np.all(lam.real < 0.0) and A[0, 0] < 0.0 and A[1, 1] < 0.0,
            float(np.max(lam.real)), "< 0")

    # --- unit fuel-flow step: dx(t) = A^-1 (e^{At} - I) B du ---
    t = np.array([i / T_PER_UNIT for i in range(int(round(T_END * T_PER_UNIT)) + 1)])
    du_fuel = [1.0, 0.0]
    dx = mc.step_response(A, B, du_fuel, t)
    gain_fuel = mc.steady_gain(A, B, du_fuel)
    gain_pitch = mc.steady_gain(A, B, [0.0, 1.0])
    aug = np.zeros((3, 3))
    aug[:2, :2], aug[:2, 2] = A, B[:, 0]
    ref = np.array([scipy.linalg.expm(aug * tk)[:2, 2] for tk in t])
    err_step = float(np.max(np.abs(dx - ref)))
    err_gain = float(np.max(np.abs(gain_fuel + np.linalg.solve(A, B[:, 0]))))
    chk.add(f"{key}: own step_response vs scipy.linalg.expm (unit fuel step, {t.size} times)",
            err_step <= TOL_LIB, err_step, TOL_LIB)
    chk.add(f"{key}: own steady_gain vs numpy.linalg.solve", err_gain <= TOL_LIB, err_gain, TOL_LIB)

    return {
        "label": f"{key} (alpha = {POINTS[key]['alpha']}) - {LABEL}",
        "alpha": POINTS[key]["alpha"], "A": A, "B": B, "inputs": list(INPUTS),
        "tr": info["tr"], "det": info["det"], "disc": info["disc"],
        "eig": {"own": lam, "numpy": lam_np, "brief": list(brief["lam"]), "max_abs_diff_brief": d_lam},
        "tau": {"slow": tau_slow, "fast": tau_fast, "brief": list(brief["tau"]), "unit": "model time units"},
        "V": V, "kappaV": {"own": kap, "numpy": kap_np, "brief": brief["kappaV"],
                           "form": "matrices as printed (normalized states), unit 2-norm eigenvector columns"},
        "norm2_A": mc.norm2(A),
        "stable": stable, "diag_dominant": mc.is_diag_dominant(A), "gs_spectral_radius": mc.gs_spectral_radius_2x2(A),
        "steady_gain": {"fuel flow": gain_fuel, "propeller pitch": gain_pitch},
        "step_fuel": {"du": 1.0, "input": "fuel flow", "t": t, "dx1": dx[:, 0], "dx2": dx[:, 1],
                      "unit_t": "model time units", "formula": "dx(t) = A^-1 (e^{At} - I) B du"},
    }


# ---------------------------------------------------------------------------
# 2. Gauss-Seidel on the cruise steady-state equation
# ---------------------------------------------------------------------------
def analyse_gauss_seidel(chk: Checks) -> dict:
    """Gauss-Seidel for A_cruise x = b, b = -B_fuel * du (du = 1): the fuel-step steady state.

    Course topic: Gauss-Seidel
    """
    A = np.array(POINTS["cruise"]["A"], dtype=float)
    B = np.array(POINTS["cruise"]["B"], dtype=float)
    b = -B[:, 0] * 1.0                               # 0 = A x + B du  ->  A x = -B du
    x0 = np.zeros(2)
    x, hist, info = mc.gauss_seidel(A, b, x0=x0, tol=1e-10, max_iter=500)
    exact = mc.steady_gain(A, B, [1.0, 0.0])         # own 2x2 inverse
    exact_np = np.linalg.solve(A, b)
    rho = mc.gs_spectral_radius_2x2(A)
    errs = [float(np.max(np.abs(h - exact))) for h in hist]
    table = []
    for k, h in enumerate(hist):
        row = {"k": k, "x1": h[0], "x2": h[1], "err_inf": errs[k],
               "err_ratio": (errs[k] / errs[k - 1]) if k >= 2 and errs[k - 1] > 0.0 else None}
        if 1 <= k <= len(BRIEF_GS):
            row["brief"] = list(BRIEF_GS[k - 1])
        table.append(row)
    for k in range(1, len(BRIEF_GS) + 1):
        d = float(np.max(np.abs(hist[k] - np.array(BRIEF_GS[k - 1]))))
        chk.add(f"Gauss-Seidel cruise, iterate {k}: {hist[k][0]:.5f} / {hist[k][1]:.5f} vs brief "
                f"{BRIEF_GS[k - 1][0]} / {BRIEF_GS[k - 1][1]} (5 d.p.)", d <= TOL_GS * (1 + 1e-9), d, TOL_GS)
    d = float(np.max(np.abs(exact - np.array(BRIEF_GS_EXACT))))
    chk.add(f"Gauss-Seidel cruise, exact solution {exact[0]:.5f} / {exact[1]:.5f} vs brief (5 d.p.)",
            d <= TOL_GS * (1 + 1e-9), d, TOL_GS)
    d = float(np.max(np.abs(exact - np.array(BRIEF_SS_CRUISE))))
    chk.add(f"cruise unit-fuel-step steady state [{exact[0]:.4f}, {exact[1]:.4f}] vs brief "
            f"{list(BRIEF_SS_CRUISE)}", d <= TOL_ABS, d, TOL_ABS)
    chk.add(f"Gauss-Seidel spectral radius |a12 a21/(a11 a22)| = {rho:.5f} vs brief {BRIEF_GS_RHO}",
            abs(rho - BRIEF_GS_RHO) <= TOL_ABS, abs(rho - BRIEF_GS_RHO), TOL_ABS)
    ratio = errs[2] / errs[1]
    chk.add("Gauss-Seidel observed error ratio per sweep equals the spectral radius",
            abs(ratio - rho) <= 1e-6, abs(ratio - rho), 1e-6)
    d = float(np.max(np.abs(x - exact_np)))
    chk.add(f"Gauss-Seidel converged ({info['iters']} sweeps, tol 1e-10) to numpy.linalg.solve",
            info["converged"] and d <= TOL_LIB, d, TOL_LIB)
    DL, U = np.tril(A), np.triu(A, 1)
    rho_np = float(np.max(np.abs(np.linalg.eigvals(-np.linalg.solve(DL, U)))))
    chk.add("own gs_spectral_radius_2x2 vs numpy eigenvalues of -(D+L)^-1 U", abs(rho - rho_np) <= TOL_LIB,
            abs(rho - rho_np), TOL_LIB)
    chk.add("cruise matrix is strictly diagonally dominant (sufficient for Gauss-Seidel)",
            mc.is_diag_dominant(A), mc.is_diag_dominant(A), True)
    return {
        "label": f"Gauss-Seidel on the cruise point - {LABEL}",
        "equation": "A_cruise x = b, b = -B_fuel * du, du = 1 (steady state of a unit fuel-flow step)",
        "A": A, "b": b, "x0": x0, "tol": 1e-10, "exact_own": exact, "exact_numpy": exact_np,
        "brief_exact": list(BRIEF_GS_EXACT), "spectral_radius": rho, "spectral_radius_numpy": rho_np,
        "brief_spectral_radius": BRIEF_GS_RHO, "observed_error_ratio": ratio,
        "diag_dominant": mc.is_diag_dominant(A), "iters": info["iters"], "converged": info["converged"],
        "residual_inf": info["residual_inf"], "table": table,
    }


# ---------------------------------------------------------------------------
# 3. real -> complex flip at idle
# ---------------------------------------------------------------------------
def analyse_flip(chk: Checks) -> dict:
    """Changing a12 of A_idle from -0.0008 to -0.0018 turns the real pair into a complex pair.

    The eigenvalues are real when disc = (a11 - a22)^2 + 4 a12 a21 >= 0; at idle disc is
    only 7.7e-4, so a change of 1e-3 in one entry makes it negative.

    Course topic: Eigenvalues
    """
    A = np.array(POINTS["idle"]["A"], dtype=float)
    A2 = A.copy()
    A2[0, 1] = BRIEF_FLIP["a12_changed"]
    E = A2 - A
    lam, info = mc.eig2x2(A)
    lam2, info2 = mc.eig2x2(A2)
    lam2_np = _np_eig_sorted(A2)
    a12_crit = -(A[0, 0] - A[1, 1]) ** 2 / (4.0 * A[1, 0])        # disc = 0
    # smallest discriminant inside the box of numbers that round to the printed entries
    # (a11, a22, a21 printed to 2 decimals: +/- 0.005; a12 printed to 4 decimals: +/- 0.00005)
    gap_min = abs(A[0, 0] - A[1, 1]) - 0.01
    disc_box_min = gap_min ** 2 + 4.0 * (A[0, 1] - 0.00005) * (A[1, 0] + 0.005)
    chk.add("idle as printed has a REAL eigenvalue pair (disc > 0)", info["disc"] > 0.0 and np.all(lam.imag == 0.0),
            info["disc"], "> 0")
    chk.add(f"idle with a12 {BRIEF_FLIP['a12_printed']} -> {BRIEF_FLIP['a12_changed']}: COMPLEX pair (disc < 0)",
            info2["disc"] < 0.0 and lam2[0].imag > 0.0, info2["disc"], "< 0")
    chk.add(f"flipped pair real part {lam2[0].real:.5f} vs brief {BRIEF_FLIP['re']}",
            abs(lam2[0].real - BRIEF_FLIP["re"]) <= TOL_ABS, abs(lam2[0].real - BRIEF_FLIP["re"]), TOL_ABS)
    chk.add(f"flipped pair imaginary part +/-{lam2[0].imag:.5f} vs brief +/-{BRIEF_FLIP['im']} "
            "(brief gives 1 significant figure)", abs(lam2[0].imag - BRIEF_FLIP["im"]) <= 5e-4,
            abs(lam2[0].imag - BRIEF_FLIP["im"]), 5e-4)
    err = float(np.max(np.abs(lam2 - lam2_np)))
    chk.add("flipped pair: own eig2x2 vs numpy.linalg.eigvals", err <= TOL_LIB, err, TOL_LIB)
    return {
        "label": f"real-to-complex flip at idle - {LABEL}",
        "A_printed": A, "A_changed": A2, "entry": "a12", "a12_printed": BRIEF_FLIP["a12_printed"],
        "a12_changed": BRIEF_FLIP["a12_changed"], "a12_critical_disc_zero": a12_crit,
        "E_norm2": mc.norm2(E), "E_norm2_over_A_norm2": mc.norm2(E) / mc.norm2(A),
        "eig_printed": lam, "eig_changed": lam2, "eig_changed_numpy": lam2_np,
        "disc_printed": info["disc"], "disc_changed": info2["disc"], "eig_shift": mc.eig_shift(A, E),
        "brief": {"re": BRIEF_FLIP["re"], "im": BRIEF_FLIP["im"]},
        "disc_min_in_printed_rounding_box": disc_box_min,
        "note": ("The change |delta a12| = 1e-3 is smaller than the rounding half-width (5e-3) of the three "
                 "entries printed to 2 decimals, but it is LARGER than the rounding half-width (5e-5) of a12 "
                 "itself as printed (-0.0008). disc_min_in_printed_rounding_box is the smallest discriminant "
                 "over the (closed) box of all values that round to the printed matrix; it is about -1e-6, "
                 "reached only at the extreme corner of that box, so the printed digits only just decide "
                 "that the idle eigenvalue pair is real."),
    }


# ---------------------------------------------------------------------------
# 4. Monte Carlo eigenvalue shifts and Bauer-Fike bounds
# ---------------------------------------------------------------------------
def analyse_monte_carlo(cfg: StudyConfig, chk: Checks) -> dict:
    """Worst eigenvalue shift under random perturbations with ||E||_2 = cfg.E_norm.

    For every seed the SAME perturbation matrices are applied to the three points
    (common random numbers), so the idle / cruise ratio compares like with like.

    Course topic: Matrix norms and conditioning
    """
    e_norm = float(cfg.E_norm)
    n = int(cfg.n_mc)
    mats = {k: np.array(POINTS[k]["A"], dtype=float) for k in ORDER}
    kappa = {k: mc.kappa2(mc.eigvec2x2(mats[k], mc.eig2x2(mats[k])[0])) for k in ORDER}
    bound = {k: kappa[k] * e_norm for k in ORDER}
    brief_applies = abs(e_norm - BRIEF_MC["E_norm"]) <= 1e-15
    per_seed = []
    quantiles = (0.5, 0.9, 0.99)
    for seed in cfg.mc_seeds:
        rng = np.random.default_rng(int(seed))
        Es = mc.random_perturbations(rng, n, e_norm)
        norm_dev = float(np.max(np.abs(np.linalg.norm(Es, 2, axis=(1, 2)) / e_norm - 1.0)))
        row = {"seed": int(seed), "n": n, "E_norm": e_norm, "E_norm_max_rel_dev": norm_dev,
               "worst_shift": {}, "mean_shift": {}, "quantiles": {}, "fraction_complex": {},
               "worst_over_bound": {}}
        err_np = 0.0
        m = min(n, 200)
        for k in ORDER:
            shifts = mc.eig_shift_batch(mats[k], Es)
            mu = mc.eig2x2_batch(mats[k][None] + Es)
            row["worst_shift"][k] = float(np.max(shifts))
            row["mean_shift"][k] = float(np.mean(shifts))
            row["quantiles"][k] = {f"q{int(q * 100)}": float(np.quantile(shifts, q)) for q in quantiles}
            row["fraction_complex"][k] = float(np.mean(mu[:, 0].imag != 0.0))
            row["worst_over_bound"][k] = row["worst_shift"][k] / bound[k]
            w = np.linalg.eigvals(mats[k])
            mu_np = np.linalg.eigvals(mats[k][None] + Es[:m])
            ref = np.abs(mu_np[:, :, None] - w[None, None, :]).min(axis=2).max(axis=1)
            err_np = max(err_np, float(np.max(np.abs(shifts[:m] - ref))))
        row["ratio_idle_over_cruise"] = row["worst_shift"]["idle"] / row["worst_shift"]["cruise"]
        per_seed.append(row)

        tag = f"Monte Carlo seed {seed}, n = {n}, ||E||_2 = {e_norm:g}"
        chk.add(f"{tag}: every perturbation has ||E||_2 == E_norm (numpy 2-norm, rel. dev.)",
                norm_dev <= 1e-12, norm_dev, 1e-12)
        chk.add(f"{tag}: own eig_shift_batch vs numpy.linalg.eigvals ({m} samples x 3 points)",
                err_np <= TOL_LIB, err_np, TOL_LIB)
        ok_bf = all(row["worst_shift"][k] <= bound[k] * (1.0 + 1e-9) for k in ORDER)
        chk.add(f"{tag}: Bauer-Fike holds, worst shift <= kappa(V) ||E||_2 at idle/mid/cruise "
                "(largest shift/bound)", ok_bf, max(row["worst_over_bound"].values()), 1.0)
        if brief_applies:
            for k in ("idle", "cruise"):
                rel = abs(row["worst_shift"][k] - BRIEF_MC[k]) / BRIEF_MC[k]
                chk.add(f"{tag}: worst eigenvalue shift at {k} {row['worst_shift'][k]:.5f} vs brief "
                        f"{BRIEF_MC[k]} (relative diff)", rel <= TOL_MC_REL, rel, TOL_MC_REL)
            r = row["ratio_idle_over_cruise"]
            chk.add(f"{tag}: idle/cruise worst-shift ratio {r:.2f} (brief: about {BRIEF_MC['ratio']:g}x)",
                    MC_RATIO_RANGE[0] <= r <= MC_RATIO_RANGE[1], r,
                    f"between {MC_RATIO_RANGE[0]:g} and {MC_RATIO_RANGE[1]:g}")
    if not brief_applies:
        chk.add("Monte Carlo comparison with the brief skipped: cfg.E_norm differs from the brief's 1e-3",
                True, e_norm, BRIEF_MC["E_norm"])
    worst_all = {k: max(r["worst_shift"][k] for r in per_seed) for k in ORDER}
    spread = {k: (max(r["worst_shift"][k] for r in per_seed) - min(r["worst_shift"][k] for r in per_seed))
              / worst_all[k] for k in ORDER}
    return {
        "label": f"Monte Carlo eigenvalue shift - {LABEL}",
        "E_norm": e_norm, "n": n, "seeds": [int(s) for s in cfg.mc_seeds],
        "shift_definition": "max_i min_j |mu_i(A+E) - lambda_j(A)|",
        "kappaV": kappa, "bauer_fike_bound": bound, "bound_formula": "kappa2(V) * ||E||_2",
        "per_seed": per_seed, "worst_over_seeds": worst_all, "seed_spread_rel": spread,
        "ratio_idle_over_cruise_over_seeds": worst_all["idle"] / worst_all["cruise"],
        "bound_ratio_idle_over_cruise": bound["idle"] / bound["cruise"],
        "brief": {"E_norm": BRIEF_MC["E_norm"], "worst_shift_idle": BRIEF_MC["idle"],
                  "worst_shift_cruise": BRIEF_MC["cruise"], "ratio": BRIEF_MC["ratio"],
                  "compared": brief_applies},
    }


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------
def run(cfg: StudyConfig, out_path: str | Path | None = None) -> dict:
    """Run the published-example anchor and write `outputs/data/anchors.json`.

    `out_path` (optional, used by the tests) redirects the JSON file. Returns the
    same dictionary that is written (plain JSON types). The file contains no time
    stamps or paths, so two runs give byte-identical files.

    Course topic: Eigenvalues
    """
    print(f"{MODULE}: {DISCLAIMER}")
    chk = Checks()
    rounding_only: list = []
    points = {key: analyse_point(key, chk, rounding_only) for key in ORDER}
    result = {
        "module": MODULE,
        "label": LABEL,
        "engine": LABEL,
        "disclaimer": DISCLAIMER,
        "source": SOURCE,
        "matrix_provenance": "matrices as printed in the source, transcribed from the project brief",
        "config": {"fast": bool(cfg.fast), "n_mc": int(cfg.n_mc), "mc_seeds": [int(s) for s in cfg.mc_seeds],
                   "E_norm": float(cfg.E_norm)},
        "tolerances": {"eig_tau_abs": TOL_ABS, "tau_printed_rounding": TOL_TAU_PRINTED,
                       "kappa_rel": TOL_KAPPA_REL, "gauss_seidel_abs": TOL_GS, "mc_rel": TOL_MC_REL,
                       "mc_ratio_range": list(MC_RATIO_RANGE), "own_vs_library": TOL_LIB},
        "order": list(ORDER),
        "points": points,
        "gauss_seidel": analyse_gauss_seidel(chk),
        "flip": analyse_flip(chk),
        "monte_carlo": analyse_monte_carlo(cfg, chk),
        "matched_to_printed_rounding_only": rounding_only,
        "checks": chk.items,
    }
    result = _plain(result)
    path = Path(out_path) if out_path is not None else DATA / f"{MODULE}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(result, fh, indent=1, sort_keys=True)
        fh.write("\n")
    n_fail = sum(1 for c in chk.items if not c["pass"])
    print(f"{MODULE}: {len(chk.items) - n_fail} PASS, {n_fail} FAIL -> {path.name}")
    return result


def main(argv: list[str] | None = None) -> int:
    """Command-line entry: `python -m eigenmap.anchors [--fast]`; exit code 1 on any FAIL.

    Course topic: (none - project plumbing)
    """
    result = run(from_argv(argv))
    return 1 if any(not c["pass"] for c in result["checks"]) else 0


if __name__ == "__main__":
    sys.exit(main())
