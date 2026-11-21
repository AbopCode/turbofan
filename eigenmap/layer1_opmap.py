"""LAYER 1 - the operating-point map: eigenvalues of A from ground idle to max power.

At every steady operating point (x0, u0) the engine is linearized,
    d(dx)/dt = A dx + B du,   A = dF/dx,   B = dF/du,
and the 2x2 matrix A is analysed with the project's own formulas:
    characteristic equation  lam^2 - tr(A) lam + det(A) = 0,
    stability  tr(A) < 0 and det(A) > 0,   time constants tau = -1/Re(lam),
    eigenvectors (mode shapes) and the condition number kappa(V) of the eigenvector matrix.

Tests hypothesis H1: "the engine responds more slowly at low power".

Run: python -m eigenmap.layer1_opmap [--fast]
Course topic: Jacobians / Taylor series (linearization); Eigenvalues; Stability
"""
from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr

from . import config
from . import mathcore as mc
from .common import Checks, eig_summary, finish, write_json
from .linearize import normalized_A, reference_AB
from .opoints import standard_points
from .steady_state import solve_steady


def point_record(engine, p) -> dict:
    """Linear model and eigen-analysis of one operating point (raw and normalized state forms)."""
    ab = reference_AB(engine, p.sp, tol=1e-13, jac=p.jac)
    A, B, An = ab["A"], ab["B"], ab["A_n"]
    raw, nrm = eig_summary(A), eig_summary(An)
    lam_np = np.linalg.eigvals(A)
    own = np.array(raw["lam"])
    err_np = max(min(abs(l - m) for m in lam_np) for l in own) / max(abs(lam_np))
    sp = p.sp
    return {
        "id": p.id, "label": p.label, "group": p.group, "power": p.power,
        "x0": sp.x0, "z0": sp.z0, "Wf": sp.ctrl.Wf, "VAFN": sp.ctrl.VAFN, "VBV": sp.ctrl.VBV,
        "Fnet": sp.out["Fnet"], "Nc_fan": sp.out["Nc_fan"], "N_fan": sp.out["N_fan"], "Ps3": sp.out["Ps3"],
        "T4": sp.out["T4"], "TSFC": sp.out["TSFC"], "W_in": sp.out["W_in"], "BPR": sp.out["BPR"],
        "SM_fan": sp.out["SM_fan"], "SM_lpc": sp.out["SM_lpc"], "SM_hpc": sp.out["SM_hpc"],
        "env": [sp.env.alt_ft, sp.env.MN, sp.env.dT],
        "steady_resnorm": sp.info["resnorm"], "F0_rel": float(np.max(np.abs(ab["F0"]) / sp.x0)),
        "A": A, "B": B.reshape(-1), "A_n": An, "B_n": (B.reshape(-1) / sp.x0) * sp.ctrl.Wf,
        "ref_rel_diff_A": ab["rel_diff_A"], "ref_rel_diff_B": ab["rel_diff_B"], "ref_method": ab["method"],
        "lam": raw["lam"], "tau_slow": raw["tau_slow"], "tau_fast": raw["tau_fast"],
        "tr": raw["tr"], "det": raw["det"], "disc": raw["disc"], "complex_pair": raw["complex_pair"],
        "stable": raw["stable"], "kappaV_raw": raw["kappaV"], "V_raw": raw["V"],
        "kappaV": nrm["kappaV"], "V_n": nrm["V"], "disc_n": nrm["disc"],
        "own_vs_numpy_eig": float(err_np),
        "gs_rho": mc.gs_spectral_radius_2x2(A), "diag_dominant": bool(mc.is_diag_dominant(A)),
    }


def h1_verdict(recs: list) -> dict:
    """H1 'slower at low power', decided by an explicit rule from the numbers."""
    groups = {}
    for g in ("SLS", "CRUISE"):
        rs = [r for r in recs if r["group"] == g]
        pw = [r["power"] for r in rs]
        ts = [r["tau_slow"] for r in rs]
        tf = [r["tau_fast"] for r in rs]
        rho = float(spearmanr(pw, ts).statistic)
        groups[g] = {"spearman_power_vs_tau_slow": rho,
                     "spearman_power_vs_tau_fast": float(spearmanr(pw, tf).statistic),
                     "tau_slow_idle": ts[0], "tau_slow_max": ts[-1], "ratio_idle_over_max": ts[0] / ts[-1],
                     "tau_fast_idle": tf[0], "tau_fast_max": tf[-1]}
    ok = [groups[g]["spearman_power_vs_tau_slow"] < -0.7 and groups[g]["ratio_idle_over_max"] > 1.5 for g in groups]
    bad = [groups[g]["ratio_idle_over_max"] < 1.0 for g in groups]
    verdict = "supported" if all(ok) else ("not supported" if all(bad) else "mixed")
    return {"verdict": verdict, "numbers": groups,
            "rule": "supported if, in BOTH the SLS and the cruise sweep, Spearman rank correlation(power, tau_slow) < -0.7 "
                    "AND tau_slow(idle)/tau_slow(max) > 1.5; not supported if the ratio < 1 in both; otherwise mixed"}


def run(cfg: config.StudyConfig) -> dict:
    ck = Checks("layer1_opmap")
    engine, pts = standard_points(cfg)
    recs = [point_record(engine, p) for p in pts]
    for r in recs:
        ck.add("%s steady residual (max|r|, |dN/dt|/N)" % r["id"], r["steady_resnorm"] < 1e-8, r["steady_resnorm"], 1e-8)
        ck.add("%s reference Jacobian self-consistent (two step sizes, rel)" % r["id"], max(r["ref_rel_diff_A"], r["ref_rel_diff_B"]) < 1e-6,
               max(r["ref_rel_diff_A"], r["ref_rel_diff_B"]), 1e-6)
        ck.add("%s own 2x2 eigenvalues vs numpy (rel)" % r["id"], r["own_vs_numpy_eig"] < 1e-8, r["own_vs_numpy_eig"], 1e-8)
        A = np.array(r["A"])
        ck.add("%s diagonal of A negative and stable (max of a11, a22, Re lam)" % r["id"],
               A[0, 0] < 0 and A[1, 1] < 0 and r["stable"], float(max(A[0, 0], A[1, 1], max(np.real(r["lam"])))), 0.0)
    n_sls = sum(1 for r in recs if r["group"] == "SLS")
    ck.add("number of SLS power settings (ground idle .. max)", n_sls >= (4 if cfg.fast else 8), n_sls, 4 if cfg.fast else 8)
    op = [r for r in recs if r["group"] != "ANCHOR"]
    h1 = h1_verdict(op)
    flips = [r["id"] for r in op if r["complex_pair"]]
    # smallest change of a12 (normalized form) that turns a real pair into a complex pair or back: disc = 0
    flip_margin = []
    for r in op:
        An = np.array(r["A_n"])
        a12_crit = -(An[0, 0] - An[1, 1]) ** 2 / (4.0 * An[1, 0]) if An[1, 0] != 0 else float("inf")
        flip_margin.append({"id": r["id"], "a12_n": An[0, 1], "a12_n_at_disc0": a12_crit,
                            "rel_change_needed": abs(a12_crit - An[0, 1]) / abs(An[0, 1])})
    # ---- multiple equilibria at idle: same corrected fan speed, same VAFN/VBV, solved from NASA's default start ----
    branches = {}
    for p in pts:
        if p.id not in ("SLS_00", "CRUISE_00"):
            continue
        sp = p.sp
        alt = solve_steady(engine, sp.env, target=("Nc_fan", float(sp.out["Nc_fan"])), VAFN=sp.ctrl.VAFN, VBV=sp.ctrl.VBV, tol=1e-12)
        b = {"study_branch": {"Wf": sp.ctrl.Wf, "x0": sp.x0, "z0": sp.z0, "T4": sp.out["T4"]}, "other_found": False}
        if alt.info["converged"] and abs(alt.ctrl.Wf / sp.ctrl.Wf - 1.0) > 1e-4:
            ab = reference_AB(engine, alt, tol=1e-13)
            e2 = eig_summary(ab["A"])
            en = eig_summary(ab["A_n"])
            b.update({"other_found": True, "other": {"Wf": alt.ctrl.Wf, "x0": alt.x0, "z0": alt.z0, "T4": alt.out["T4"], "A": ab["A"],
                                                     "lam": e2["lam"], "tau_slow": e2["tau_slow"], "tau_fast": e2["tau_fast"],
                                                     "kappaV": en["kappaV"], "stable": e2["stable"], "resnorm": alt.info["resnorm"]}})
        branches[p.id] = b
    # does H1 change if the other idle branch were used?
    h1_alt = {}
    for g, pid in (("SLS", "SLS_00"), ("CRUISE", "CRUISE_00")):
        if branches.get(pid, {}).get("other_found"):
            n = h1["numbers"][g]
            h1_alt[g] = {"tau_slow_idle_other_branch": branches[pid]["other"]["tau_slow"],
                         "ratio_idle_over_max_other_branch": branches[pid]["other"]["tau_slow"] / n["tau_slow_max"]}
    h1["idle_branch_sensitivity"] = h1_alt
    res = {"points": recs, "verdicts": {"H1_slower_at_low_power": h1}, "complex_pair_points": flips, "idle_branches": branches,
           "idle_branch_note": "At low power the AGTF30 model (faithful port too, checked) has more than one steady state for the same "
                               "corrected fan speed, VAFN and VBV. The study uses the branch reached by continuation from max power along "
                               "NASA's VAFN/VBV schedules (what a decelerating engine follows); 'other' is the steady state found from "
                               "NASA's default initial guess. Each operating point is reproducible from its tabulated z0 (warm start). "
                               "The flow-matching equations themselves can also have more than one root at a fixed (x, u): from NASA's "
                               "default start the inner solver can converge to a spurious root (e.g. LPT pressure ratio ~150 at SLS idle) that "
                               "is not the engine's state; all F evaluations of the study therefore start from the operating point's own z0.",
           "flip_margin": flip_margin, "kappa_form_note": "kappaV = normalized states (fraction of operating speed), "
           "A_n = D^-1 A D, D = diag(x0); kappaV_raw = rpm states. Eigenvalues are the same in both forms.",
           "held_constant": ["VAFN", "VBV", "ambient", "theta"], "fast": cfg.fast}
    print("H1 slower at low power: %s  (SLS tau_slow %.3f s -> %.3f s; cruise %.3f s -> %.3f s)" % (
        h1["verdict"], h1["numbers"]["SLS"]["tau_slow_idle"], h1["numbers"]["SLS"]["tau_slow_max"],
        h1["numbers"]["CRUISE"]["tau_slow_idle"], h1["numbers"]["CRUISE"]["tau_slow_max"]))
    write_json("layer1_opmap", res, ck)
    return ck


if __name__ == "__main__":
    c = config.from_argv()
    c.ensure_dirs()
    finish(run(c))
