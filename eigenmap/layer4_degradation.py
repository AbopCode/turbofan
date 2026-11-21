"""LAYER 4 - engine wear: how do the eigenvalues move when components lose efficiency or flow capacity?

Health parameters theta (CONTRACT section 3) multiply the map efficiency (X_eff) and map corrected
flow (X_flow) of the fan, booster (LPC), HPC, HPT and LPT; theta = 0.97 means 3 % degradation.
Each parameter is reduced alone, and all ten together ("combined"), from 0 to 5 %, under two
comparison rules:
    same_fuel : the degraded engine burns the SAME fuel flow W_f as the healthy one;
    same_fan  : the controller restores the SAME corrected fan speed (thrust surrogate), so the
                degraded engine needs MORE fuel (the fuel penalty).
At each case the engine is re-trimmed and re-linearized; we record the eigenvalue and time-constant
shifts and the slope d(lam)/d(theta). The sign is measured, never assumed.

Tests hypothesis H4: "degradation makes the engine respond more slowly (tau_slow increases)".

Run: python -m eigenmap.layer4_degradation [--fast]
Course topic: Eigenvalues (sensitivity to parameters); Jacobians
"""
from __future__ import annotations

import numpy as np

from . import config
from . import mathcore as mc
from .common import Checks, finish, write_json
from .engine_model import HEALTH_KEYS
from .linearize import reference_AB
from .opoints import named_points, standard_points
from .steady_state import solve_steady

CASES = list(HEALTH_KEYS) + ["combined"]


def theta_for(case: str, level: float) -> dict:
    """Health-parameter dict for one case at degradation `level` (fraction)."""
    if case == "combined":
        return {k: 1.0 - level for k in HEALTH_KEYS}
    return {case: 1.0 - level}


def lam_slow(A: np.ndarray) -> complex:
    lam, _ = mc.eig2x2(A)
    return complex(min(lam, key=lambda l: abs(l.real)))


def sweep(engine, p, case: str, rule: str, levels) -> dict:
    """Re-trim and re-linearize along the degradation levels for one case and one comparison rule."""
    sp0 = p.sp
    target = ("Wf", float(sp0.ctrl.Wf)) if rule == "same_fuel" else ("Nc_fan", float(sp0.out["Nc_fan"]))
    rows = []
    guess = sp0
    for lev in levels:
        th = theta_for(case, lev) if lev > 0 else {}
        try:
            sp = solve_steady(engine, sp0.env, target=target, theta=th, tol=1e-12, guess=guess)
            ok = sp.info["converged"]
        except Exception:
            ok = False
        if not ok:
            rows.append({"level": lev, "converged": False})
            continue
        guess = sp
        ab = reference_AB(engine, sp, tol=1e-13)
        lam, _ = mc.eig2x2(ab["A"])
        ts, tf = mc.time_constants(lam)
        rows.append({"level": lev, "converged": True, "resnorm": sp.info["resnorm"], "x0": sp.x0, "Wf": sp.ctrl.Wf,
                     "Fnet": sp.out["Fnet"], "TSFC": sp.out["TSFC"], "T4": sp.out["T4"], "Nc_fan": sp.out["Nc_fan"],
                     "SM_hpc": sp.out["SM_hpc"], "SM_lpc": sp.out["SM_lpc"], "SM_fan": sp.out["SM_fan"],
                     "A": ab["A"], "lam": lam, "lam_slow": lam_slow(ab["A"]), "tau_slow": ts, "tau_fast": tf,
                     "ref_rel_diff": ab["rel_diff_A"]})
    base = rows[0]
    for r in rows:
        if r.get("converged") and base.get("converged"):
            r["dtau_slow_pct"] = (r["tau_slow"] / base["tau_slow"] - 1.0) * 100.0
            r["dtau_fast_pct"] = (r["tau_fast"] / base["tau_fast"] - 1.0) * 100.0
            r["fuel_penalty_pct"] = (r["Wf"] / base["Wf"] - 1.0) * 100.0
            r["thrust_change_pct"] = (r["Fnet"] / base["Fnet"] - 1.0) * 100.0
            r["T4_change_R"] = r["T4"] - base["T4"]
    conv = [r for r in rows if r.get("converged")]
    slope = None
    if len(conv) >= 2:     # least-squares slope of Re(lam_slow) vs theta (theta = 1 - level)
        th = np.array([1.0 - r["level"] for r in conv])
        y = np.array([r["lam_slow"].real for r in conv])
        slope = float(np.polyfit(th, y, 1)[0])
    return {"case": case, "rule": rule, "rows": rows, "dRe_lam_slow_dtheta": slope,
            "all_converged": all(r.get("converged") for r in rows)}


def h4_verdict(results: list) -> dict:
    """H4 'degradation slows the engine', from the combined 5 % case at every point under both rules."""
    signs = {}
    for res in results:
        if res["case"] != "combined":
            continue
        last = [r for r in res["rows"] if r.get("converged")][-1]
        signs["%s/%s" % (res["point"], res["rule"])] = {"level": last["level"], "dtau_slow_pct": last.get("dtau_slow_pct"),
                                                        "dtau_fast_pct": last.get("dtau_fast_pct")}
    vals = [v["dtau_slow_pct"] for v in signs.values() if v["dtau_slow_pct"] is not None]
    verdict = "supported" if vals and all(v > 0 for v in vals) else ("not supported" if vals and all(v < 0 for v in vals) else "mixed")
    return {"verdict": verdict, "numbers": signs,
            "rule": "combined degradation (all 10 health parameters reduced together) at the largest converged level: supported "
                    "if tau_slow increases at every Layer-4 point under BOTH comparison rules; not supported if it decreases "
                    "everywhere; otherwise mixed"}


def run(cfg: config.StudyConfig) -> Checks:
    ck = Checks("layer4_degradation")
    engine, _ = standard_points(cfg)
    named = named_points(cfg)
    results = []
    for key in ("idle", "cruise", "max"):
        p = named[key]
        for rule in ("same_fuel", "same_fan"):
            for case in CASES:
                res = sweep(engine, p, case, rule, cfg.deg_grid)
                res["point"], res["point_name"] = p.id, key
                results.append(res)
                ok = res["all_converged"]
                ck.add("%s %s %s: every level 0..%g%% trimmed and linearized" % (p.id, rule, case, 100 * cfg.deg_max), ok,
                       sum(1 for r in res["rows"] if r.get("converged")), len(res["rows"]))
                conv = [r for r in res["rows"] if r.get("converged")]
                if conv:
                    worst = max(r["resnorm"] for r in conv)
                    stable = all(np.all(np.real(r["lam"]) < 0) for r in conv)
                    ck.add("%s %s %s: steady residuals and stability" % (p.id, rule, case), worst < 1e-8 and stable, worst, 1e-8)
        print("layer4: %s done" % key, flush=True)
    h4 = h4_verdict(results)
    # component ranking: |d tau_slow| at the largest level, per point and rule
    ranking = {}
    for res in results:
        last = [r for r in res["rows"] if r.get("converged")]
        if len(last) > 1:
            ranking.setdefault("%s/%s" % (res["point"], res["rule"]), []).append((res["case"], last[-1]["dtau_slow_pct"]))
    for k in ranking:
        ranking[k] = sorted(ranking[k], key=lambda t: -abs(t[1]))
    print("H4 degradation slows the engine: %s  %s" % (h4["verdict"], {k: round(v["dtau_slow_pct"], 2) for k, v in h4["numbers"].items()}))
    write_json("layer4_degradation", {"results": results, "ranking_by_abs_dtau_slow_pct": ranking, "levels": cfg.deg_grid,
                                      "verdicts": {"H4_degradation_direction": h4}, "fast": cfg.fast,
                                      "rules": {"same_fuel": "W_f held at the healthy value; VAFN/VBV on NASA schedules",
                                                "same_fan": "corrected fan speed held at the healthy value (W_f re-trimmed); VAFN/VBV on schedules"}}, ck)
    return ck


if __name__ == "__main__":
    c = config.from_argv()
    c.ensure_dirs()
    finish(run(c))
