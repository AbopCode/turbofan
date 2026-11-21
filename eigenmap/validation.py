"""Validation of the Python port against NASA's published AGTF30 performance (Table 1 of
Chapman & Litt, NASA/TM-2017-219569) and against NASA's own steady-state initial conditions.

No NASA A/B matrices exist in the repository and none were supplied (nasa_truth/nasa_lin.mat
absent), so the linear model itself cannot be validated against NASA numbers; this module
validates the nonlinear steady-state model that the linearization is built on.

Run: python -m eigenmap.validation [--fast]
Course topic: Engine model (not a course topic) - model validation
"""
from __future__ import annotations

import json
import re

import numpy as np

from . import config
from .common import Checks, finish, write_json
from .engine_model import Env, load_engine
from .opoints import ENV_CRUISE, ENV_TAKEOFF_HOT, NC_MAX_CRUISE, NC_MAX_SLS, VAFN_CRUISE_ANCHOR, VAFN_TAKEOFF_ANCHOR
from .steady_state import solve_steady

# NASA/TM-2017-219569, Table 1 (opened 2026-11-14). (value, printed decimals)
TABLE1 = {
    "cruise": {"W_in": (813.0, 0), "TSFC": (0.4637, 4), "BPR": (24.0, 0), "OPR": (55.0, 0), "T4": (3150.0, 0), "Fnet": (6073.0, 0)},
    "takeoff": {"W_in": (1724.0, 0), "TSFC": (0.1751, 4), "BPR": (28.0, 0), "OPR": (37.0, 0), "T4": (3170.0, 0), "Fnet": (28622.0, 0)},
}
# NASA's steady-state initial conditions in setup_InputsSS.m (outputs of NASA's own MATLAB trim)
NASA_IC = {
    "cruise": {"Wf": 0.7837, "N_LP": 6772.0, "N_HP": 20796.0,
               "z": [815.48, 4.0955, 1.994, 2.0007, 23.19937, 2.19956, 11.0178]},
    "takeoff": {"Wf": 1.3896, "N_LP": 6079.0, "N_HP": 21583.2,
                "z": [1723.81, 4.17, 2.0325, 1.75, 27.5078, 1.8740, 7.15]},
}
Z_NAMES = ["W", "HPT_PR", "HPC_Rline", "FAN_Rline", "BPR", "LPC_Rline", "LPT_PR"]


def _anchor(engine, which: str, tol: float):
    if which == "cruise":
        return solve_steady(engine, ENV_CRUISE, target=("Nc_fan", NC_MAX_CRUISE), VAFN=VAFN_CRUISE_ANCHOR, VBV=0.0, tol=tol)
    return solve_steady(engine, ENV_TAKEOFF_HOT, target=("Nc_fan", NC_MAX_SLS), VAFN=VAFN_TAKEOFF_ANCHOR, VBV=0.0, tol=tol)


def _quantity(sp, key: str) -> float:
    """Model value for a Table-1 quantity. OPR = product of the three compressor pressure ratios (ASSUMED definition)."""
    if key == "OPR":
        return sp.out["PR_fan"] * sp.out["PR_lpc"] * sp.out["PR_hpc"]
    return float(sp.out[key])


def run(cfg: config.StudyConfig) -> dict:
    ck = Checks("validation")
    res = {"source": "Chapman & Litt, NASA/TM-2017-219569, Table 1; repo setup_InputsSS.m",
           "nasa_lin_available": False, "modes": {}}
    for mode in ("nasa", "study"):
        eng = load_engine(mode=mode)
        tol = 1e-12 if mode == "study" else 1e-7      # faithful mode carries 1e-6-level internal loop noise
        rows = {}
        for which in ("cruise", "takeoff"):
            sp = _anchor(eng, which, tol)
            row = {"converged": sp.info["converged"], "resnorm": sp.info["resnorm"], "table1": {}, "ic": {}}
            for key, (pub, dec) in TABLE1[which].items():
                val = _quantity(sp, key)
                rel = val / pub - 1.0
                lim = 0.01 + (0.5 * 10.0 ** (-dec)) / pub          # 1 % plus half a unit of the printed digit
                row["table1"][key] = {"model": val, "nasa": pub, "rel_err": rel, "limit": lim}
                if mode == "nasa":
                    ck.add("%s %s vs NASA Table 1 (%.5g vs %.5g): |rel err| <= 1%% + half a unit of Table 1's last printed digit" % (which, key, val, pub),
                           abs(rel) <= lim, abs(rel), lim)
            row["OPR_Pt3_over_Pt2"] = float(sp.out["OPR"])
            ic = NASA_IC[which]
            row["ic"] = {"Wf": [sp.ctrl.Wf, ic["Wf"]], "N_LP": [sp.x0[0], ic["N_LP"]], "N_HP": [sp.x0[1], ic["N_HP"]]}
            for name, a, b in zip(Z_NAMES, sp.z0, ic["z"]):
                row["ic"][name] = [float(a), b]
            row["ic_max_rel"] = max(abs(a / b - 1.0) for a, b in row["ic"].values())
            row["x0"], row["Wf"], row["z0"] = sp.x0, sp.ctrl.Wf, sp.z0
            row["out"] = {k: sp.out[k] for k in ("Fnet", "TSFC", "W_in", "BPR", "OPR", "T4", "T45", "Ps3", "Nc_fan",
                                                  "N_fan", "SM_fan", "SM_lpc", "SM_hpc", "PR_fan", "PR_lpc", "PR_hpc")}
            if mode == "nasa":
                ck.add("%s steady state converged (max residual)" % which, sp.info["resnorm"] < 1e-6, sp.info["resnorm"], 1e-6)
                ck.add("%s shaft speeds within 1%% of NASA's own trim (setup_InputsSS.m)" % which,
                       max(abs(row["ic"]["N_LP"][0] / row["ic"]["N_LP"][1] - 1), abs(row["ic"]["N_HP"][0] / row["ic"]["N_HP"][1] - 1)) < 0.01,
                       max(abs(row["ic"]["N_LP"][0] / row["ic"]["N_LP"][1] - 1), abs(row["ic"]["N_HP"][0] / row["ic"]["N_HP"][1] - 1)), 0.01)
            rows[which] = row
        res["modes"][mode] = rows
    # how far the differentiable 'study' engine moves the anchors relative to the faithful port
    dev = 0.0
    for which in ("cruise", "takeoff"):
        for key in TABLE1[which]:
            a = res["modes"]["study"][which]["table1"][key]["model"]
            b = res["modes"]["nasa"][which]["table1"][key]["model"]
            dev = max(dev, abs(a / b - 1.0))
    res["study_vs_nasa_mode_max_rel"] = dev
    ck.add("study-mode engine (C1 maps, tight loops) within 1% of the faithful port at both anchors", dev < 0.01, dev, 0.01)
    worst = max(abs(v["rel_err"]) for m in ("nasa",) for w in ("cruise", "takeoff")
                for v in res["modes"][m][w]["table1"].values())
    res["headline_max_rel_err_nasa_mode"] = worst
    # Gate-1 summary, if the log exists
    g1 = config.ROOT / "handoff" / "gate1_output.txt"
    if g1.exists():
        txt = g1.read_text()
        n_pass = len(re.findall(r"^PASS", txt, re.M))
        n_fail = len(re.findall(r"^FAIL", txt, re.M))
        worst_g1 = max([float(m) for m in re.findall(r"\(([0-9.e+-]+) vs 1\.0e-06\)", txt)] or [float("nan")])
        res["gate1"] = {"passed": n_pass, "failed": n_fail, "worst_rel_disagreement": worst_g1}
        ck.add("Gate 1: port vs blind replica agree (worst rel. disagreement)", n_fail == 0 and worst_g1 < 1e-6, worst_g1, 1e-6)
    write_json("validation", res, ck)
    return res


if __name__ == "__main__":
    c = config.from_argv()
    c.ensure_dirs()
    finish_ck = run(c)
    data = json.loads((config.DATA / "validation.json").read_text())
    print("headline: worst |rel err| vs NASA Table 1 (faithful port) = %.2f %%" % (100 * data["headline_max_rel_err_nasa_mode"]))
    import sys
    sys.exit(1 if any(not c["pass"] for c in data["checks"]) else 0)
