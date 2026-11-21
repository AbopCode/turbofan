"""GATE 1: the port (agent A) against the blind replica (agent A2).

Run:  cd <project root> && python -m tests.gate1_compare

Two engineers ported the NASA AGTF30 engine from the T-MATS source without seeing
each other's code. If both ports agree to 1e-6 on raw residuals, torques and
steady-state points, a porting error would have to be made identically twice.

Prints PASS/FAIL lines; exit code 1 on any FAIL.
"""
from __future__ import annotations

import itertools
import sys

import numpy as np

from eigenmap import engine_model as A
from eigenmap import steady_state as AS
from replica import engine_replica as R
from replica import steady_replica as RS

TOL = 1e-6
SEED = 20261114
fails = 0


def rel(a: float, b: float, floor: float = 1e-12) -> float:
    return abs(a - b) / max(abs(a), abs(b), floor)


def check(name: str, value: float, limit: float) -> None:
    global fails
    ok = value <= limit
    fails += (not ok)
    print("%s gate1: %s (%.3e vs %.1e)" % ("PASS" if ok else "FAIL", name, value, limit))


def best_perm(ra: np.ndarray, rb: np.ndarray):
    """Residual ordering may legitimately differ; find the signed permutation that matches."""
    best = (1e99, None, None)
    for perm in itertools.permutations(range(7)):
        rp = rb[list(perm)]
        for_sign = np.where(np.abs(ra - rp) <= np.abs(ra + rp), 1.0, -1.0)
        d = float(np.max(np.abs(ra - for_sign * rp)))
        if d < best[0]:
            best = (d, perm, for_sign)
    return best


def main() -> int:
    ea, er = A.load_engine(mode="nasa"), R.load_engine()
    rng = np.random.default_rng(SEED)
    nominal = [
        # env (alt, MN, dT), x, z (NASA ICcom order), Wf, VAFN, VBV
        ((35000.0, 0.8, 0.0), [6772.0, 20796.0], [815.48, 4.0955, 1.994, 2.0007, 23.19937, 2.19956, 11.0178], 0.7837, 4775.22, 0.0),
        ((0.0, 0.0, 27.0), [6079.0, 21583.2], [1723.81, 4.17, 2.0325, 1.75, 27.5078, 1.8740, 7.15], 1.3896, 6314.82, 0.0),
    ]
    print("--- random states: raw residuals and torques (no solver involved) ---")
    for k in range(5):
        env_t, x, z, wf, vafn, vbv = nominal[k % 2]
        x = [v * (1 + 0.03 * rng.uniform(-1, 1)) for v in x]
        z = [v * (1 + 0.02 * rng.uniform(-1, 1)) for v in z]
        wf = wf * (1 + 0.05 * rng.uniform(-1, 1))
        vafn = vafn * (1 + 0.03 * rng.uniform(-1, 1))
        vbv = float(rng.uniform(0, 0.3)) if k >= 2 else 0.0
        theta = None if k < 4 else {"hpc_eff": 0.98, "hpt_flow": 0.99, "fan_eff": 0.97, "lpt_eff": 0.985, "lpc_flow": 0.98}
        ra, oa = ea.residuals(z, x, A.Controls(wf, vafn, vbv), A.Env(*env_t), theta)
        rb, ob = er.residuals(z, x, R.Controls(wf, vafn, vbv), R.Env(*env_t), theta)
        ra, rb = np.asarray(ra, float), np.asarray(rb, float)
        d_same = float(np.max(np.abs(ra - rb)))
        if d_same > TOL:
            d, perm, sign = best_perm(ra, rb)
            print("  state %d: residual order/sign differs; best signed permutation %s %s -> max abs diff %.3e" % (k, perm, sign, d))
        else:
            d = d_same
        check("state %d residuals max abs diff (normalized residuals)" % k, d, TOL)
        for key in ("Trq_LP", "Trq_HP", "Ndot_LP", "Ndot_HP", "Fnet", "T4", "Ps3", "W_vbv", "Nc_fan", "TSFC", "OPR", "T45", "SM_hpc", "SM_lpc", "SM_fan"):
            if key in oa and key in ob:
                check("state %d %s rel diff (A=%.8g, A2=%.8g)" % (k, key, oa[key], ob[key]), rel(oa[key], ob[key], 1e-6), TOL)
            else:
                print("  note: key %s missing in %s" % (key, "A" if key not in oa else "A2"))

    print("--- steady-state points (each port's own Newton-Raphson solver) ---")
    cases = [
        ("cruise, Wf fixed", (35000.0, 0.8, 0.0), ("Wf", 0.7837), 4775.22, 0.0),
        ("takeoff hot day, Wf fixed", (0.0, 0.0, 27.0), ("Wf", 1.3896), 6314.82, 0.0),
        ("cruise, Nc_fan target", (35000.0, 0.8, 0.0), ("Nc_fan", 2359.983186), 4775.22, 0.0),
        ("SLS part power, scheduled VAFN/VBV", (0.0, 0.0, 0.0), ("Nc_fan", 1500.0), None, None),
        ("SLS low power, scheduled VAFN/VBV", (0.0, 0.0, 0.0), ("Nc_fan", 1200.0), None, None),
        ("SLS degraded (hpc_eff 0.97), Wf fixed", (0.0, 0.0, 0.0), ("Wf", 1.0), 6300.0, 0.0),
    ]
    for name, env_t, target, vafn, vbv in cases:
        theta = {"hpc_eff": 0.97} if "degraded" in name else None
        try:
            sa = AS.solve_steady(ea, A.Env(*env_t), target=target, theta=theta, VAFN=vafn, VBV=vbv, tol=1e-12)
            sb = RS.solve_steady(er, R.Env(*env_t), target=target, theta=theta, VAFN=vafn, VBV=vbv, tol=1e-12)
        except Exception as exc:  # report, do not hide
            global fails
            fails += 1
            print("FAIL gate1: %s raised %r" % (name, exc))
            continue
        check("%s: N_LP rel diff (A=%.6f, A2=%.6f)" % (name, sa.x0[0], sb.x0[0]), rel(sa.x0[0], sb.x0[0]), TOL)
        check("%s: N_HP rel diff (A=%.6f, A2=%.6f)" % (name, sa.x0[1], sb.x0[1]), rel(sa.x0[1], sb.x0[1]), TOL)
        check("%s: Wf rel diff (A=%.8f)" % (name, sa.ctrl.Wf), rel(sa.ctrl.Wf, sb.ctrl.Wf), TOL)
        check("%s: z max rel diff" % name, max(rel(a, b) for a, b in zip(sa.z0, sb.z0)), TOL)
        for key in ("Fnet", "W_in", "TSFC", "T4", "OPR", "BPR"):
            check("%s: %s rel diff (A=%.8g)" % (name, key, sa.out[key]), rel(sa.out[key], sb.out[key]), TOL)
    print("GATE 1 RESULT: %s (%d failures)" % ("PASS" if fails == 0 else "FAIL", fails))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
