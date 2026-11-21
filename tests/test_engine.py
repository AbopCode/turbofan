"""Acceptance tests for the AGTF30 port (engine_model.py, steady_state.py).

Run:  python -m tests.test_engine        (about 30 s)
Anchors vs NASA Table 1 are checked in `eigenmap.validation`; this file checks solver behaviour,
steady-state residuals across the power range, physical sanity, and speed.
"""
from __future__ import annotations

import sys
import time

import numpy as np
from scipy.integrate import solve_ivp

from eigenmap.engine_model import Controls, Env, load_engine
from eigenmap.steady_state import solve_steady
from eigenmap.linearize import inner_jacobian, reference_AB

fails = 0


def check(name: str, ok: bool, value, limit) -> None:
    global fails
    fails += (not ok)
    print("%s engine: %s (%.3g vs %.3g)" % ("PASS" if ok else "FAIL", name, value, limit), flush=True)


def main() -> int:
    eng = load_engine()                       # study mode (the engine the study uses)
    sls = Env(0.0, 0.0, 0.0)
    cruise = Env(35000.0, 0.8, 0.0)
    # 1. steady residuals across the SLS power range, scheduled VAFN/VBV, continuation from max power down
    sp = solve_steady(eng, sls, target=("Nc_fan", 1911.586), tol=1e-12)
    pts = [sp]
    for nc in (1800, 1650, 1500, 1350, 1200, 1050, 900):
        sp = solve_steady(eng, sls, target=("Nc_fan", float(nc)), tol=1e-12, guess=sp)
        pts.append(sp)
    for p in pts:
        check("SLS Nc_fan %.0f steady: max|r| and |dN/dt|/N" % p.out["Nc_fan"], p.info["converged"] and p.info["resnorm"] < 1e-10,
              p.info["resnorm"], 1e-10)
    # 2. the inner solver honours tol, and looser tol -> larger error in F
    base = solve_steady(eng, cruise, target=("Nc_fan", 2359.983186), tol=1e-13)
    x = base.x0 * np.array([1.01, 0.995])
    J = inner_jacobian(eng, base)
    F_ref = eng.F(x, base.ctrl, cruise, tol=1e-14, z0=base.z0, jac=J)
    prev = 0.0
    errs = []
    for tol in (1e-13, 1e-9, 1e-6, 1e-4):
        z, out, info = eng.solve_flow(x, base.ctrl, cruise, tol=tol, z0=base.z0, jac=J)
        check("solve_flow resnorm < tol=%g" % tol, info["converged"] and info["resnorm"] < tol, info["resnorm"], tol)
        errs.append(float(np.max(np.abs(np.array([out["Ndot_LP"], out["Ndot_HP"]]) - F_ref)) / np.max(np.abs(F_ref))))
    check("error in F grows as tol loosens (1e-13 -> 1e-4)", errs[0] < errs[2] < errs[3] and errs[0] < 1e-9, errs[3], 1.0)
    # 3. physical sanity: diag(A) < 0 at cruise and take-off, and a fuel step up raises both speeds
    for name, env, nc in (("cruise", cruise, 2359.983186), ("take-off", Env(0, 0, 27), 1911.586)):
        p = solve_steady(eng, env, target=("Nc_fan", nc), tol=1e-13)
        ab = reference_AB(eng, p)
        A = ab["A"]
        check("%s dF_LP/dN_LP < 0" % name, A[0, 0] < 0, A[0, 0], 0.0)
        check("%s dF_HP/dN_HP < 0" % name, A[1, 1] < 0, A[1, 1], 0.0)
        print("     %s A = %s, eig = %s" % (name, np.array2string(A, precision=4), np.round(np.linalg.eigvals(A), 4)))
        ctrl = Controls(p.ctrl.Wf * 1.02, p.ctrl.VAFN, p.ctrl.VBV)
        jac = inner_jacobian(eng, p)
        sol = solve_ivp(lambda t, xx: eng.F(xx, ctrl, env, tol=1e-10, z0=p.z0, jac=jac), (0, 1.0), p.x0, rtol=1e-8, atol=1e-6)
        d = sol.y[:, -1] - p.x0
        check("%s 2%% fuel step raises N_LP (rpm after 1 s)" % name, d[0] > 0, d[0], 0.0)
        check("%s 2%% fuel step raises N_HP (rpm after 1 s)" % name, d[1] > 0, d[1], 0.0)
    # 4. timing
    t = time.perf_counter()
    for _ in range(200):
        eng.residuals(base.z0, base.x0, base.ctrl, cruise)
    ms = (time.perf_counter() - t) / 200 * 1e3
    check("one residuals() call [ms]", ms < 2.0, ms, 2.0)
    t = time.perf_counter()
    for k in range(50):
        eng.F(base.x0 * (1 + 1e-4 * k), base.ctrl, cruise, tol=1e-12, z0=base.z0, jac=J)
    ms = (time.perf_counter() - t) / 50 * 1e3
    check("one warm-started F() call at tol 1e-12 [ms]", ms < 20.0, ms, 20.0)
    print("engine: %d failed" % fails)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
