"""Item 6 - step response and new steady state after a small fuel step, at every operating point.

Linear model d(dx)/dt = A dx + B du with a constant fuel step du (1 % of W_f):
    step response   dx(t) = A^-1 (e^{At} - I) B du            (own diagonalization A = P D P^-1)
    cross-checked against scipy.linalg.expm and scipy.integrate.solve_ivp on the same linear ODE;
    new steady state dx_ss = -A^-1 B du  (own 2x2 inverse) and by own Gauss-Seidel iteration,
    with the Gauss-Seidel spectral radius |a12 a21 / (a11 a22)| and the diagonal-dominance test.

Run: python -m eigenmap.steps [--fast]
Course topic: Linear ODE systems; Diagonalization; Gauss-Seidel
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from . import config
from . import mathcore as mc
from .common import Checks, finish, write_json
from .layer1_opmap import point_record
from .opoints import standard_points

STEP_FRAC = 0.01        # fuel step = 1 % of the operating-point fuel flow


def analyse(r: dict) -> dict:
    """Step response and Gauss-Seidel steady state for one operating point record."""
    A = np.array(r["A"])
    B = np.array(r["B"]).reshape(2, 1)
    du = STEP_FRAC * r["Wf"]
    t = np.linspace(0.0, 5.0 * r["tau_slow"], 101)
    own = mc.step_response(A, B, du, t)                                   # own diagonalization
    Ainv = np.linalg.inv(A)
    lib = np.array([(Ainv @ (expm(A * ti) - np.eye(2)) @ B).ravel() * du for ti in t])
    sol = solve_ivp(lambda _t, y: A @ y + B.ravel() * du, (0.0, t[-1]), [0.0, 0.0], method="DOP853",
                    t_eval=t, rtol=1e-12, atol=1e-12 * max(1.0, np.max(np.abs(lib))))
    ivp = sol.y.T
    scale = np.max(np.abs(lib))
    b = -(B.ravel() * du)                                                 # A dx_ss = -B du
    x_direct = mc.steady_gain(A, B, du)
    x_np = np.linalg.solve(A, b)
    x_gs, hist, info = mc.gauss_seidel(A, b, tol=1e-12, max_iter=500)
    errs = [float(np.max(np.abs(h - x_np)) / np.max(np.abs(x_np))) for h in hist]
    return {
        "id": r["id"], "du": du, "t": t, "dx_own": own, "dx_expm": lib, "dx_ivp": ivp,
        "own_vs_expm": float(np.max(np.abs(own - lib)) / scale), "own_vs_ivp": float(np.max(np.abs(own - ivp)) / scale),
        "dx_ss_direct": x_direct, "dx_ss_numpy": x_np, "dx_ss_gs": x_gs, "gs_history": hist, "gs_rel_err": errs,
        "gs_iters": info.get("iterations", len(hist) - 1), "gs_converged": bool(info.get("converged", True)),
        "gs_rho": mc.gs_spectral_radius_2x2(A), "diag_dominant": bool(mc.is_diag_dominant(A)),
        "diag_dominant_normalized": bool(mc.is_diag_dominant(np.array(r["A_n"]))),
        "gs_vs_direct": float(np.max(np.abs(x_gs - x_direct)) / np.max(np.abs(x_direct))),
        "direct_vs_numpy": float(np.max(np.abs(x_direct - x_np)) / np.max(np.abs(x_np))),
        "dx_ss_pct": (np.asarray(x_direct) / np.array(r["x0"]) * 100.0),
        "final_vs_ss": float(np.max(np.abs(own[-1] - x_direct)) / np.max(np.abs(x_direct))),
    }


def run(cfg: config.StudyConfig) -> Checks:
    ck = Checks("steps")
    engine, pts = standard_points(cfg)
    out = []
    for p in pts:
        r = point_record(engine, p)
        s = analyse(r)
        out.append(s)
        ck.add("%s step response: own diagonalization vs expm (rel)" % s["id"], s["own_vs_expm"] < 1e-8, s["own_vs_expm"], 1e-8)
        ck.add("%s step response: own vs solve_ivp (rel)" % s["id"], s["own_vs_ivp"] < 1e-8, s["own_vs_ivp"], 1e-8)
        ck.add("%s Gauss-Seidel converged, matches direct -A^-1 B du (rel)" % s["id"], s["gs_converged"] and s["gs_vs_direct"] < 1e-8,
               s["gs_vs_direct"], 1e-8)
        ck.add("%s own 2x2 inverse vs numpy.linalg.solve (rel)" % s["id"], s["direct_vs_numpy"] < 1e-8, s["direct_vs_numpy"], 1e-8)
        ck.add("%s Gauss-Seidel spectral radius < 1" % s["id"], s["gs_rho"] < 1.0, s["gs_rho"], 1.0)
    res = {"step_frac": STEP_FRAC, "points": out,
           "note": "Gauss-Seidel solves A dx_ss = -B du with A in rpm units; the spectral radius |a12 a21/(a11 a22)| is the "
                   "same for the normalized form (diagonal similarity), diagonal dominance is not."}
    write_json("steps", res, ck)
    return ck


if __name__ == "__main__":
    c = config.from_argv()
    c.ensure_dirs()
    finish(run(c))
