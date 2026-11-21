"""Steady-state (trim) solver for the AGTF30 port: find x0, u0 with F(x0, u0) = 0.

A steady operating point satisfies 9 equations in 9 unknowns:

    unknowns  y = [z (7 flow-matching variables), N_LP, N_HP]
    equations     r_i(z; x, u) = 0   (7 flow residuals)
                  dN_LP/dt = 0,  dN_HP/dt = 0   (torque balance on each shaft)

for a given fuel flow W_f. If the operating point is specified by something
else (corrected fan speed, burner pressure, thrust ...), W_f becomes a 10th
unknown and "output - target = 0" a 10th equation. The system is solved by
Newton-Raphson with a finite-difference Jacobian, like T-MATS's own
steady-state solver block.

Course topic: Jacobians / partial derivatives (Newton's method uses the Jacobian of the residual vector)
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .engine_model import Controls, Engine, Env, DEFAULT_Z


@dataclass
class SteadyPoint:
    """One converged steady operating point."""
    env: Env
    ctrl: Controls
    x0: np.ndarray            # [N_LP, N_HP] rpm
    z0: np.ndarray            # 7 flow-matching variables
    out: dict                 # engine outputs at the point
    theta: dict = field(default_factory=dict)
    info: dict = field(default_factory=dict)


# NASA's initial guesses for the shaft speeds and fuel flow (setup_InputsSS.m)
_X_CRUISE, _WF_CRUISE = (6772.0, 20796.0), 0.7837
_X_TAKEOFF, _WF_TAKEOFF = (6079.0, 21583.2), 1.3896


def _default_start(env: Env) -> tuple:
    if env.alt_ft > 17500.0:
        return list(DEFAULT_Z(env)), list(_X_CRUISE), _WF_CRUISE
    return list(DEFAULT_Z(env)), list(_X_TAKEOFF), _WF_TAKEOFF


def _controls(engine: Engine, env: Env, N_LP: float, Wf: float, VAFN, VBV) -> Controls:
    """Controls for the current iterate: scheduled VAFN/VBV follow the iterate's own corrected fan speed."""
    if VAFN is None or VBV is None:
        nc = engine.Nc_fan_of(N_LP, env)
    return Controls(Wf=Wf,
                    VAFN=engine.schedule_VAFN(nc, env.MN) if VAFN is None else float(VAFN),
                    VBV=engine.schedule_VBV(nc, env.MN) if VBV is None else float(VBV))


def _newton(engine: Engine, env: Env, target: tuple, theta, VAFN, VBV, tol: float,
            y: list, max_iter: int = 60):
    """Damped Newton-Raphson on the 9 (or 10) steady-state equations. Returns (y, converged, iters, resnorm)."""
    key, value = target
    free_wf = key != "Wf"
    n = 10 if free_wf else 9

    def resid(yv):
        wf = yv[9] if free_wf else value
        ctrl = _controls(engine, env, yv[7], wf, VAFN, VBV)
        r, nd_lp, nd_hp, out = engine._eval(yv[:7], [yv[7], yv[8]], ctrl, env, theta, free_wf)
        f = list(r) + [nd_lp / yv[7], nd_hp / yv[8]]
        if free_wf:
            f.append((out[key] - value) / abs(value))
        return f

    def ok(yv):
        return (yv[0] > 0 and yv[1] > 1.0 and yv[6] > 1.0 and yv[4] > 0 and yv[7] > 100 and yv[8] > 100
                and (not free_wf or yv[9] > 0))

    f = resid(y)
    res = max(abs(v) for v in f)
    iters = 0
    while res >= tol and iters < max_iter:
        J = np.empty((n, n))
        for j in range(n):
            dy = 1e-7 * max(abs(y[j]), 1e-3)
            yp = list(y)
            yp[j] += dy
            fp = resid(yp)
            for i in range(n):
                J[i, j] = (fp[i] - f[i]) / dy
        try:
            step = np.linalg.solve(J, -np.array(f))
        except np.linalg.LinAlgError:
            break
        lam, accepted = 1.0, False
        for _ in range(10):
            yt = [y[i] + lam * float(step[i]) for i in range(n)]
            if ok(yt):
                ft = resid(yt)
                rest = max(abs(v) for v in ft)
                if rest == rest and (rest < res or rest < 10.0 * tol):
                    accepted = True
                    break
            lam *= 0.5
        iters += 1
        if not accepted:
            break
        y, f, res = yt, ft, rest
    return y, res < tol, iters, res


def solve_steady(engine: Engine, env: Env, *, target: tuple, theta: dict | None = None,
                 VAFN: float | None = None, VBV: float | None = None, tol: float = 1e-12,
                 guess: SteadyPoint | None = None) -> SteadyPoint:
    """Find the steady operating point for `target` = ('Wf', value) or (any engine output name, value).

    VAFN / VBV = None follow NASA's baseline schedules (functions of corrected fan speed and Mach,
    evaluated at the solution itself); floats are held fixed.
    If the direct Newton solve fails, the target is approached in smaller steps from the starting
    point (continuation). Converged means max|r| < tol and |dN/dt|/N < tol (1/s) on both shafts.
    Course topic: Jacobians / partial derivatives (Newton-Raphson for nonlinear systems)
    """
    key, value = target
    free_wf = key != "Wf"
    if guess is not None:
        z, x, wf = list(guess.z0), list(guess.x0), float(guess.ctrl.Wf)
    else:
        z, x, wf = _default_start(env)
    y0 = [float(v) for v in z] + [float(x[0]), float(x[1])] + ([wf] if free_wf else [])
    y, conv, iters, res = _newton(engine, env, target, theta, VAFN, VBV, tol, list(y0))
    if not conv:
        # continuation: first converge at the start point's own value of the target, then walk
        ctrl0 = _controls(engine, env, y0[7], wf, VAFN, VBV)
        if free_wf:
            _r, _a, _b, out0 = engine._eval(y0[:7], y0[7:9], ctrl0, env, theta, True)
            v_start = out0[key]
        else:
            v_start = wf
        yc, okc, it0, _res0 = _newton(engine, env, (key, v_start), theta, VAFN, VBV, tol, list(y0))
        iters += it0
        if okc:
            v_now, step = v_start, value - v_start
            fails = 0
            while abs(value - v_now) > 1e-14 * max(1.0, abs(value)) and fails < 12:
                v_try = value if abs(value - v_now) <= abs(step) else v_now + step
                yt, okt, it1, res1 = _newton(engine, env, (key, v_try), theta, VAFN, VBV, tol, list(yc))
                iters += it1
                if okt:
                    yc, v_now, res = yt, v_try, res1
                    if v_try == value:
                        conv, y = True, yt
                        break
                    step *= 1.5
                else:
                    step *= 0.5
                    fails += 1
    wf_sol = y[9] if free_wf else float(value)
    ctrl = _controls(engine, env, y[7], wf_sol, VAFN, VBV)
    r, nd_lp, nd_hp, out = engine._eval(y[:7], [y[7], y[8]], ctrl, env, theta, True)
    resnorm = max(max(abs(v) for v in r), abs(nd_lp / y[7]), abs(nd_hp / y[8]))
    info = {"converged": bool(conv and resnorm < max(tol, 1e-9)), "iters": iters, "resnorm": float(resnorm),
            "flow_resnorm": float(max(abs(v) for v in r)), "target": [key, float(value)]}
    return SteadyPoint(env=env, ctrl=ctrl, x0=np.array([y[7], y[8]]), z0=np.array(y[:7]), out=out,
                       theta=dict(theta or {}), info=info)
