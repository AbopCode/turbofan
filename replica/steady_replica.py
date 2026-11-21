"""Steady-state solver of the A2 replica: Newton-Raphson on [z(7), N_LP, N_HP (, Wf)].

Equations: 7 flow residuals, 2 torque balances (Ndot_i / N_i = 0) and, when the target is
not fuel flow, (out[key] - value)/|value| = 0. This mirrors AGTF30SysSS.mdl, whose
dependents are [7 engine errors, Ndot(2), Nc1 error, VBV error, VAFN error].
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .engine_replica import NASA_IC, Controls, Engine, Env, _NUMERIC_ERRORS


@dataclass
class SteadyPoint:
    """A converged (or best-effort) steady operating point."""
    env: Env
    ctrl: Controls
    x0: np.ndarray
    z0: np.ndarray
    out: dict
    theta: dict
    info: dict = field(default_factory=dict)


class _Problem:
    """Residual function G(u) of one steady-state problem; u = [z(7), N_LP, N_HP, Wf]."""

    def __init__(self, engine: Engine, env: Env, key: str, value: float, theta, VAFN, VBV):
        self.e, self.env, self.key, self.value = engine, env, key, float(value)
        self.theta, self.VAFN, self.VBV = theta, VAFN, VBV
        self.n = 9 if key == "Wf" else 10
        self.Tt2 = engine.ambient(env)[1]

    def controls(self, u: np.ndarray) -> Controls:
        """Controls at iterate u; schedules use this iterate's own fan corrected speed."""
        Wf = self.value if self.key == "Wf" else float(u[9])
        Nc = (u[7] / self.e.GR) / math.sqrt(self.Tt2 / self.e.TSTD)
        vafn = self.e.schedule_VAFN(Nc, self.env.MN) if self.VAFN is None else float(self.VAFN)
        vbv = self.e.schedule_VBV(Nc, self.env.MN) if self.VBV is None else float(self.VBV)
        return Controls(Wf, vafn, vbv)

    def G(self, u: np.ndarray) -> tuple[np.ndarray, dict, Controls]:
        ctrl = self.controls(u)
        r, out = self.e.residuals(u[:7], u[7:9], ctrl, self.env, self.theta)
        g = np.empty(self.n)
        g[:7] = r
        g[7] = out["Ndot_LP"] / u[7]
        g[8] = out["Ndot_HP"] / u[8]
        if self.n == 10:
            g[9] = (out[self.key] - self.value) / max(abs(self.value), 1e-12)
        return g, out, ctrl


def _newton(prob, u0: np.ndarray, tol: float, max_iter: int = 60):
    """Damped Newton with a fresh forward-difference Jacobian every iteration."""
    u = np.array(u0[: prob.n], dtype=float)
    info = dict(converged=False, iters=0, resnorm=float("inf"), nfev=0)
    try:
        g, out, ctrl = prob.G(u)
    except _NUMERIC_ERRORS:
        return u, {}, None, info
    res = float(np.max(np.abs(g)))
    for it in range(max_iter + 1):
        info.update(iters=it, resnorm=res)
        if res < tol:
            info["converged"] = True
            break
        if it == max_iter:
            break
        Jm = np.empty((prob.n, prob.n))
        try:
            for j in range(prob.n):
                h = 1e-6 * max(abs(u[j]), 1e-3)
                up = u.copy()
                up[j] += h
                Jm[:, j] = (prob.G(up)[0] - g) / h
            du = np.linalg.solve(Jm, -g)
        except (np.linalg.LinAlgError,) + _NUMERIC_ERRORS:
            break
        info["nfev"] += prob.n
        lam, ok = 1.0, False
        while lam >= 1.0 / 1024:
            un = u + lam * du
            try:
                gn, on, cn = prob.G(un)
                resn = float(np.max(np.abs(gn)))
                if math.isfinite(resn) and resn < res:
                    ok = True
                    break
            except _NUMERIC_ERRORS:
                pass
            lam *= 0.5
        if not ok:
            break
        u, g, out, ctrl, res = un, gn, on, cn, resn
    return u, out, ctrl, info


class _Blend:
    """Homotopy between a known start point (s = 0) and the requested problem (s = 1)."""

    def __init__(self, prob: _Problem, start: dict, s: float):
        self.p, self.st, self.s, self.n, self.e = prob, start, s, 10, prob.e
        e0, e1 = start["env"], prob.env
        self.env = Env(e0.alt_ft + s * (e1.alt_ft - e0.alt_ft), e0.MN + s * (e1.MN - e0.MN),
                       e0.dT + s * (e1.dT - e0.dT))
        self.Tt2 = self.e.ambient(self.env)[1]

    def G(self, u: np.ndarray):
        s, p, st = self.s, self.p, self.st
        Nc = (u[7] / self.e.GR) / math.sqrt(self.Tt2 / self.e.TSTD)
        vafn = self.e.schedule_VAFN(Nc, self.env.MN) if p.VAFN is None else float(p.VAFN)
        vbv = self.e.schedule_VBV(Nc, self.env.MN) if p.VBV is None else float(p.VBV)
        ctrl = Controls(float(u[9]), st["VAFN"] + s * (vafn - st["VAFN"]), st["VBV"] + s * (vbv - st["VBV"]))
        r, out = self.e.residuals(u[:7], u[7:9], ctrl, self.env, p.theta)
        g = np.empty(10)
        g[:7] = r
        g[7] = out["Ndot_LP"] / u[7]
        g[8] = out["Ndot_HP"] / u[8]
        tgt = (ctrl.Wf if p.key == "Wf" else out[p.key]) - p.value
        g[9] = (1 - s) * (u[9] - st["Wf"]) / st["Wf"] + s * tgt / max(abs(p.value), 1e-12)
        return g, out, ctrl


def solve_steady(engine: Engine, env: Env, *, target: tuple[str, float], theta: dict | None = None,
                 VAFN: float | None = None, VBV: float | None = None, tol: float = 1e-12,
                 guess: SteadyPoint | None = None) -> SteadyPoint:
    """Steady point with out[target[0]] == target[1] (or Wf fixed when target[0] == 'Wf').

    VAFN / VBV = None -> NASA schedules evaluated at the solution's own Nc_fan and MN.
    Converged when max|r| < tol and |Ndot_i|/N_i < tol (and the relative target error < tol).
    """
    key, value = target
    prob = _Problem(engine, env, key, value, theta, VAFN, VBV)
    if guess is not None:
        start = dict(env=guess.env, VAFN=guess.ctrl.VAFN, VBV=guess.ctrl.VBV, Wf=guess.ctrl.Wf,
                     u=np.concatenate([guess.z0, guess.x0, [guess.ctrl.Wf]]), solved=True)
    else:
        a = NASA_IC["cruise"] if env.alt_ft >= 20000.0 else NASA_IC["takeoff"]
        start = dict(env=Env(a["alt"], a["MN"], a["dT"]), VAFN=a["VAFN"], VBV=a["VBV"], Wf=a["Wf"],
                     u=np.array(list(a["z"]) + [a["N2"], a["N3"], a["Wf"]]), solved=False)
    u, out, ctrl, info = _newton(prob, start["u"], tol)
    if not info["converged"]:
        # continuation: first make the start point an exact steady point, then walk s: 0 -> 1
        us = start["u"].copy()
        ok = True
        if not start["solved"]:
            p0 = _Problem(engine, start["env"], "Wf", start["Wf"], theta, start["VAFN"], start["VBV"])
            u9, _, _, i0 = _newton(p0, us, 1e-9)
            ok = i0["converged"]
            us[:9] = u9
        s, ds = 0.0, 0.25
        while ok and s < 1.0 and ds > 1e-4:
            sn = min(1.0, s + ds)
            ut, _, _, it_ = _newton(_Blend(prob, start, sn), us, 1e-9)
            if it_["converged"]:
                s, us, ds = sn, ut, min(2 * ds, 0.5)
            else:
                ds *= 0.5
        if ok and s >= 1.0:
            u2, out2, ctrl2, info2 = _newton(prob, us, tol)
            if info2["converged"] or info2["resnorm"] < info["resnorm"]:
                u, out, ctrl, info = u2, out2, ctrl2, info2
                info["continuation"] = True
    if ctrl is None:
        ctrl = prob.controls(u)
    return SteadyPoint(env=env, ctrl=ctrl, x0=np.array(u[7:9]), z0=np.array(u[:7]), out=out,
                       theta=dict(theta or {}), info=info)
