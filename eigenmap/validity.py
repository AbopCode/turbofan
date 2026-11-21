"""How far can the linear model be trusted?  Nonlinear vs linear response to fuel steps.

Nonlinear: dx/dt = F(x, W_f0 (1 + s)) integrated with scipy solve_ivp (every F solves the engine).
Linear   : dx(t) = A^-1 (e^{At} - I) B du with du = s W_f0 (own diagonalization).
Steps s = +1, +3, +10, +25 % at ground idle, max cruise and SLS max power. VAFN and VBV are held at
their steady values in both models (same assumption as the linearization).
Error = max over time and both shafts of |dx_nonlin - dx_lin| / max|dx_nonlin| (per shaft).

Run: python -m eigenmap.validity [--fast]
Course topic: Taylor series (range of validity of a linearization); Linear ODE systems
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from . import config
from . import mathcore as mc
from .common import Checks, finish, write_json
from .engine_model import Controls
from .linearize import reference_AB
from .opoints import named_points, standard_points


def simulate(engine, p, frac: float, t: np.ndarray, tau_fast: float) -> np.ndarray:
    """Nonlinear response dx(t) to a fuel step of `frac`, warm-starting each F from the last flow solution."""
    sp = p.sp
    ctrl = Controls(sp.ctrl.Wf * (1.0 + frac), sp.ctrl.VAFN, sp.ctrl.VBV)
    state = {"z": np.array(sp.z0)}

    def rhs(_t, x):
        z, out, _info = engine.solve_flow(x, ctrl, sp.env, sp.theta, tol=1e-11, z0=state["z"], jac=p.jac)
        state["z"] = z
        return [out["Ndot_LP"], out["Ndot_HP"]]

    # step sizes bounded by the fast time constant: RK45's automatic first step can otherwise jump so far
    # from the operating point that the inner flow solver fails (seen at ground idle)
    sol = solve_ivp(rhs, (0.0, float(t[-1])), sp.x0, t_eval=t, rtol=1e-8, atol=1e-6, method="RK45",
                    first_step=0.01 * tau_fast, max_step=0.5 * tau_fast)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y.T - sp.x0


def run(cfg: config.StudyConfig) -> Checks:
    ck = Checks("validity")
    engine, _ = standard_points(cfg)
    named = named_points(cfg)
    fracs = (0.01, 0.10) if cfg.fast else cfg.step_fracs
    out = []
    for key in ("idle", "cruise", "max"):
        p = named[key]
        ab = reference_AB(engine, p.sp, tol=1e-13, jac=p.jac)
        lam, _ = mc.eig2x2(ab["A"])
        ts, tf = mc.time_constants(lam)
        t = np.linspace(0.0, 5.0 * ts, 61 if cfg.fast else 151)
        rec = {"id": p.id, "name": key, "label": p.label, "x0": p.sp.x0, "Wf": p.sp.ctrl.Wf, "t": t, "steps": []}
        for s in fracs:
            lin = mc.step_response(ab["A"], ab["B"], s * p.sp.ctrl.Wf, t)
            try:
                nl = simulate(engine, p, s, t, tf)
            except Exception as exc:
                rec["steps"].append({"frac": s, "ok": False, "error": str(exc)})
                continue
            per_shaft = [float(np.max(np.abs(nl[:, i] - lin[:, i])) / np.max(np.abs(nl[:, i]))) for i in range(2)]
            rec["steps"].append({"frac": s, "ok": True, "dx_nonlin": nl, "dx_lin": lin, "err_per_shaft": per_shaft,
                                 "err": max(per_shaft), "final_dx_pct": nl[-1] / p.sp.x0 * 100.0})
        good = [st for st in rec["steps"] if st["ok"]]
        rec["valid_up_to_2pct"] = max([st["frac"] for st in good if st["err"] <= 0.02], default=0.0)
        rec["valid_up_to_5pct"] = max([st["frac"] for st in good if st["err"] <= 0.05], default=0.0)
        one = [st for st in good if st["frac"] == 0.01]
        ck.add("%s linear vs nonlinear agree for a 1%% fuel step (max rel error)" % p.id, bool(one) and one[0]["err"] <= 0.02,
               one[0]["err"] if one else float("nan"), 0.02)
        for st in good:
            print("   %s step %+.0f%%: max rel error %.4f" % (p.id, 100 * st["frac"], st["err"]))
        out.append(rec)
    write_json("validity", {"points": out, "fracs": fracs, "fast": cfg.fast,
                            "error_definition": "max_t |dx_nl - dx_lin| / max_t |dx_nl|, worst of the two shafts",
                            "held_constant": ["VAFN", "VBV"]}, ck)
    return ck


if __name__ == "__main__":
    c = config.from_argv()
    c.ensure_dirs()
    finish(run(c))
