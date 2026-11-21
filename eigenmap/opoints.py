"""The standard set of steady operating points used by every layer of the study.

Groups
------
SLS     sea-level static, standard day: `cfg.n_sls` power settings from ground idle to max power
CRUISE  35,000 ft, Mach 0.8, standard day: `cfg.n_cruise` power settings from flight idle to max cruise
ANCHOR  the two points published by NASA (used for validation of the port)

How a "power setting" is defined (DERIVED from NASA's baseline controller, see ASSUMPTIONS.md):
  * power is set by corrected fan speed Nc_fan (NASA's controller uses it as the thrust surrogate);
  * idle  = the point where burner pressure Ps3 equals NASA's minimum-Ps3 limiter
            (120 psia at sea level, falling linearly to 40 psia at 40,000 ft -> 50 psia at 35,000 ft);
  * max   = Nc_fan of NASA's published max-thrust point (0.81 * s_Nc) for SLS and
            Nc_fan of the published cruise point (1.0 * s_Nc) for the cruise sweep;
  * VAFN and VBV follow NASA's baseline schedules at every point.

Owned by the orchestrator (written at Gate 1 on the frozen engine).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from .config import StudyConfig
from .engine_model import Env, load_engine
from .steady_state import SteadyPoint, solve_steady
from .linearize import inner_jacobian

ENV_SLS = Env(alt_ft=0.0, MN=0.0, dT=0.0)
ENV_CRUISE = Env(alt_ft=35000.0, MN=0.8, dT=0.0)
ENV_TAKEOFF_HOT = Env(alt_ft=0.0, MN=0.0, dT=27.0)

S_NC_FAN = 2359.983186            # MWS.FAN.s_Nc (rpm), from setup_FAN.m
NC_MAX_SLS = 0.810 * S_NC_FAN     # NASA max-thrust point
NC_MAX_CRUISE = 1.0 * S_NC_FAN    # NASA cruise point
VAFN_CRUISE_ANCHOR = 4775.22      # in^2, setup_InputsSS.m
VAFN_TAKEOFF_ANCHOR = 6314.82


def ps3_min_limit(alt_ft: float) -> float:
    """NASA's minimum-Ps3 limiter (psia): MWS.Cntrl.Ps3mn, linear in altitude, 120 -> 40 over 0 -> 40,000 ft."""
    return float(np.interp(alt_ft, [0.0, 40000.0], [120.0, 40.0]))


@dataclass
class OpPoint:
    """One steady operating point plus the inner-solver Jacobian that is re-used at that point."""
    id: str            # e.g. 'SLS_00'
    label: str         # human-readable
    group: str         # 'SLS' | 'CRUISE' | 'ANCHOR'
    sp: SteadyPoint
    power: float       # net thrust as a fraction of the group's max-power thrust
    jac: np.ndarray    # 7x7 dr/dz at the point


def _sweep(engine, env: Env, nc_max: float, n: int, tol: float) -> list[SteadyPoint]:
    """Steady points from max power down to idle, equally spaced in corrected fan speed.

    Walks downward using each converged point as the starting guess for the next
    (continuation), then finds idle as the point where Ps3 = NASA's minimum limit.
    """
    top = solve_steady(engine, env, target=("Nc_fan", nc_max), tol=tol)
    ps3_idle = ps3_min_limit(env.alt_ft)
    # walk down in small steps towards idle so the idle solve starts close
    guess = top
    nc = nc_max
    idle = None
    while idle is None:
        nc_next = nc - 60.0
        try:
            nxt = solve_steady(engine, env, target=("Nc_fan", nc_next), tol=tol, guess=guess)
        except Exception:
            nxt = None
        if nxt is None or not nxt.info.get("converged", False):
            idle = guess                      # lowest point the model converges to
            break
        if nxt.out["Ps3"] <= ps3_idle:
            try:
                idle = solve_steady(engine, env, target=("Ps3", ps3_idle), tol=tol, guess=nxt)
                if not idle.info.get("converged", False):
                    idle = nxt
            except Exception:
                idle = nxt
            break
        guess, nc = nxt, nc_next
    nc_idle = float(idle.out["Nc_fan"])
    targets = np.linspace(nc_max, nc_idle, n)
    pts = [top]
    guess = top
    for t in targets[1:-1]:
        # continuation in sub-steps of <= 60 rpm
        nsub = max(1, int(np.ceil(abs(float(guess.out["Nc_fan"]) - t) / 60.0)))
        for s in np.linspace(float(guess.out["Nc_fan"]), t, nsub + 1)[1:]:
            guess = solve_steady(engine, env, target=("Nc_fan", float(s)), tol=tol, guess=guess)
        pts.append(guess)
    pts.append(idle)
    return pts[::-1]                          # idle first


def _name(group: str, i: int, n: int) -> str:
    if i == 0:
        return "idle"
    if i == n - 1:
        return "max power" if group == "SLS" else "max cruise"
    return "part power %d" % i


@lru_cache(maxsize=4)
def _points_cached(fast: bool) -> tuple:
    cfg = StudyConfig(fast=fast)
    engine = load_engine()
    tol = cfg.tol_tight
    out: list[OpPoint] = []
    for group, env, nc_max, n in (("SLS", ENV_SLS, NC_MAX_SLS, cfg.n_sls),
                                  ("CRUISE", ENV_CRUISE, NC_MAX_CRUISE, cfg.n_cruise)):
        sps = _sweep(engine, env, nc_max, n, tol)
        fmax = float(sps[-1].out["Fnet"])
        for i, sp in enumerate(sps):
            where = "SLS" if group == "SLS" else "35k ft / M0.8"
            out.append(OpPoint(id="%s_%02d" % (group, i),
                               label="%s %s (Nc_fan %.0f rpm)" % (where, _name(group, i, n), sp.out["Nc_fan"]),
                               group=group, sp=sp, power=float(sp.out["Fnet"]) / fmax,
                               jac=inner_jacobian(engine, sp)))
    a_cr = solve_steady(engine, ENV_CRUISE, target=("Nc_fan", NC_MAX_CRUISE),
                        VAFN=VAFN_CRUISE_ANCHOR, VBV=0.0, tol=tol)
    a_to = solve_steady(engine, ENV_TAKEOFF_HOT, target=("Nc_fan", NC_MAX_SLS),
                        VAFN=VAFN_TAKEOFF_ANCHOR, VBV=0.0, tol=tol)
    out.append(OpPoint("ANCHOR_CRUISE", "NASA cruise anchor (35k ft, M0.8)", "ANCHOR", a_cr, 1.0,
                       inner_jacobian(engine, a_cr)))
    out.append(OpPoint("ANCHOR_TAKEOFF", "NASA max-thrust anchor (SLS, +27 F)", "ANCHOR", a_to, 1.0,
                       inner_jacobian(engine, a_to)))
    return engine, tuple(out)


def standard_points(cfg: StudyConfig) -> tuple:
    """Return (engine, list of OpPoint). Cached per process; deterministic."""
    engine, pts = _points_cached(cfg.fast)
    return engine, list(pts)


def named_points(cfg: StudyConfig) -> dict:
    """The three reference points used by the detailed studies: idle (SLS), cruise (max cruise), max power (SLS)."""
    _, pts = standard_points(cfg)
    sls = [p for p in pts if p.group == "SLS"]
    cr = [p for p in pts if p.group == "CRUISE"]
    return {"idle": sls[0], "cruise_idle": cr[0], "cruise": cr[-1], "max": sls[-1]}
