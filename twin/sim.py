"""Simulation core of the virtual lab (an EXTENSION, outside the report objectives O1-O5).

* `envelope`, `trim_at` : steady operating points at ANY flight condition (altitude, Mach, temperature offset),
                          power setting and engine health, found by continuation from NASA's two anchor points.
* `LabController`       : NASA's baseline fuel-control STRUCTURE as published (Chapman & Litt, NASA/TM-2017-219569):
                          gain-scheduled PI on corrected fan speed + limiter loops (max N1, N3, Ps3, T45; min Ps3)
                          combined by min/max selection, Wf/Ps3 acceleration/deceleration schedules. Gains and limits are
                          read from NASA's setup_Controller.m. The lab adds user multipliers on the gains and on/off
                          switches for each limiter (for experiments; defaults = NASA's values, all limiters on).
* `LiveSim`             : the live, open-ended simulation behind the "Live engine" tab (Start/Stop, no time limit),
                          with every input changeable while it runs.
* `simulate`            : fixed-length runs for the experiment tabs.

Honest simplifications: the controller is our reading of NASA's structure in velocity (incremental) form, not a
block-by-block port of NASA's Simulink controller; no sensor or actuator dynamics apart from simple rate limits;
explicit Euler time step of 0.02 s for the shaft speeds (the flow-matching equations are solved at every step).
"Max power" away from NASA's two anchor points is ASSUMED: corrected fan speed interpolated linearly in altitude
between 0.81 x s_Nc (sea level) and 1.0 x s_Nc (35,000 ft).

Course topic: (extension) Linear ODE systems / stability / Newton's method
"""
from __future__ import annotations

import math
import threading
import time
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np

from eigenmap import mathcore as mc
from eigenmap.engine_model import Controls, Env, HEALTH_KEYS, Map2D, load_engine
from eigenmap.linearize import inner_jacobian, normalized_A, reference_AB
from eigenmap.opoints import ps3_min_limit
from eigenmap.steady_state import SteadyPoint, solve_steady

DT = 0.02                     # s, integration step of the shaft speeds
S_NC = 2359.983186            # MWS.FAN.s_Nc
ANCHORS = {"SLS": (0.0, 0.0, 0.0, 0.810), "CRUISE": (35000.0, 0.8, 0.0, 1.0)}   # alt, MN, dT, Nc_max / s_Nc
LIMITERS = ("N1max", "N3max", "Ps3max", "T45max", "Ps3min", "RUmax", "RUmin")
LIMITER_LABELS = {"N1max": "max fan speed N1", "N3max": "max core speed N3", "Ps3max": "max burner pressure Ps3",
                  "T45max": "max LPT inlet temp T45", "Ps3min": "min burner pressure Ps3 (idle)",
                  "RUmax": "Wf/Ps3 acceleration limit", "RUmin": "Wf/Ps3 deceleration limit"}


def _v(x) -> list:
    return [float(t) for t in np.asarray(x, dtype=float).reshape(-1)]


def health_dict(level_pct: float = 0.0, per_component: dict | None = None) -> dict:
    """theta from a combined degradation level [%] and/or per-component levels {key: %}; 1.0 = healthy."""
    th = {}
    for k in HEALTH_KEYS:
        lv = (per_component or {}).get(k, 0.0) + (level_pct or 0.0)
        if lv:
            th[k] = 1.0 - lv / 100.0
    return th


def nc_max(env: Env) -> float:
    """Max-power corrected fan speed at this flight condition (ASSUMED interpolation between NASA's two anchors)."""
    return S_NC * float(np.interp(env.alt_ft, [0.0, 35000.0], [0.810, 1.0]))


def _key(env: Env, theta: dict) -> tuple:
    return (round(env.alt_ft, 1), round(env.MN, 4), round(env.dT, 2), tuple(sorted((k, round(v, 6)) for k, v in theta.items())))


@lru_cache(maxsize=64)
def _top_cached(key: tuple) -> SteadyPoint:
    alt, MN, dT, th = key
    env, theta = Env(alt, MN, dT), dict(th)
    eng = load_engine()
    # continuation from the nearest NASA anchor in (altitude, Mach, temperature)
    name = "SLS" if alt < 17500.0 else "CRUISE"
    a_alt, a_mn, a_dt, a_nc = ANCHORS[name]
    sp = solve_steady(eng, Env(a_alt, a_mn, a_dt), target=("Nc_fan", a_nc * S_NC), tol=1e-11)
    dist = max(abs(alt - a_alt) / 2500.0, abs(MN - a_mn) / 0.05, abs(dT - a_dt) / 5.0)
    n = max(1, int(math.ceil(dist)))
    for k in range(1, n + 1):
        f = k / n
        e = Env(a_alt + f * (alt - a_alt), a_mn + f * (MN - a_mn), a_dt + f * (dT - a_dt))
        sp = solve_steady(eng, e, target=("Nc_fan", nc_max(e)), tol=1e-11, guess=sp)
        if not sp.info["converged"]:
            raise RuntimeError("no steady state found at alt %.0f ft, Mach %.2f, dT %+.0f F (continuation step %d/%d)"
                               % (e.alt_ft, e.MN, e.dT, k, n))
    if theta:   # degrade at constant fan speed, in steps
        for f in (0.25, 0.5, 0.75, 1.0):
            th_f = {k: 1.0 - f * (1.0 - v) for k, v in theta.items()}
            sp = solve_steady(eng, env, target=("Nc_fan", nc_max(env)), theta=th_f, tol=1e-11, guess=sp)
        if not sp.info["converged"]:
            raise RuntimeError("no steady state with this degradation at max power")
    return sp


@lru_cache(maxsize=64)
def _idle_cached(key: tuple) -> SteadyPoint:
    alt, MN, dT, th = key
    env, theta = Env(alt, MN, dT), dict(th)
    eng = load_engine()
    sp = _top_cached(key)
    lim = ps3_min_limit(alt)
    nc = float(sp.out["Nc_fan"])
    while sp.out["Ps3"] > lim:
        nxt = solve_steady(eng, env, target=("Nc_fan", nc - 60.0), theta=theta, tol=1e-11, guess=sp)
        if not nxt.info["converged"]:
            break
        if nxt.out["Ps3"] <= lim:
            idle = solve_steady(eng, env, target=("Ps3", lim), theta=theta, tol=1e-11, guess=nxt)
            return idle if idle.info["converged"] else nxt
        sp, nc = nxt, nc - 60.0
    return sp


def envelope(env: Env, theta: dict | None = None) -> tuple:
    """(idle point, max point) at this flight condition and health. Idle = NASA's minimum-Ps3 limit."""
    k = _key(env, theta or {})
    return _idle_cached(k), _top_cached(k)


def nc_of_power(env: Env, power: float, theta: dict | None = None) -> float:
    idle, top = envelope(env, theta)
    return float(idle.out["Nc_fan"] + power * (top.out["Nc_fan"] - idle.out["Nc_fan"]))


@dataclass
class Trim:
    sp: SteadyPoint
    jac: np.ndarray
    A: np.ndarray
    B: np.ndarray
    power: float
    label: str = ""


@lru_cache(maxsize=256)
def _trim_cached(key: tuple, power: float) -> Trim:
    alt, MN, dT, th = key
    env, theta = Env(alt, MN, dT), dict(th)
    eng = load_engine()
    idle, top = _idle_cached(key), _top_cached(key)
    target = float(idle.out["Nc_fan"] + power * (top.out["Nc_fan"] - idle.out["Nc_fan"]))
    sp = top if power >= 0.5 else idle
    start = float(sp.out["Nc_fan"])
    n = max(1, int(abs(start - target) / 60.0) + 1)
    for s in np.linspace(start, target, n + 1)[1:]:
        sp = solve_steady(eng, env, target=("Nc_fan", float(s)), theta=theta, tol=1e-11, guess=sp)
    if not sp.info["converged"]:
        raise RuntimeError("no steady state at %.0f %% power" % (100 * power))
    jac = inner_jacobian(eng, sp)
    ab = reference_AB(eng, sp, tol=1e-13, jac=jac)
    return Trim(sp, jac, ab["A"], ab["B"], power,
                "alt %.0f ft, M%.2f, dT %+.0f F, %.0f %% power" % (alt, MN, dT, 100 * power))


def trim_at(env: Env, power: float, theta: dict | None = None) -> Trim:
    """Steady operating point and its linear model at `power` (0 = idle, 1 = max) for this flight condition and health."""
    return _trim_cached(_key(env, theta or {}), round(float(power), 4))


def eig_info(A: np.ndarray, x0) -> dict:
    """Eigenvalues (own 2x2 formula), time constants, kappa(V) of the normalized form."""
    lam, info = mc.eig2x2(A)
    An = normalized_A(A, x0)
    lamn, _ = mc.eig2x2(An)
    ts, tf = mc.time_constants(lam)
    return {"lam": [complex(l) for l in lam], "tau_slow": ts, "tau_fast": tf, "kappaV": mc.kappa2(mc.eigvec2x2(An, lamn)),
            "tr": info["tr"], "det": info["det"], "disc": info["disc"], "stable": bool(mc.is_stable2x2(A))}


# --------------------------------------------------------------------------------------------- controller
class LabController:
    """NASA baseline fuel control (structure + gains + limits from setup_Controller.m), velocity form, with lab switches."""

    def __init__(self, engine, env: Env, kp_mult: float = 1.0, ki_mult: float = 1.0, enabled: dict | None = None):
        m = engine.m
        self.env = env
        self.kp_mult, self.ki_mult = kp_mult, ki_mult
        self.enabled = {k: True for k in LIMITERS}
        self.enabled.update(enabled or {})
        self.kp = Map2D(_v(m["Cntrl.Kp_NcVec"]), _v(m["Cntrl.Kp_AltVec"]), m["Cntrl.Kp_sch"])
        self.ki = Map2D(_v(m["Cntrl.Ki_NcVec"]), _v(m["Cntrl.Ki_AltVec"]), m["Cntrl.Ki_sch"])
        self.kix = Map2D(_v(m["Cntrl.Ki_MNVecX"]), _v(m["Cntrl.Ki_AltVecX"]), m["Cntrl.KiX_sch"])
        self.lim = {"N1max": (float(m["Cntrl.N1_MAX"]), float(m["Cntrl.N1mx_Kp"]), float(m["Cntrl.N1mx_Ki"])),
                    "N3max": (float(m["Cntrl.N3_MAX"]), float(m["Cntrl.N3mx_Kp"]), float(m["Cntrl.N3mx_Ki"])),
                    "Ps3max": (float(m["Cntrl.Ps3_MAX"]), float(m["Cntrl.Ps3mx_Kp"]), float(m["Cntrl.Ps3mx_Ki"])),
                    "T45max": (float(m["Cntrl.T45_MAX"]), float(m["Cntrl.T45mx_Kp"]), float(m["Cntrl.T45mx_Ki"])),
                    "Ps3min": (0.0, float(m["Cntrl.Ps3mn_Kp"]), float(m["Cntrl.Ps3mn_Ki"]))}
        self.ru_nc, self.ru_mx, self.ru_mn = _v(m["Cntrl.RU_Nc"]), _v(m["Cntrl.RUmx"]), float(m["Cntrl.RUmn"])
        self.prev: dict = {}

    def limits(self) -> dict:
        return {"N1max": self.lim["N1max"][0], "N3max": self.lim["N3max"][0], "Ps3max": self.lim["Ps3max"][0],
                "T45max": self.lim["T45max"][0], "Ps3min": ps3_min_limit(self.env.alt_ft)}

    def gains(self, nc: float) -> tuple:
        kp = self.kp(nc, self.env.alt_ft) * self.kp_mult
        ki = self.ki(nc, self.env.alt_ft) * self.kix(self.env.MN, self.env.alt_ft) * self.ki_mult
        return kp, ki

    def _rate(self, name: str, err: float, kp: float, ki: float, dt: float) -> float:
        de = (err - self.prev.get(name, err)) / dt
        self.prev[name] = err
        return kp * de + ki * err

    def step(self, wf: float, nc_dem: float, out: dict, n_hp: float, dt: float) -> tuple:
        """One update: returns (new W_f, active loop name)."""
        kp, ki = self.gains(out["Nc_fan"])
        req = {"fan speed": self._rate("fan", nc_dem - out["Nc_fan"], kp, ki, dt)}
        for name, meas in (("N1max", out["N_fan"]), ("N3max", n_hp), ("Ps3max", out["Ps3"]), ("T45max", out["T45"])):
            lim, kpl, kil = self.lim[name]
            r = self._rate(name, lim - meas, kpl, kil, dt)
            if self.enabled[name]:
                req[name] = r
        active = min(req, key=req.get)
        r = req[active]
        _l, kpl, kil = self.lim["Ps3min"]
        r_min = self._rate("Ps3min", ps3_min_limit(self.env.alt_ft) - out["Ps3"], kpl, kil, dt)
        if self.enabled["Ps3min"] and r_min > r:
            r, active = r_min, "Ps3min"
        wf_new = wf + r * dt
        if self.enabled["RUmax"]:
            cap = float(np.interp(out["Nc_fan"], self.ru_nc, self.ru_mx)) * out["Ps3"]
            if wf_new > cap:
                wf_new, active = cap, "RUmax"
        if self.enabled["RUmin"] and wf_new < self.ru_mn * out["Ps3"]:
            wf_new, active = self.ru_mn * out["Ps3"], "RUmin"
        return max(wf_new, 0.02), active


# --------------------------------------------------------------------------------------------- fixed-length runs
def simulate(env: Env, power0: float, schedule: list, t_end: float, theta: dict | None = None, kp_mult: float = 1.0,
             ki_mult: float = 1.0, enabled: dict | None = None, demand_rate: float = 1.0) -> dict:
    """Closed-loop run: throttle schedule [(t, power 0..1)], demand moves at most `demand_rate` x (max-idle) per second."""
    eng = load_engine()
    theta = theta or {}
    tr = trim_at(env, power0, theta)
    idle, top = envelope(env, theta)
    nc_lo, nc_hi = float(idle.out["Nc_fan"]), float(top.out["Nc_fan"])
    x, z, wf = np.array(tr.sp.x0, float), np.array(tr.sp.z0), float(tr.sp.ctrl.Wf)
    ctl = LabController(eng, env, kp_mult, ki_mult, enabled)
    jac = tr.jac
    nc_dem = float(tr.sp.out["Nc_fan"])
    keys = ("t", "N_LP", "N_HP", "N_fan", "Nc_fan", "Nc_dem", "Wf", "Fnet", "T4", "T45", "Ps3", "VAFN", "VBV",
            "SM_hpc", "SM_lpc", "SM_fan", "W_in", "loop", "RU_cap")
    rec = {k: [] for k in keys}
    n = int(round(t_end / DT))
    for k in range(n + 1):
        t = k * DT
        pw = [p for (ts, p) in schedule if ts <= t + 1e-9][-1] if schedule else power0
        target = nc_lo + pw * (nc_hi - nc_lo)
        step = demand_rate * (nc_hi - nc_lo) * DT
        nc_dem += max(-step, min(step, target - nc_dem))
        ncf = eng.Nc_fan_of(x[0], env)
        ctrl = Controls(wf, eng.schedule_VAFN(ncf, env.MN), eng.schedule_VBV(ncf, env.MN))
        z, out, info = eng.solve_flow(x, ctrl, env, theta, tol=1e-9, z0=z, jac=jac)
        if not info["converged"]:
            z, out, info = eng.solve_flow(x, ctrl, env, theta, tol=1e-9, z0=z)
            if not info["converged"]:
                rec["error"] = "the engine model could not be solved at t = %.2f s (state left the map range)" % t
                break
        wf_next, loop = ctl.step(wf, nc_dem, out, x[1], DT)
        row = dict(out, t=t, N_LP=x[0], N_HP=x[1], Nc_dem=nc_dem, loop=loop,
                   RU_cap=float(np.interp(out["Nc_fan"], ctl.ru_nc, ctl.ru_mx)) * out["Ps3"])
        for key in keys:
            rec[key].append(row[key])
        x = x + DT * np.array([out["Ndot_LP"], out["Ndot_HP"]])
        wf = wf_next
    rec["limits"] = ctl.limits()
    rec["nc_idle"], rec["nc_max"] = nc_lo, nc_hi
    return rec


# --------------------------------------------------------------------------------------------- live simulation
DEFAULTS = {
    "mode": "auto",            # "auto" = NASA controller follows the power lever; "manual" = you set the fuel flow
    "power": 60.0,             # power lever, % (0 = idle, 100 = max)
    "wf_manual": 0.8,          # lbm/s, manual fuel flow
    "alt": 0.0, "mach": 0.0, "dT": 0.0,
    "vafn_mode": "schedule", "vafn": 6000.0, "vbv_mode": "schedule", "vbv": 0.0,
    "health": {k: 0.0 for k in HEALTH_KEYS},  # % degradation per parameter
    "kp_mult": 1.0, "ki_mult": 1.0, "limiters": list(LIMITERS), "speed": 1.0,
}


class LiveSim:
    """The open-ended live simulation. One instance per app process; all methods are thread-safe."""

    WINDOW = 60.0             # seconds of history kept for the strip charts

    def __init__(self):
        self.lock = threading.RLock()
        self.eng = load_engine()
        self.settings = dict(DEFAULTS, health=dict(DEFAULTS["health"]))
        self.running = False
        self.message = ""
        self.reset()

    # ---------------- helpers
    def _env_target(self) -> Env:
        s = self.settings
        return Env(float(s["alt"]), float(s["mach"]), float(s["dT"]))

    def _theta(self) -> dict:
        return health_dict(0.0, self.settings["health"])

    def reset(self) -> str:
        """Re-trim at the current settings (steady start), clear history. Returns a status message."""
        with self.lock:
            s = self.settings
            env = self._env_target()
            theta = self._theta()
            try:
                if s["mode"] == "auto":
                    tr = trim_at(env, s["power"] / 100.0, theta)
                else:
                    tr = trim_at(env, 0.6, theta)
                self.message = ""
            except Exception as exc:            # fall back to sea level so the app always works
                self.message = "Could not start at these settings (%s); reset to sea level." % exc
                s.update(alt=0.0, mach=0.0, dT=0.0)
                env = Env(0.0, 0.0, 0.0)
                tr = trim_at(env, s["power"] / 100.0, theta)
            sp = tr.sp
            self.env = env
            self.theta = theta
            self.x = np.array(sp.x0, float)
            self.z = np.array(sp.z0)
            self.wf = float(sp.ctrl.Wf)
            self.vafn, self.vbv = float(sp.ctrl.VAFN), float(sp.ctrl.VBV)
            self.jac = tr.jac
            self.ctl = LabController(self.eng, env, s["kp_mult"], s["ki_mult"], {k: (k in s["limiters"]) for k in LIMITERS})
            self.nc_dem = float(sp.out["Nc_fan"])
            self.loop = "fan speed"
            self.t = 0.0
            self.out = dict(sp.out)
            self.hist = {k: [] for k in ("t", "Nfan_pct", "NHP_pct", "Nc", "Nc_dem", "Wf", "Fnet", "T45", "Ps3", "T4")}
            self.new_rows = []
            self.env_lo_hi = self._envelope_nc(env, theta)
            if s["mode"] == "manual":
                s["wf_manual"] = round(self.wf, 3)
            self._record()
            self.last_wall = None
            return self.message or "Ready: steady state at %s." % tr.label

    def _envelope_nc(self, env: Env, theta: dict) -> tuple:
        try:
            idle, top = envelope(Env(round(env.alt_ft, -2), round(env.MN, 2), round(env.dT)), theta)
            return float(idle.out["Nc_fan"]), float(top.out["Nc_fan"])
        except Exception:
            return 0.4 * S_NC, nc_max(env)

    def update(self, new: dict) -> str:
        """Apply new settings while running. Flight condition and health change smoothly (slew), others at once."""
        with self.lock:
            s = self.settings
            old_env = (s["alt"], s["mach"], s["dT"])
            s.update({k: v for k, v in new.items() if v is not None})
            self.ctl.kp_mult, self.ctl.ki_mult = float(s["kp_mult"]), float(s["ki_mult"])
            self.ctl.enabled = {k: (k in s["limiters"]) for k in LIMITERS}
            if (s["alt"], s["mach"], s["dT"]) != old_env:
                # envelope for the new flight condition (cached after the first time)
                try:
                    self.env_lo_hi_target = self._envelope_nc(self._env_target(), self._theta())
                    self.message = ""
                except Exception as exc:
                    self.message = "Flight condition not reachable (%s); kept the previous one." % exc
                    s["alt"], s["mach"], s["dT"] = old_env
            if s["mode"] == "manual" and new.get("mode") == "manual" and new.get("wf_manual") is None:
                s["wf_manual"] = round(self.wf, 3)
            return self.message

    def fuel_pulse(self, pct: float) -> None:
        """Manual mode: change the fuel flow by `pct` percent at once (a step input)."""
        with self.lock:
            self.settings["mode"] = "manual"
            self.settings["wf_manual"] = round(self.wf * (1 + pct / 100.0), 4)

    def _slew(self, dt: float) -> None:
        s, e = self.settings, self.env
        def mv(cur, tgt, rate):
            return cur + max(-rate * dt, min(rate * dt, tgt - cur))
        alt = mv(e.alt_ft, float(s["alt"]), 2500.0)          # ft/s (fast, to keep experiments short)
        mach = mv(e.MN, float(s["mach"]), 0.05)               # 1/s
        dT = mv(e.dT, float(s["dT"]), 5.0)                    # F/s
        if (alt, mach, dT) != (e.alt_ft, e.MN, e.dT):
            self.env = Env(alt, mach, dT)
            self.ctl.env = self.env
            if (alt, mach, dT) == (float(s["alt"]), float(s["mach"]), float(s["dT"])):
                self.env_lo_hi = getattr(self, "env_lo_hi_target", self.env_lo_hi)
        th_t = self._theta()
        keys = set(th_t) | set(self.theta)
        th = {k: mv(self.theta.get(k, 1.0), th_t.get(k, 1.0), 0.02) for k in keys}  # 2 % per second
        self.theta = {k: v for k, v in th.items() if v != 1.0}

    def _record(self) -> None:
        o = self.out
        row = {"t": self.t, "Nfan_pct": 100 * o["N_fan"] / 2300.0, "NHP_pct": 100 * self.x[1] / 22500.0, "Nc": o["Nc_fan"],
               "Nc_dem": self.nc_dem, "Wf": self.wf, "Fnet": o["Fnet"], "T45": o["T45"], "Ps3": o["Ps3"], "T4": o["T4"]}
        for k, v in row.items():
            self.hist[k].append(v)
        if len(self.hist["t"]) > int(self.WINDOW / DT) + 50:
            for k in self.hist:
                del self.hist[k][:50]
        self.new_rows.append(row)

    def _one_step(self) -> bool:
        s = self.settings
        self._slew(DT)
        env = self.env
        ncf = self.eng.Nc_fan_of(self.x[0], env)
        # actuators: NASA schedules or manual, with rate limits
        vafn_t = self.eng.schedule_VAFN(ncf, env.MN) if s["vafn_mode"] == "schedule" else float(s["vafn"])
        vbv_t = self.eng.schedule_VBV(ncf, env.MN) if s["vbv_mode"] == "schedule" else float(s["vbv"])
        self.vafn += max(-3000 * DT, min(3000 * DT, vafn_t - self.vafn))
        self.vbv += max(-1.0 * DT, min(1.0 * DT, vbv_t - self.vbv))
        ctrl = Controls(self.wf, self.vafn, self.vbv)
        z, out, info = self.eng.solve_flow(self.x, ctrl, env, self.theta, tol=1e-9, z0=self.z, jac=self.jac)
        if not info["converged"]:
            self.jac = self.eng.flow_jacobian(self.z, self.x, ctrl, env, self.theta)
            z, out, info = self.eng.solve_flow(self.x, ctrl, env, self.theta, tol=1e-9, z0=self.z, jac=self.jac)
        if not info["converged"]:
            self.running = False
            self.message = ("Paused: the engine model cannot be solved at this state (it has left NASA's component maps). "
                            "Undo the last change or press Reset.")
            return False
        self.z, self.out = z, out
        lo, hi = self.env_lo_hi
        if s["mode"] == "auto":
            target = lo + float(s["power"]) / 100.0 * (hi - lo)
            rate = (hi - lo) * 1.0 * DT                     # full lever travel takes 1 s
            self.nc_dem += max(-rate, min(rate, target - self.nc_dem))
            self.wf, self.loop = self.ctl.step(self.wf, self.nc_dem, out, self.x[1], DT)
        else:
            self.nc_dem = out["Nc_fan"]
            self.wf += max(-1.0 * DT, min(1.0 * DT, float(s["wf_manual"]) - self.wf))   # fuel valve 1 lbm/s per s
            self.loop = "manual fuel"
        self.x = self.x + DT * np.array([out["Ndot_LP"], out["Ndot_HP"]])
        self.t += DT
        if self.x[0] < 0.25 * 6772.0 or self.x[1] < 0.55 * 20796.0:
            self.running = False
            self.message = "Paused: shaft speeds fell below the model's valid range (no sub-idle or start-up model). Press Reset."
            return False
        self._record()
        return True

    def tick(self, budget_s: float = 0.08) -> list:
        """Advance by the wall time since the last tick x speed (bounded by a compute budget). Returns new rows."""
        with self.lock:
            now = time.perf_counter()
            want = 0.1 if self.last_wall is None else min(now - self.last_wall, 0.3)
            self.last_wall = now
            self.new_rows = []
            if not self.running:
                return []
            n = max(1, int(round(want * float(self.settings["speed"]) / DT)))
            t0 = time.perf_counter()
            for _ in range(n):
                if not self._one_step() or time.perf_counter() - t0 > budget_s:
                    break
            rows, self.new_rows = self.new_rows, []
            return rows

    def snapshot(self) -> dict:
        """Everything the animated engine view and the gauges need."""
        with self.lock:
            o = self.out
            lim = self.ctl.limits()
            return {
                "t": self.t, "running": self.running, "message": self.message, "loop": self.loop, "mode": self.settings["mode"],
                "N_fan": o["N_fan"], "N_LP": float(self.x[0]), "N_HP": float(self.x[1]), "Nc_fan": o["Nc_fan"], "Nc_dem": self.nc_dem,
                "Wf": self.wf, "Fnet": o["Fnet"], "Fg_byp": o["Fg_byp"], "Fg_core": o["Fg_core"], "W_in": o["W_in"],
                "W_byp": o["W_byp"], "W_core": o["W_core"], "W_vbv": o["W_vbv"], "W5": o["W5"], "BPR": o["BPR"],
                "Tt2": o["Tt2"], "Tt21": o["Tt21"], "Tt25": o["Tt25"], "Tt3": o["Tt3"], "T4": o["T4"], "T45": o["T45"], "Tt5": o["Tt5"],
                "Pt2": o["Pt2"], "Pt3": o["Pt3"], "Ps3": o["Ps3"], "Pamb": o["Pamb"], "OPR": o["OPR"], "TSFC": o["TSFC"],
                "SM_fan": o["SM_fan"], "SM_lpc": o["SM_lpc"], "SM_hpc": o["SM_hpc"], "VAFN": self.vafn, "VBV": self.vbv,
                "alt": self.env.alt_ft, "mach": self.env.MN, "dT": self.env.dT,
                "health": {k: round(100 * (1 - v), 3) for k, v in self.theta.items()},
                "limits": lim, "power": float(self.settings["power"]),
                "dNdt_rel": float(max(abs(o["Ndot_LP"]) / self.x[0], abs(o["Ndot_HP"]) / self.x[1])),
            }

    def linearize_here(self) -> dict:
        """Linear model at the CURRENT state (inputs frozen). Meaningful when the engine is steady."""
        with self.lock:
            sp = SteadyPoint(env=self.env, ctrl=Controls(self.wf, self.vafn, self.vbv), x0=self.x.copy(), z0=self.z.copy(),
                             out=dict(self.out), theta=dict(self.theta), info={})
            ab = reference_AB(self.eng, sp, tol=1e-12)
            res = eig_info(ab["A"], sp.x0)
            res.update(A=ab["A"], B=ab["B"].ravel(), dNdt_rel=float(np.max(np.abs(ab["F0"]) / sp.x0)), t=self.t)
            return res
