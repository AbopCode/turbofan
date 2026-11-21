"""Python port of NASA's AGTF30 geared-turbofan engine model (T-MATS components).

What is modelled (wiring read from `vendor/agtf30/AGTF30/AGTF30_eng.mdl`):

    Ambient -> Inlet -> FAN -> Splitter --(bypass)--> [+VBV bleed] -> duct17 -> bypass nozzle (area = VAFN)
                                  \\--(core)--> duct2 -> LPC -> VBV valve -> duct25 -> HPC -> Burner
                                               -> HPT -> duct45 -> LPT -> duct5 -> core nozzle
    FAN --gearbox 3.1:1--> LP shaft (LPC + LPT)          HP shaft (HPC + HPT, 350 hp power off-take)

The component equations are ported from the T-MATS C sources
(`vendor/tmats/TMATS_Library/MEX/*_TMATS_body.c`) and library block diagrams.
The maps are the real AGTF30 maps (see `agtf30_data.py`).

The model has NO volume dynamics: at every instant the seven "flow-matching"
equations r(z) = 0 are solved for the seven unknowns

    z = [W, HPT_PR, HPC_Rline, FAN_Rline, BPR, LPC_Rline, LPT_PR]      (NASA's `ICcom` order)

and only then are the two shaft accelerations known:

    I_LP (pi/30) dN_LP/dt = net LP torque,     I_HP (pi/30) dN_HP/dt = net HP torque.

So dx/dt = F(x, u; theta) with x = [N_LP, N_HP], u = W_f is an *implicit*
function: every evaluation of F runs the inner Newton solver `solve_flow`,
whose tolerance `tol` is a study variable.

Two switches (see `load_engine`: mode='nasa' = faithful port, mode='study' = differentiable version):
  * map_interp='linear' : bilinear map interpolation as in T-MATS (derivatives jump at every map
    breakpoint - and NASA's design point sits exactly on one);  map_interp='smooth' : C1 bicubic
    interpolation through the same map points (DERIVED deviation).
Iteration modes:
  * iter_mode='tmats' : every inner loop of the C code (gas-property inversions, duct static
    pressure search, nozzle Mach-1 search) is reproduced with NASA's thresholds and iteration
    caps. Faithful, but F then carries a small non-smooth error of its own (about 1e-6 relative).
  * iter_mode='tight' (default): the same equations, with those loops converged to machine
    precision (DERIVED deviation), so that F is smooth and `tol` is the only solver noise.
The Ambient block always uses the faithful loop (it depends on the flight condition only).

Course topic: Engine model (not a course topic)
"""
from __future__ import annotations

import math
from bisect import bisect_right
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import tmats_thermo as th
from .agtf30_data import load_params

# ---- T-MATS constants (constants_TMATS.h) ---------------------------------
C_PSTD = 14.695951            # psia
C_TSTD = 518.67               # deg R
C_BTU_PER_SECtoHP = 1.414853
C_HP_PER_RPMtoFT_LBF = 5252.113
C_GRAVITY = 32.174049
JOULES_CONST = 778.169233
C_PSItoPSF = 144.0
C_GAMMA = 1.4
RT_GAS = 0.0686               # BTU/(lbm R): the Rt table of every AGTF30 block is constant
GJ2 = 2.0 * C_GRAVITY * JOULES_CONST
SHAFT_PWR2TRQ = -5252.11      # gain inside the T-MATS Shaft block (power off-take -> torque)
G_AMBIENT = 32.1740           # constant "g" inside the Ambient block (ram drag)
HP_POWER_DEBIT = 350.0        # hp, constant "PwrDebit [HP]" on the HP shaft in AGTF30_eng.mdl

HEALTH_KEYS = ("fan_eff", "fan_flow", "lpc_eff", "lpc_flow", "hpc_eff", "hpc_flow",
               "hpt_eff", "hpt_flow", "lpt_eff", "lpt_flow")
Z_NAMES = ("W", "HPT_PR", "HPC_Rline", "FAN_Rline", "BPR", "LPC_Rline", "LPT_PR")
R_NAMES = ("FAN", "LPC", "HPC", "HPT", "LPT", "CoreNoz", "BypNoz")    # order of the 7 residuals


@dataclass(frozen=True)
class Env:
    """Flight condition: altitude [ft], Mach number, ambient temperature offset dT [deg F = deg R]."""
    alt_ft: float
    MN: float
    dT: float = 0.0


@dataclass(frozen=True)
class Controls:
    """Engine inputs: fuel flow Wf [lbm/s], bypass nozzle throat area VAFN [in^2], bleed valve VBV [0..1]."""
    Wf: float
    VAFN: float
    VBV: float


# ---------------------------------------------------------------------------
# Small helpers ported from functions_TMATS.c / interp*Ac_TMATS.c
# ---------------------------------------------------------------------------
def _divby(x: float) -> float:
    """1/x with T-MATS's protection against division by ~0 (divby())."""
    if -1e-10 < x < 1e-10:
        return 1e10 if x >= 0 else -1e10
    return 1.0 / x


def _sqrt(x: float) -> float:
    """sqrt with T-MATS's protection against negative arguments (sqrtT())."""
    return math.sqrt(x) if x > 0.0 else 0.0


def _interp1(X: list, Y: list, xi: float) -> float:
    """1-D linear interpolation, clamped at the table ends (interp1Ac())."""
    n = len(X)
    if xi < X[0] or not (xi >= X[0]):
        xi = X[0]
    elif xi > X[n - 1]:
        xi = X[n - 1]
    i = bisect_right(X, xi) - 1
    if i > n - 2:
        i = n - 2
    return Y[i] + (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]) * (xi - X[i])


def _interp1_extrap(X: list, Y: list, xi: float) -> float:
    """1-D linear interpolation with linear extrapolation (Simulink n-D Lookup, used by the Inlet)."""
    n = len(X)
    i = bisect_right(X, xi) - 1
    i = 0 if i < 0 else (n - 2 if i > n - 2 else i)
    return Y[i] + (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]) * (xi - X[i])


def _node_slopes(X: list, V: list) -> list:
    """Slopes at the table nodes for C1 cubic interpolation (three-point / Catmull-Rom formula).

    Interior node: weighted mean of the two neighbouring secant slopes; end nodes: the end secant.
    The formula is LINEAR in the data, so the 2-D interpolant built from it is C1 everywhere.
    """
    n = len(X)
    sec = [(V[k + 1] - V[k]) / (X[k + 1] - X[k]) for k in range(n - 1)]
    d = [sec[0]] + [0.0] * (n - 2) + [sec[-1]]
    for k in range(1, n - 1):
        h0, h1 = X[k] - X[k - 1], X[k + 1] - X[k]
        d[k] = (sec[k - 1] * h1 + sec[k] * h0) / (h0 + h1)
    return d


def _hermite(p0: float, p1: float, m0: float, m1: float, t: float, h: float) -> float:
    """Cubic Hermite polynomial on one cell: values p0, p1 and slopes m0, m1 at the ends, t in [0, 1]."""
    t2 = t * t
    t3 = t2 * t
    return ((2.0 * t3 - 3.0 * t2 + 1.0) * p0 + (t3 - 2.0 * t2 + t) * h * m0
            + (3.0 * t2 - 2.0 * t3) * p1 + (t3 - t2) * h * m1)


class Map2D:
    """Table Z[row = y index][col = x index], clamped at the edges.

    smooth=False : bilinear interpolation, exactly T-MATS's interp2Ac() (kinks at every breakpoint).
    smooth=True  : C1 bicubic Hermite interpolation through the SAME table points (DERIVED deviation),
                   so that partial derivatives of the engine model exist everywhere inside the map.
    """

    __slots__ = ("X", "Y", "Z", "nx", "ny", "smooth", "DX", "DY", "DXY")

    def __init__(self, X, Y, Z, smooth: bool = False):
        self.X = [float(v) for v in X]
        self.Y = [float(v) for v in Y]
        self.Z = [[float(v) for v in row] for row in Z]
        self.nx, self.ny = len(self.X), len(self.Y)
        if len(self.Z) != self.ny or len(self.Z[0]) != self.nx:
            raise ValueError("map shape %dx%d does not match axes %d,%d" % (len(self.Z), len(self.Z[0]), self.ny, self.nx))
        self.smooth = smooth
        if smooth:
            self.DX = [_node_slopes(self.X, row) for row in self.Z]                       # dZ/dx at nodes
            cols = [_node_slopes(self.Y, [self.Z[j][i] for j in range(self.ny)]) for i in range(self.nx)]
            self.DY = [[cols[i][j] for i in range(self.nx)] for j in range(self.ny)]       # dZ/dy at nodes
            self.DXY = [_node_slopes(self.X, row) for row in self.DY]                      # d2Z/dxdy at nodes

    def locate(self, xi: float, yi: float):
        """Clamp (xi, yi) and return the cell indices and fractions; shared by several maps with the same axes."""
        X, Y = self.X, self.Y
        if xi < X[0] or not (xi >= X[0]):
            xi = X[0]
        elif xi > X[-1]:
            xi = X[-1]
        if yi < Y[0] or not (yi >= Y[0]):
            yi = Y[0]
        elif yi > Y[-1]:
            yi = Y[-1]
        i = bisect_right(X, xi) - 1
        if i > self.nx - 2:
            i = self.nx - 2
        j = bisect_right(Y, yi) - 1
        if j > self.ny - 2:
            j = self.ny - 2
        return i, j, (xi - X[i]) / (X[i + 1] - X[i]), (yi - Y[j]) / (Y[j + 1] - Y[j])

    def at(self, loc) -> float:
        i, j, fx, fy = loc
        r0, r1 = self.Z[j], self.Z[j + 1]
        if not self.smooth:
            z1 = r0[i] + (r0[i + 1] - r0[i]) * fx
            z2 = r1[i] + (r1[i + 1] - r1[i]) * fx
            return z1 + (z2 - z1) * fy
        hx = self.X[i + 1] - self.X[i]
        hy = self.Y[j + 1] - self.Y[j]
        dx0, dx1, dy0, dy1, c0, c1 = self.DX[j], self.DX[j + 1], self.DY[j], self.DY[j + 1], self.DXY[j], self.DXY[j + 1]
        v0 = _hermite(r0[i], r0[i + 1], dx0[i], dx0[i + 1], fx, hx)       # value along row j
        v1 = _hermite(r1[i], r1[i + 1], dx1[i], dx1[i + 1], fx, hx)       # value along row j+1
        g0 = _hermite(dy0[i], dy0[i + 1], c0[i], c0[i + 1], fx, hx)       # dZ/dy along row j
        g1 = _hermite(dy1[i], dy1[i + 1], c1[i], c1[i + 1], fx, hx)       # dZ/dy along row j+1
        return _hermite(v0, v1, g0, g1, fy, hy)

    def __call__(self, xi: float, yi: float) -> float:
        return self.at(self.locate(xi, yi))


def _vec(v) -> list:
    return [float(x) for x in np.asarray(v, dtype=float).reshape(-1)]


class _CompMap:
    """One compressor: maps Wc, PR, Eff = f(Rline, NcMap), map scalars, surge line."""

    def __init__(self, m: dict, name: str, smooth: bool):
        g = lambda k: m["%s.%s" % (name, k)]
        R, Nc = _vec(g("RVec")), _vec(g("NcVec"))
        self.Wc = Map2D(R, Nc, g("WcArray"), smooth)
        self.PR = Map2D(R, Nc, g("PRArray"), smooth)
        self.Eff = Map2D(R, Nc, g("EffArray"), smooth)
        self.s_Nc, self.s_Wc, self.s_PR, self.s_Eff = (float(g(k)) for k in ("s_Nc", "s_Wc", "s_PR", "s_Eff"))
        self.SM_Wc, self.SM_PR = _vec(g("SM_WcVec")), _vec(g("SM_PRVec"))


class _TurbMap:
    """One turbine: maps Wc, Eff = f(PRmap, NcMap), map scalars, cooling-flow positions."""

    def __init__(self, m: dict, name: str, smooth: bool):
        g = lambda k: m["%s.%s" % (name, k)]
        PR, Nc = _vec(g("PRVec")), _vec(g("NcVec"))
        self.Wc = Map2D(PR, Nc, g("WcArray"), smooth)
        self.Eff = Map2D(PR, Nc, g("effArray"), smooth)
        self.s_Nc, self.s_Wc, self.s_PR, self.s_Eff = (float(g(k)) for k in ("s_Nc", "s_Wc", "s_PR", "s_Eff"))
        self.pos = _vec(g("CWPos"))


class Engine:
    """The AGTF30 engine: flow-matching residuals, inner solver and shaft dynamics.

    Course topic: Engine model (not a course topic)
    """

    def __init__(self, params: dict, iter_mode: str = "tight", map_interp: str = "smooth"):
        if iter_mode not in ("tight", "tmats"):
            raise ValueError("iter_mode must be 'tight' or 'tmats'")
        if map_interp not in ("smooth", "linear"):
            raise ValueError("map_interp must be 'smooth' or 'linear'")
        self.iter_mode = iter_mode
        self.map_interp = map_interp
        self.tight = iter_mode == "tight"
        sm = map_interp == "smooth"
        m = params["MWS"]
        self.m = m
        self.fan, self.lpc, self.hpc = _CompMap(m, "FAN", sm), _CompMap(m, "LPC", sm), _CompMap(m, "HPC", sm)
        self.hpt, self.lpt = _TurbMap(m, "HPT", sm), _TurbMap(m, "LPT", sm)
        self.hpc_FBD, self.hpc_BFht, self.hpc_BFPt = _vec(m["HPC.FBD"]), _vec(m["HPC.BFht"]), _vec(m["HPC.BFPt"])
        self.GR = float(m["GearBox.GearRatio"])
        self.gb_eff = float(m["GearBox.Eff"])
        self.I_LP = float(m["Shaft.LPS_Inertia"])
        self.I_HP = float(m["Shaft.HPS_Inertia"])
        self.lps_eff = float(m["Shaft.LPS_Eff"])
        self.k_LP = 60.0 / (2.0 * math.pi * self.I_LP)      # rpm/s per ft lbf  (T-MATS Shaft gain)
        self.k_HP = 60.0 / (2.0 * math.pi * self.I_HP)
        self.eRam = float(m["Inlet.eRam"])
        self.eRamVec, self.eRamTBL = _vec(m["Inlet.eRamVec"]), _vec(m["Inlet.eRamTBL"])
        self.duct = {k: (float(m["Duct.D" + k]), float(m["Duct.MNdes" + k]), float(m["Duct.Ath" + k]))
                     for k in ("2", "25", "45", "5", "17")}
        self.hpc_Ath = float(m["HPC.Ath"])
        self.burn_LHV, self.burn_eff, self.burn_dP = float(m["Burn.LHV"]), float(m["Burn.Eff"]), float(m["Burn.dPqP_dmd"])
        self.noz_core_Ath, self.noz_core_Cfg = float(m["NozCor.Ath"]), float(m["NozCor.Cfg"])
        self.noz_byp_Cfg = float(m["NozByp.Cfg"])
        self.vbv_PR, self.vbv_flow = _vec(m["VBV.PRVec"]), _vec(m["VBV.Flowe"])
        self.vbv_Ae, self.vbv_open, self.vbv_closed = float(m["VBV.Ae"]), float(m["VBV.fullopen"]), float(m["VBV.fullclosed"])
        amb = params["ambient"]
        self.alt_vec, self.ts_vec, self.ps_vec = _vec(amb["AltVec"]), _vec(amb["TsVec"]), _vec(amb["PsVec"])
        # NASA baseline schedules VAFN, VBV = f(corrected fan speed, Mach)
        self._vafn = Map2D(_vec(m["Cntrl.VAFN_MN"]), _vec(m["Cntrl.VAFN_Nc1"]), m["Cntrl.VAFN_sch"])
        self._vbv = Map2D(_vec(m["Cntrl.VBV_MN"]), _vec(m["Cntrl.VBV_Nc1"]), m["Cntrl.VBV_sch"])
        self._amb_cache: dict = {}
        self._h2t = th.h2t_tight if self.tight else th.h2t

    # ------------------------------------------------------------------
    # gas-path building blocks
    # ------------------------------------------------------------------
    def _sp2t(self, S: float, P: float, fa: float, guess: float) -> float:
        """Temperature from entropy and pressure (mode dependent)."""
        return th.sp2t_tight(S, P, fa, guess) if self.tight else th.sp2t(S, P, fa)

    def _pcalc(self, S: float, Ps: float, Tt: float, Pt: float, ht: float, fa: float):
        """Static state for a guessed static pressure (port of PcalcStat()): returns Ts, rho_s, V."""
        Ts = self._sp2t(S, Ps, fa, Tt * (Ps / Pt) ** 0.2857 if Ps < Pt else Tt)
        if Ts > Tt:
            Ts = Tt
        hs = th.t2h(Ts, fa)
        if hs > ht:
            hs = ht
        rhos = Ps * C_PSItoPSF * _divby(RT_GAS * Ts * JOULES_CONST)
        return Ts, rhos, _sqrt((ht - hs) * GJ2)

    def ambient(self, env: Env) -> tuple:
        """Ambient block (port of Ambient_TMATS_body.c): returns (Tt0, Pt0, Ps0, Ts0, Veng). Cached per env."""
        key = (env.alt_ft, env.MN, env.dT)
        hit = self._amb_cache.get(key)
        if hit is not None:
            return hit
        fa = 0.0
        Ts = _interp1(self.alt_vec, self.ts_vec, env.alt_ft) + env.dT
        Ps = _interp1(self.alt_vec, self.ps_vec, env.alt_ft)
        S = th.pt2s(Ps, Ts, fa)
        hs = th.t2h(Ts, fa)
        Ttg = Ts * (1.0 + env.MN * env.MN * (C_GAMMA - 1.0) / 2.0)
        Ptg = Ps * _divby((Ts * _divby(Ttg)) ** (C_GAMMA / (C_GAMMA - 1.0)))
        Ttg = th.sp2t(S, Ptg, fa)
        htg = th.t2h(Ttg, fa)
        Vsg = _sqrt(C_GAMMA * RT_GAS * Ts * C_GRAVITY * JOULES_CONST)
        er = env.MN - _sqrt((htg - hs) * GJ2) * _divby(Vsg)
        Ptg_new = Ptg + 0.05
        it = 0
        while abs(er) > 0.001 and it < 15:
            er_old, Ptg_old = er, Ptg
            Ptg = Ptg + 0.05 if abs(Ptg - Ptg_new) < 0.03 else Ptg_new
            Ttg = th.sp2t(S, Ptg, fa)
            htg = th.t2h(Ttg, fa)
            er = env.MN - _sqrt((htg - hs) * GJ2) * _divby(Vsg)
            if abs(er) > 0.001:
                Ptg_new = Ptg - er * (Ptg - Ptg_old) * _divby(er - er_old)
            it += 1
        out = (Ttg, Ptg, Ps, Ts, Vsg * env.MN)
        self._amb_cache[key] = out
        return out

    def _compressor(self, W, Tt, Pt, fa, N, Rline, cm: _CompMap, th_flow, th_eff, bleeds):
        """Compressor (port of Compressor_TMATS_body.c, off-design branch, fractional bleeds).

        Returns (W_out, Tt_out, Pt_out, torque, NErr, stall margin %, bleed list, info tuple).
        """
        htin = th.t2h(Tt, fa)
        S = th.pt2s(Pt, Tt, fa)
        theta_c = Tt / C_TSTD
        Wcin = W * _sqrt(theta_c) * _divby(Pt / C_PSTD)
        Nc = N * _divby(_sqrt(theta_c))
        loc = cm.Wc.locate(Rline, Nc * _divby(cm.s_Nc))
        WcMap = cm.Wc.at(loc)
        WcCalc = WcMap * cm.s_Wc * th_flow          # health scalar on flow capacity
        PRMap = cm.PR.at(loc)
        PR = cm.s_PR * (PRMap - 1.0) + 1.0
        Eff = cm.Eff.at(loc) * cm.s_Eff * th_eff    # health scalar on efficiency
        PtOut = Pt * PR
        TtId = self._sp2t(S, PtOut, fa, Tt * PR ** 0.2857)
        htOut = (th.t2h(TtId, fa) - htin) * _divby(Eff) + htin
        TtOut = self._h2t(htOut, fa)
        Wbleeds = 0.0
        PwrBld = 0.0
        blds = []
        if bleeds is not None:
            for frac, fht, fpt in bleeds:
                if frac <= 0.0:
                    blds.append((0.0, 0.0, 0.0, 0.0))
                    continue
                Wb = frac * W
                Wbleeds += Wb
                hb = htin + fht * (htOut - htin)
                blds.append((Wb, self._h2t(hb, fa), Pt + fpt * (PtOut - Pt), fa))
                PwrBld += Wb * (hb - htOut) * C_BTU_PER_SECtoHP
        Pwr = W * (htin - htOut) * C_BTU_PER_SECtoHP - PwrBld
        Trq = C_HP_PER_RPMtoFT_LBF * Pwr * _divby(N)
        NErr = 100.0 if W == 0 else (Wcin - WcCalc) * _divby(Wcin)
        SPR = cm.s_PR * (_interp1(cm.SM_Wc, cm.SM_PR, WcMap) - 1.0) + 1.0
        SM = (SPR - PR) * _divby(PR) * 100.0
        return W - Wbleeds, TtOut, PtOut, Trq, NErr, SM, blds, (Nc, PR, Eff, Wcin)

    def _turbine(self, W, Tt, Pt, fa, N, PRin, tm: _TurbMap, th_flow, th_eff, cool):
        """Turbine (port of Turbine_TMATS_body.c, off-design, NPSS-style corrected values, cooling flows).

        cool = list of (W, Tt, Pt, FAR) streams; tm.pos[i] = 0 -> enters before the first rotor (does work),
        1 -> enters at the exit (no work). Returns (W_out, Tt_out, Pt_out, FAR_out, torque, NErr, info).
        """
        Wcs1 = Wcout = Wfs1 = Wfout = dHs1 = dHout = 0.0
        for (Wc, Ttc, _Ptc, fac), pos in zip(cool, tm.pos):
            hc = th.t2h(Ttc, fac)
            Wcs1 += Wc * (1.0 - pos)
            Wcout += Wc
            Wfs1 += fac * Wc * (1.0 - pos) * _divby(1.0 + fac)
            Wfout += fac * Wc * _divby(1.0 + fac)
            dHs1 += hc * Wc * (1.0 - pos)
            dHout += hc * Wc * pos
        Ws1 = W + Wcs1
        WOut = W + Wcout
        air = W * _divby(1.0 + fa)
        FARs1 = (fa * air + Wfs1) * _divby(air + Wcs1 - Wfs1)
        FARout = (fa * air + Wfout) * _divby(air + Wcout - Wfout)
        htin = th.t2h(Tt, fa)
        hts1 = (htin * W + dHs1) * _divby(Ws1)
        Tts1 = self._h2t(hts1, FARs1)
        Ss1 = th.pt2s(Pt, Tts1, FARs1)
        Nc = N * _divby(_sqrt(Tt))                      # NPSS convention: no standard-day reference
        PRmap = (PRin - 1.0) * _divby(tm.s_PR) + 1.0
        PtOut = Pt * _divby(PRin)
        loc = tm.Wc.locate(PRmap, Nc * _divby(tm.s_Nc))
        WcCalc = tm.Wc.at(loc) * tm.s_Wc * th_flow      # health scalar on flow capacity
        Wcin = W * _sqrt(Tt) * _divby(Pt)
        Eff = tm.Eff.at(loc) * tm.s_Eff * th_eff        # health scalar on efficiency
        TtId = self._sp2t(Ss1, PtOut, FARout, Tts1 * (1.0 / PRin) ** 0.25)
        htId = th.t2h(TtId, FARout)
        Pwr = (hts1 - htId) * Eff * Ws1 * C_BTU_PER_SECtoHP
        htOut = (((htId - hts1) * Eff + hts1) * Ws1 + dHout) * _divby(WOut)
        TtOut = self._h2t(htOut, FARout)
        Trq = C_HP_PER_RPMtoFT_LBF * Pwr * _divby(N)
        NErr = 100.0 if Ws1 == 0 else (Wcin - WcCalc) * _divby(Wcin)
        return WOut, TtOut, PtOut, FARout, Trq, NErr, (Nc, Eff, Wcin)

    def _static_area(self, W, Tt, Pt, fa, Ath, MN_guess):
        """Static conditions in a duct of known area (port of StaticCalc_TMATS_body.c, SolveType 0).

        Returns (Ts, Ps, MN). 'tmats' mode copies NASA's creeping secant search (stops at 1e-4 area
        error, steps limited to 0.1-0.2 %); 'tight' mode solves the same equation to machine precision.
        """
        S = th.pt2s(Pt, Tt, fa)
        ht = th.t2h(Tt, fa)
        Tsg = Tt * _divby(1.0 + MN_guess * MN_guess * 0.2)
        Psg = Pt * (Tsg * _divby(Tt)) ** 3.5
        if self.tight:
            return self._static_area_tight(W, Tt, Pt, fa, Ath, S, ht, Psg)
        Tsg, rhosg, Vg = self._pcalc(S, Psg, Tt, Pt, ht, fa)
        Acalc = W * _divby(Vg * rhosg / C_PSItoPSF)
        MNg = Vg * _divby(_sqrt(C_GAMMA * RT_GAS * Tsg * C_GRAVITY * JOULES_CONST))
        erA = (Ath - Acalc) * _divby(Ath)
        it = 0
        Psg_new = Psg + 0.05
        while abs(erA) > 1e-4 and it < 1000:
            erA_old, Psg_old = erA, Psg
            Psg = Psg + 0.0005 if abs(Psg - Psg_new) < 0.0003 else Psg_new
            Tsg, rhosg, Vg = self._pcalc(S, Psg, Tt, Pt, ht, fa)
            MNg = Vg * _divby(_sqrt(C_GAMMA * RT_GAS * Tsg * C_GRAVITY * JOULES_CONST))
            if Vg > 0.0001:
                erA = (Ath - W * _divby(Vg * rhosg / C_PSItoPSF)) * _divby(Ath)
            else:
                erA, Psg, Tsg = 0.0, Pt, Tt
            if abs(erA) > 1e-4:
                Psg_new = Psg - erA * (Psg - Psg_old) * _divby(erA - erA_old)
                if Psg_new > 1.001 * Psg:
                    Psg_new = 1.002 * Psg
                elif Psg_new < 0.999 * Psg:
                    Psg_new = 0.998 * Psg
            it += 1
        return Tsg, Psg, MNg

    def _static_area_tight(self, W, Tt, Pt, fa, Ath, S, ht, Ps0):
        """Subsonic root of  Ath - W/(rho V) = 0  in Ps, by a safeguarded secant iteration."""
        def err(Ps):
            Ts, rhos, V = self._pcalc(S, Ps, Tt, Pt, ht, fa)
            if V <= 1e-4:
                return -1e30, Ts, V
            return (Ath - W * C_PSItoPSF / (V * rhos)) / Ath, Ts, V
        lo, hi = 0.54 * Pt, Pt * (1.0 - 1e-13)           # subsonic branch: err(lo) > 0 > err(hi)
        f_lo = err(lo)[0]
        if f_lo <= 0.0:                                   # duct would be choked: return the sonic state
            Ts, _r, V = self._pcalc(S, lo, Tt, Pt, ht, fa)
            return Ts, lo, V * _divby(_sqrt(C_GAMMA * RT_GAS * Ts * C_GRAVITY * JOULES_CONST))
        a = min(max(Ps0, lo), hi)
        fa_, Ts, V = err(a)
        b = min(a * 1.01, hi) if fa_ > 0 else max(a * 0.99, lo)
        fb, Ts, V = err(b)
        best = (abs(fb), b, Ts, V)
        for _ in range(60):
            if fa_ > 0:
                lo = max(lo, a)
            else:
                hi = min(hi, a)
            if fb > 0:
                lo = max(lo, b)
            else:
                hi = min(hi, b)
            c = b - fb * (b - a) / (fb - fa_) if fb != fa_ else 0.5 * (lo + hi)
            if not (lo < c < hi):
                c = 0.5 * (lo + hi)
            fc, Ts, V = err(c)
            if abs(fc) < best[0]:
                best = (abs(fc), c, Ts, V)
            a, fa_, b, fb = b, fb, c, fc
            if abs(fc) < 2e-15 or abs(b - a) <= 4e-16 * abs(b) or hi - lo <= 4e-16 * hi:
                break
        _e, Ps, Ts, V = best
        return Ts, Ps, V * _divby(_sqrt(C_GAMMA * RT_GAS * Ts * C_GRAVITY * JOULES_CONST))

    def _duct(self, W, Tt, Pt, fa, key: str) -> float:
        """Duct pressure loss: Pt_out = Pt_in * (1 - dP * (MN/MNdes)^2), MN from the static calculation."""
        dP, MNdes, Ath = self.duct[key]
        MN = self._static_area(W, Tt, Pt, fa, Ath, 0.45)[2]
        return Pt * (1.0 - dP * (MN / MNdes) ** 2)

    def _nozzle(self, W, Tt, Pt, fa, Pamb, Ath, Cfg):
        """Convergent nozzle (port of Nozzle_TMATS_body.c, Cfg thrust method): returns (Fg, NErr)."""
        ht = th.t2h(Tt, fa)
        Ptin = Pt if Pt > Pamb else Pamb + 0.1
        # PcalcStat() recomputes the entropy from the pressure it is given, i.e. from the
        # substituted Ptin (found at Gate 1: the blind replica disagreed here when Pt <= Pamb)
        S = th.pt2s(Ptin, Tt, fa)
        Ts_s, rhos_s, V_s = self._pcalc(S, Pamb, Tt, Ptin, ht, fa)
        # throat static pressure at Mach 1
        TsM = Tt / 1.2
        PsM = Ptin * (TsM / Tt) ** 3.5
        kk = C_GAMMA * RT_GAS * C_GRAVITY * JOULES_CONST
        if self.tight:
            def er_of(Ps):
                Ts, rho, V = self._pcalc(S, Ps, Tt, Ptin, ht, fa)
                return 1.0 - V * _divby(_sqrt(kk * Ts)), Ts, rho
            a = PsM
            ea, TsM, rhoM = er_of(a)
            b = a * 1.001
            eb, TsM, rhoM = er_of(b)
            for _ in range(40):
                if eb == ea:
                    break
                c = b - eb * (b - a) / (eb - ea)
                a, ea = b, eb
                b = c
                eb, TsM, rhoM = er_of(b)
                if abs(eb) < 4e-16 or abs(b - a) <= 4e-16 * abs(b):
                    break
            PsM = b
        else:
            TsM, rhoM, Vg = self._pcalc(S, PsM, Tt, Ptin, ht, fa)
            er = 1.0 - Vg * _divby(_sqrt(kk * TsM))
            Ps_new = PsM + 0.05
            it = 0
            while abs(er) > 0.001 and it < 200:
                er_old, Ps_old = er, PsM
                PsM = PsM + 0.005 if abs(PsM - Ps_new) < 0.003 else Ps_new
                TsM, rhoM, Vg = self._pcalc(S, PsM, Tt, Ptin, ht, fa)
                er = 1.0 - Vg * _divby(_sqrt(kk * TsM))
                if abs(er) > 0.001:
                    Ps_new = PsM - er * (PsM - Ps_old) * _divby(er - er_old)
                it += 1
        if PsM < Pamb:                       # not choked: throat is at ambient pressure
            Vth, rhoth = V_s, rhos_s
        else:                                # choked: Mach 1 at the throat
            Vth, rhoth = _sqrt(kk * TsM), rhoM
        Wcalc = Ath * rhoth * Vth / C_PSItoPSF
        Fg = (W / C_GRAVITY) * V_s * Cfg
        NErr = 100.0 if W == 0 else (W - Wcalc) * _divby(W)
        return Fg, NErr

    # ------------------------------------------------------------------
    # the whole engine for given z, x, u
    # ------------------------------------------------------------------
    def _eval(self, z, x, ctrl: Controls, env: Env, theta, want_out: bool):
        """Run the gas path once. Returns (r list of 7, Ndot_LP, Ndot_HP, out dict or None)."""
        W, HPT_PR, HPC_RL, FAN_RL, BPR, LPC_RL, LPT_PR = z
        N_LP, N_HP = x
        t = theta or {}
        g = t.get
        Tt0, Pt0, Pamb, _Ts0, Veng = self.ambient(env)
        # Inlet: ram recovery
        Pt2 = Pt0 * self.eRam * _interp1_extrap(self.eRamVec, self.eRamTBL, Pt0 / Pamb)
        Tt2 = Tt0
        # Fan (on its own shaft, N_fan = N_LP / gear ratio)
        N_fan = N_LP / self.GR
        W21, Tt21, Pt21, Trq_fan, e_fan, SM_fan, _b, i_fan = self._compressor(
            W, Tt2, Pt2, 0.0, N_fan, FAN_RL, self.fan, g("fan_flow", 1.0), g("fan_eff", 1.0), None)
        # Splitter
        bpr = BPR if BPR > 0.0 else 0.0
        W_core = W21 / (bpr + 1.0)
        W_byp = bpr * W_core
        # duct2 -> LPC
        Pt23 = self._duct(W_core, Tt21, Pt21, 0.0, "2")
        W24, Tt24, Pt24, Trq_lpc, e_lpc, SM_lpc, _b, i_lpc = self._compressor(
            W_core, Tt21, Pt23, 0.0, N_LP, LPC_RL, self.lpc, g("lpc_flow", 1.0), g("lpc_eff", 1.0), None)
        # Variable bleed valve: LPC exit -> bypass stream (port of Valve_TMATS_body.c)
        frac = (ctrl.VBV - self.vbv_closed) * _divby(self.vbv_open - self.vbv_closed)
        vpr = Pt24 * _divby(Pt21)
        if frac <= 0.0 or vpr <= 1.0:
            W_vbv = 0.0
        else:
            W_vbv = _interp1(self.vbv_PR, self.vbv_flow, vpr) * Pt24 * _divby(_sqrt(Tt24)) * frac * self.vbv_Ae
            if W_vbv > W24:
                W_vbv = W24
        W_byp2 = W_byp + W_vbv            # bleed joins the bypass stream (no enthalpy mixing in T-MATS Valve)
        W25 = W24 - W_vbv
        # duct25 -> HPC (three fractional bleeds: LPT exit, HPT exit, HPT inlet)
        Pt25 = self._duct(W25, Tt24, Pt24, 0.0, "25")
        W3, Tt3, Pt3, Trq_hpc, e_hpc, SM_hpc, blds, i_hpc = self._compressor(
            W25, Tt24, Pt25, 0.0, N_HP, HPC_RL, self.hpc, g("hpc_flow", 1.0), g("hpc_eff", 1.0),
            list(zip(self.hpc_FBD, self.hpc_BFht, self.hpc_BFPt)))
        # Burner (port of Burner_TMATS_body.c, LHV method)
        Wf = ctrl.Wf
        W4 = W3 + Wf
        FAR4 = Wf * _divby(W3)
        Tt4 = self._h2t((W3 * th.t2h(Tt3, 0.0) + Wf * self.burn_LHV * self.burn_eff) * _divby(W4), FAR4)
        Pt4 = (1.0 - self.burn_dP) * Pt3
        # HPT (cooling: [HPT exit bleed, HPT inlet bleed])
        W45, Tt45, Pt45, FAR45, Trq_hpt, e_hpt, i_hpt = self._turbine(
            W4, Tt4, Pt4, FAR4, N_HP, HPT_PR, self.hpt, g("hpt_flow", 1.0), g("hpt_eff", 1.0), [blds[1], blds[2]])
        # duct45 -> LPT (cooling: LPT exit bleed)
        Pt48 = self._duct(W45, Tt45, Pt45, FAR45, "45")
        W5, Tt5, Pt5, FAR5, Trq_lpt, e_lpt, i_lpt = self._turbine(
            W45, Tt45, Pt48, FAR45, N_LP, LPT_PR, self.lpt, g("lpt_flow", 1.0), g("lpt_eff", 1.0), [blds[0]])
        # duct5 -> core nozzle ; duct17 -> bypass nozzle (throat area = VAFN)
        Pt7 = self._duct(W5, Tt5, Pt5, FAR5, "5")
        Fg_core, e_cn = self._nozzle(W5, Tt5, Pt7, FAR5, Pamb, self.noz_core_Ath, self.noz_core_Cfg)
        Pt17 = self._duct(W_byp2, Tt21, Pt21, 0.0, "17")
        Fg_byp, e_bn = self._nozzle(W_byp2, Tt21, Pt17, 0.0, Pamb, ctrl.VAFN, self.noz_byp_Cfg)
        # Shafts (T-MATS Shaft block): Ndot = (sum of torques) * 60 / (2 pi I)
        Trq_fan_LP = Trq_fan / self.gb_eff / self.GR         # fan torque seen by the LP shaft through the gearbox
        Trq_LP = Trq_lpc + Trq_fan_LP + Trq_lpt * self.lps_eff
        Trq_HP = Trq_hpc + Trq_hpt + SHAFT_PWR2TRQ * HP_POWER_DEBIT / N_HP
        Ndot_LP = Trq_LP * self.k_LP
        Ndot_HP = Trq_HP * self.k_HP
        r = [e_fan, e_lpc, e_hpc, e_hpt, e_lpt, e_cn, e_bn]
        if not want_out:
            return r, Ndot_LP, Ndot_HP, None
        Ps3 = self._static_area(W3, Tt3, Pt3, 0.0, self.hpc_Ath, 0.3)[1]
        Fdrag = W * Veng / G_AMBIENT
        Fg = Fg_byp + Fg_core
        Fnet = Fg - Fdrag
        out = {
            "Trq_LP": Trq_LP, "Trq_HP": Trq_HP, "Ndot_LP": Ndot_LP, "Ndot_HP": Ndot_HP,
            "Fnet": Fnet, "Fgross": Fg, "Fdrag": Fdrag, "W_in": W, "TSFC": Wf * 3600.0 * _divby(Fnet),
            "BPR": BPR, "OPR": Pt3 / Pt2, "T4": Tt4, "T45": Tt45, "Ps3": Ps3,
            "Pt2": Pt2, "Tt2": Tt2, "Pt21": Pt21, "Tt21": Tt21, "Pt25": Pt25, "Tt25": Tt24,
            "Pt3": Pt3, "Tt3": Tt3, "Pt4": Pt4, "Pt45": Pt45, "Pt5": Pt5, "Tt5": Tt5, "Pt7": Pt7, "Pt17": Pt17,
            "Pamb": Pamb, "N_fan": N_fan, "Nc_fan": N_fan / math.sqrt(Tt2 / C_TSTD),
            "SM_fan": SM_fan, "SM_lpc": SM_lpc, "SM_hpc": SM_hpc,
            "W25": W25, "W_byp": W_byp2, "W_vbv": W_vbv, "W_core": W_core, "W3": W3, "W4": W4, "W5": W5,
            "FAR4": FAR4, "Wf": Wf, "VAFN": ctrl.VAFN, "VBV": ctrl.VBV,
            "Trq_fan": Trq_fan, "Trq_fan_LP": Trq_fan_LP, "Trq_lpc": Trq_lpc, "Trq_lpt": Trq_lpt,
            "Trq_hpc": Trq_hpc, "Trq_hpt": Trq_hpt,
            "PR_fan": i_fan[1], "PR_lpc": i_lpc[1], "PR_hpc": i_hpc[1],
            "Eff_fan": i_fan[2], "Eff_lpc": i_lpc[2], "Eff_hpc": i_hpc[2], "Eff_hpt": i_hpt[1], "Eff_lpt": i_lpt[1],
            "Nc_lpc": i_lpc[0], "Nc_hpc": i_hpc[0], "Fg_byp": Fg_byp, "Fg_core": Fg_core,
        }
        return r, Ndot_LP, Ndot_HP, out

    # ------------------------------------------------------------------
    # public interface (CONTRACT.md section 4)
    # ------------------------------------------------------------------
    def residuals(self, z, x, ctrl: Controls, env: Env, theta: dict | None = None):
        """The 7 normalized flow-matching residuals r(z; x, u) and a dict of engine outputs.

        r order: [FAN, LPC, HPC, HPT, LPT, core nozzle, bypass nozzle], each = (W_in - W_map)/W_in
        in corrected-flow form, exactly as the T-MATS blocks define their "NErr" outputs.
        Course topic: Engine model (not a course topic)
        """
        r, _a, _b, out = self._eval([float(v) for v in z], [float(x[0]), float(x[1])], ctrl, env, theta, True)
        return np.array(r), out

    def flow_jacobian(self, z, x, ctrl: Controls, env: Env, theta: dict | None = None,
                      rel_step: float = 1e-7) -> np.ndarray:
        """7x7 Jacobian dr/dz by forward differences (the inner solver's own Jacobian).

        Course topic: Jacobians / partial derivatives
        """
        z = [float(v) for v in z]
        xx = [float(x[0]), float(x[1])]
        r0 = self._eval(z, xx, ctrl, env, theta, False)[0]
        J = np.empty((7, 7))
        for j in range(7):
            dz = rel_step * max(abs(z[j]), 1e-3)
            zp = list(z)
            zp[j] += dz
            rp = self._eval(zp, xx, ctrl, env, theta, False)[0]
            for i in range(7):
                J[i, j] = (rp[i] - r0[i]) / dz
        return J

    def solve_flow(self, x, ctrl: Controls, env: Env, theta: dict | None = None, tol: float = 1e-12,
                   z0=None, max_iter: int = 60, jac: np.ndarray | None = None):
        """Inner Newton-Raphson (chord) solve of r(z) = 0 for fixed shaft speeds x.

        The ONLY stopping rule is max|r_i| < tol, so `tol` is exactly the residual left behind.
        `jac`: optional precomputed 7x7 dr/dz that is re-used (chord iterations, as T-MATS does);
        it is refreshed by finite differences when an iteration fails to halve the residual.
        Returns (z, out, info) with info = {'converged', 'iters', 'resnorm', 'nfev'}.
        Course topic: Engine model (not a course topic) - Newton's method for nonlinear systems
        """
        xx = [float(x[0]), float(x[1])]
        z = [float(v) for v in (z0 if z0 is not None else DEFAULT_Z(env))]
        nfev = 1
        r = self._eval(z, xx, ctrl, env, theta, False)[0]
        res = max(abs(v) for v in r)
        J_inv = None
        if jac is not None:
            J_inv = np.linalg.inv(np.asarray(jac, dtype=float))
        iters = 0
        fresh = False
        stalls = 0
        while res >= tol and iters < max_iter:
            if J_inv is None:
                J_inv = np.linalg.inv(self.flow_jacobian(z, xx, ctrl, env, theta))
                nfev += 8
                fresh = True
            dz = -(J_inv @ np.array(r))
            lam = 1.0
            accepted = False
            for _ in range(6):
                zt = [z[i] + lam * float(dz[i]) for i in range(7)]
                if zt[0] > 0 and zt[1] > 1.0 and zt[6] > 1.0 and zt[4] > 0:
                    rt = self._eval(zt, xx, ctrl, env, theta, False)[0]
                    nfev += 1
                    rest = max(abs(v) for v in rt)
                    if rest == rest and rest < res * (1.0 if not fresh else 2.0) or (rest == rest and rest < 0.9 * res):
                        accepted = True
                        break
                lam *= 0.5
            iters += 1
            if not accepted:
                if fresh:                      # even a fresh Jacobian cannot reduce the residual
                    stalls += 1
                    if stalls >= 2:
                        break
                J_inv = None                   # stale Jacobian: recompute and try again
                continue
            if rest > 0.5 * res and not fresh:
                J_inv = None                   # slow convergence: refresh the Jacobian next time
            fresh = False if J_inv is not None else fresh
            z, r, res = zt, rt, rest
            if J_inv is not None:
                fresh = False
        _r, _a, _b, out = self._eval(z, xx, ctrl, env, theta, True)
        info = {"converged": bool(res < tol), "iters": iters, "resnorm": float(res), "nfev": nfev + 1}
        return np.array(z), out, info

    def F(self, x, ctrl: Controls, env: Env, theta: dict | None = None, tol: float = 1e-12,
          z0=None, jac: np.ndarray | None = None) -> np.ndarray:
        """State derivative F(x, u; theta) = [dN_LP/dt, dN_HP/dt] in rpm/s.

        Solves the flow-matching equations to tolerance `tol`, then applies the shaft equations
        I (pi/30) dN/dt = net torque.
        Course topic: Linear ODE systems (this is the nonlinear right-hand side that gets linearized)
        """
        _z, out, _info = self.solve_flow(x, ctrl, env, theta, tol=tol, z0=z0, jac=jac)
        return np.array([out["Ndot_LP"], out["Ndot_HP"]])

    def schedule_VAFN(self, Nc_fan: float, MN: float) -> float:
        """NASA baseline VAFN schedule MWS.Cntrl.VAFN_sch(Nc_fan, MN) [in^2], clamped at the table edges."""
        return self._vafn(MN, Nc_fan)

    def schedule_VBV(self, Nc_fan: float, MN: float) -> float:
        """NASA baseline VBV schedule MWS.Cntrl.VBV_sch(Nc_fan, MN) [0 = closed .. 1 = open]."""
        return self._vbv(MN, Nc_fan)

    def Nc_fan_of(self, N_LP: float, env: Env) -> float:
        """Corrected fan speed N_fan / sqrt(Tt2 / 518.67) for LP-shaft speed N_LP at this flight condition."""
        return N_LP / self.GR / math.sqrt(self.ambient(env)[0] / C_TSTD)


# NASA's own initial guesses (setup_InputsSS.m), used when the caller gives no starting point
Z_CRUISE = (815.48, 4.0955, 1.994, 2.0007, 23.19937, 2.19956, 11.0178)
Z_TAKEOFF = (1723.81, 4.17, 2.0325, 1.75, 27.5078, 1.8740, 7.15)


def DEFAULT_Z(env: Env) -> tuple:
    """Starting guess for z: NASA's cruise or take-off initial conditions, whichever is nearer in altitude."""
    return Z_CRUISE if env.alt_ft > 17500.0 else Z_TAKEOFF


_ENGINES: dict = {}
MODES = {"study": ("tight", "smooth"), "nasa": ("tmats", "linear")}


def load_engine(vendor_dir: str | Path | None = None, mode: str = "study",
                iter_mode: str | None = None, map_interp: str | None = None) -> Engine:
    """Build the engine from the vendored NASA files (real AGTF30 maps). Cached per configuration.

    mode='nasa'  : the faithful port - T-MATS inner loops with NASA's thresholds and bilinear map
                   interpolation. This is what is validated against NASA's published numbers and
                   against the blind replica.
    mode='study' : (default) the same engine made differentiable - inner loops converged to machine
                   precision and C1 bicubic map interpolation through the same map points - so that
                   Jacobians are well defined and `tol` is the only solver noise (DERIVED deviation).
    `iter_mode` ('tight'|'tmats') and `map_interp` ('smooth'|'linear') override the two switches one by one.
    Course topic: Engine model (not a course topic)
    """
    im, mi = MODES[mode]
    im, mi = iter_mode or im, map_interp or mi
    key = (str(vendor_dir) if vendor_dir else None, im, mi)
    if key not in _ENGINES:
        _ENGINES[key] = Engine(load_params(vendor_dir), im, mi)
    return _ENGINES[key]
