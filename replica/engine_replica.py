"""A2 blind replica of NASA's AGTF30 geared turbofan (T-MATS), written from the NASA source only.

Every equation below is a transcription of a T-MATS C routine
(vendor/tmats/TMATS_Library/MEX/*_body.c) or of a pure-Simulink library block
(Lib_Turbo_*.slx / .mdl); the wiring is the one of AGTF30_eng.mdl (checked at load time by
replica._sources.check_wiring). The internal T-MATS iterations (h2tc, sp2tc, Ambient,
StaticCalc, Nozzle) use the SAME thresholds, iteration caps and step rules as the C code.

Units: lbm/s, psia, deg R, BTU/lbm, lbf, ft lbf, in^2, rpm, hp.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Callable, Sequence

import numpy as np

from ._sources import Sources

THETA_KEYS = ("fan_eff", "fan_flow", "lpc_eff", "lpc_flow", "hpc_eff", "hpc_flow",
              "hpt_eff", "hpt_flow", "lpt_eff", "lpt_flow")
_NUMERIC_ERRORS = (ValueError, ZeroDivisionError, OverflowError, FloatingPointError)

# NASA steady-state initial conditions (setup_InputsSS.m, cases 'cruise' and 'takeoff'):
# used ONLY as Newton starting guesses.
NASA_IC = {
    "cruise": dict(alt=35000.0, MN=0.8, dT=0.0, Wf=0.7837, N2=6772.0, N3=20796.0, VAFN=4775.22, VBV=0.0,
                   Nc_frac=1.0, z=(815.48, 4.0955, 1.994, 2.0007, 23.19937, 2.19956, 11.0178)),
    "takeoff": dict(alt=0.0, MN=0.0, dT=27.0, Wf=1.3896, N2=6079.0, N3=21583.2, VAFN=6314.82, VBV=0.0,
                    Nc_frac=0.810, z=(1723.81, 4.17, 2.0325, 1.75, 27.5078, 1.8740, 7.15)),
}


@dataclass(frozen=True)
class Env:
    """Flight condition: altitude [ft], Mach number, delta ambient temperature [deg F = deg R]."""
    alt_ft: float
    MN: float
    dT: float = 0.0


@dataclass(frozen=True)
class Controls:
    """Fuel flow [lbm/s], bypass-nozzle throat area [in^2], VBV position (0 closed .. 1 open)."""
    Wf: float
    VAFN: float
    VBV: float


# ----------------------------------------------------------------------------- T-MATS helpers
def sqrtT(x: float) -> float:
    """functions_TMATS.c sqrtT: square root, 0 for negative input."""
    return 0.0 if x < 0 else math.sqrt(x)


def divby(x: float) -> float:
    """functions_TMATS.c divby: 1/x with the +-1e10 guard near zero."""
    if -1e-10 < x < 1e-10:
        return 1e10 if x >= 0 else -1e10
    return 1.0 / x


def powT(a: float, n: float) -> float:
    """functions_TMATS.c powT: a^n with the 1e10 guard."""
    if -1e-10 < a < 1e-10 and n < 0:
        return 1e10
    return math.pow(a, n)


def interp1Ac(X: Sequence[float], Y: Sequence[float], xi: float) -> float:
    """interp1Ac_TMATS.c: linear interpolation, clamped at both table ends."""
    n = len(X)
    if xi < X[0] or not (xi >= X[0]):
        xi = X[0]
    elif xi > X[n - 1]:
        xi = X[n - 1]
    i = n - 2
    while i > 0 and xi < X[i]:
        i -= 1
    return Y[i] + (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]) * (xi - X[i])


def interp2Ac(X: Sequence[float], Y: Sequence[float], Z: Sequence[Sequence[float]],
              xi: float, yi: float) -> float:
    """interp2Ac_TMATS.c: Z[row = Y index][col = X index]; first along X, then along Y; clamped."""
    a, b = len(X), len(Y)
    if xi < X[0] or not (xi >= X[0]):
        xi = X[0]
    elif xi > X[a - 1]:
        xi = X[a - 1]
    if yi < Y[0] or not (yi >= Y[0]):
        yi = Y[0]
    elif yi > Y[b - 1]:
        yi = Y[b - 1]
    ii = a - 2
    while ii > 0 and xi < X[ii]:
        ii -= 1
    jj = b - 2
    while jj > 0 and yi < Y[jj]:
        jj -= 1
    r0, r1 = Z[jj], Z[jj + 1]
    dx = X[ii + 1] - X[ii]
    z1 = r0[ii] + (r0[ii + 1] - r0[ii]) / dx * (xi - X[ii])
    z2 = r1[ii] + (r1[ii + 1] - r1[ii]) / dx * (xi - X[ii])
    return z1 + (z2 - z1) / (Y[jj + 1] - Y[jj]) * (yi - Y[jj])


def lookup_clip(X: Sequence[float], Y: Sequence[float], xi: float) -> float:
    """Simulink 1-D Lookup Table, linear interpolation, ExtrapMethod = Clip."""
    return interp1Ac(X, Y, xi)


def make_thermo(tab: dict[str, list[float]]) -> tuple[Callable, Callable, Callable, Callable]:
    """Build t2hc, h2tc, pt2sc, sp2tc from the tables embedded in the T-MATS C files."""
    ITAB = [int(v) for v in tab["ITAB"]]
    TT = tab["TTTAB"]
    AH, BH, CH, DH = tab["AHAIR"], tab["BHAIR"], tab["CHAIR"], tab["DHAIR"]
    AS, BS, CS, DS = tab["AHSTOC"], tab["BHSTOC"], tab["CHSTOC"], tab["DHSTOC"]
    AP, BP, CP, DP = tab["APAIR"], tab["BPAIR"], tab["CPAIR"], tab["DPAIR"]
    AQ, BQ, CQ, DQ = tab["APSTOC"], tab["BPSTOC"], tab["CPSTOC"], tab["DPSTOC"]
    fmod, log = math.fmod, math.log

    def mix(fa: float) -> tuple[float, float, float, float]:
        if fa == 0:
            return 4.7642, 0.0, 0.2098988288, 0.0345194683
        tmls = 4.7642 + fa * 4.721362582
        return 4.7642 - fa * 69.69056873, fa * 74.411931335, 1 / tmls, tmls / (138.0148721 * (1 + fa))

    def tidx(T: float) -> int:
        temp = 0.01 * (T - fmod(T, 100))
        k = int(temp) if temp > 1 else 1
        return ITAB[(k if k < 60 else 60) - 1] - 1

    def t2hc(T: float, fa: float) -> float:
        """t2hc_TMATS.c: enthalpy [BTU/lbm] from temperature [R] and fuel-air ratio."""
        zmea, zmsp, tmlsr, zmwtr = mix(fa)
        it = tidx(T)
        dl = T - TT[it]
        hgea = ((DH[it] * dl + CH[it]) * dl + BH[it]) * dl + AH[it]
        hgsp = ((DS[it] * dl + CS[it]) * dl + BS[it]) * dl + AS[it] if fa > 0 else 0.0
        return (hgsp * zmsp + hgea * zmea) * (tmlsr * zmwtr)

    def h2tc(H: float, fa: float) -> float:
        """h2tc_TMATS.c: temperature from enthalpy; secant, |dh| > 1e-3, at most 11 passes."""
        zmea, zmsp, tmlsr, zmwtr = mix(fa)
        tg, tgo, hgo = 3.55 * H, 1.0, 0.282
        zz = tmlsr * zmwtr
        ii, hh = 0, 99.99
        while abs(hh) > 1e-3 and ii < 11:
            ii += 1
            it = tidx(tg)
            d1 = tg - TT[it]
            hgea = ((DH[it] * d1 + CH[it]) * d1 + BH[it]) * d1 + AH[it]
            hgsp = ((DS[it] * d1 + CS[it]) * d1 + BS[it]) * d1 + AS[it] if fa > 0 else 0.0
            hg = (hgsp * zmsp + hgea * zmea) * zz
            hh = H - hg
            d = (tg - tgo) / (hg - hgo)
            tgo, hgo = tg, hg
            tg = tg + d * hh
        return tg

    def phi(T: float, fa: float, zmea: float, zmsp: float, tmlsr: float) -> float:
        it = tidx(T)
        dl = T - TT[it]
        phiea = ((DP[it] * dl + CP[it]) * dl + BP[it]) * dl + AP[it]
        phisp = ((DQ[it] * dl + CQ[it]) * dl + BQ[it]) * dl + AQ[it] if fa > 0 else 0.0
        return (phisp * zmsp + phiea * zmea) * tmlsr

    def pt2sc(P: float, T: float, fa: float) -> float:
        """pt2sc_TMATS.c: entropy from pressure [psia] and temperature [R]."""
        zmea, zmsp, tmlsr, zmwtr = mix(fa)
        rcas = 1.98587 * zmwtr
        return (phi(T, fa, zmea, zmsp, tmlsr) * 0.5035576347 - 23.0258509) * rcas - 0.1841304 \
            - rcas * log(P / 14.696)

    def sp2tc(S: float, P: float, fa: float) -> float:
        """sp2tc_TMATS.c: temperature from entropy and pressure; secant, |dS| >= 1e-4, 2..10 passes."""
        zmea, zmsp, tmlsr, zmwtr = mix(fa)
        rcas = 1.98587 * zmwtr
        lp = rcas * log(P / 14.696)
        Tg, jj, Sg = 1000.0, 0, 1e3
        Sg1 = Tg1 = 0.0
        while (abs(S - Sg) >= 1e-4 and jj < 10) or jj < 2:
            Sg = (phi(Tg, fa, zmea, zmsp, tmlsr) * 0.5035576347 - 23.0258509) * rcas - 0.1841304 - lp
            if jj == 0:
                Sg1, Tg1 = Sg, Tg
                Tg = Tg1 + 50
            else:
                dTdS = (Tg - Tg1) / (Sg - Sg1)
                Tg1, Sg1 = Tg, Sg
                Tg = Tg1 + (S - Sg1) * dTdS
            jj += 1
        return Tg

    return t2hc, h2tc, pt2sc, sp2tc


class _CompMap:
    """Compressor map + scalars + bleed set-up of one Compressor block."""

    def __init__(self, src: Sources, block: str):
        g = lambda k: src.vec(block, k)
        self.R, self.Nc = g("X_C_RlineVec_M"), g("Y_C_Map_NcVec_M")
        arr = lambda k: [list(map(float, row)) for row in np.atleast_2d(src.mask(block, k))]
        self.Wc, self.PR, self.Eff = arr("T_C_Map_WcArray_M"), arr("T_C_Map_PRArray_M"), arr("T_C_Map_EffArray_M")
        self.s_Nc, self.s_Wc = src.mask(block, "s_C_Nc_M"), src.mask(block, "s_C_Wc_M")
        self.s_PR, self.s_Eff = src.mask(block, "s_C_PR_M"), src.mask(block, "s_C_Eff_M")
        self.SMWc, self.SMPR = g("X_C_Map_WcSurgeVec_M"), g("T_C_Map_PRSurgeVec_M")
        self.fbld_en = src.mask(block, "FBLDEN_M") > 0.5
        self.FBD = g("C_FBD_M") if self.fbld_en else []
        self.BFht = g("C_BFht_M") if self.fbld_en else []
        self.BFPt = g("C_BFPt_M") if self.fbld_en else []
        assert len(g("Z_C_AlphaVec_M")) == 1 and src.mask(block, "iDesign_M") == 2.0
        assert src.mask(block, "SMNEn_M") < 0.5 and src.mask(block, "CBLDEN_M") < 0.5
        assert len(self.Wc) == len(self.Nc) and len(self.Wc[0]) == len(self.R)


class _TurbMap:
    """Turbine map + scalars + cooling-flow positions of one 'Turbine, NPSS' block."""

    def __init__(self, src: Sources, block: str):
        g = lambda k: src.vec(block, k)
        self.PRv, self.Nc = g("X_T_Map_PRVec_M"), g("Y_T_Map_NcVec_M")
        arr = lambda k: [list(map(float, row)) for row in np.atleast_2d(src.mask(block, k))]
        self.Wc, self.Eff = arr("T_T_Map_WcArray_M"), arr("T_T_Map_EffArray_M")
        self.s_Nc, self.s_Wc = src.mask(block, "s_T_Nc_M"), src.mask(block, "s_T_Wc_M")
        self.s_PR, self.s_Eff = src.mask(block, "s_T_PR_M"), src.mask(block, "s_T_Eff_M")
        self.pos = g("T_CoolFlwPos_M")
        assert src.mask(block, "ConfigNPSS_M") > 0.5 and src.mask(block, "CFlwEn_M") > 0.5
        assert src.mask(block, "iDesign_M") == 2.0
        assert len(self.Wc) == len(self.Nc) and len(self.Wc[0]) == len(self.PRv)


class Engine:
    """AGTF30 engine: 7 flow-matching residuals, shaft torques, steady flow solve."""

    def __init__(self, src: Sources):
        if src.wiring_missing:
            raise RuntimeError("AGTF30_eng.mdl wiring differs from the replica: " + "; ".join(src.wiring_missing))
        self.src = src
        M, c = src.MWS, src.const
        self.t2hc, self.h2tc, self.pt2sc, self.sp2tc = make_thermo(src.tables)
        self.PSTD, self.TSTD = c["C_PSTD"], c["C_TSTD"]
        self.HP, self.TRQ = c["C_BTU_PER_SECtoHP"], c["C_HP_PER_RPMtoFT_LBF"]
        self.G, self.J, self.GAMMA = c["C_GRAVITY"], c["JOULES_CONST"], c["C_GAMMA"]
        self.PSF, self.SFT = c["C_PSItoPSF"], c["C_SINtoSFT"]
        # --- Ambient block
        A = "Engine/Ambient"
        self.amb_alt, self.amb_Ts, self.amb_Ps = src.vec(A, "X_A_AltVec_M"), src.vec(A, "T_A_TsVec_M"), src.vec(A, "T_A_PsVec_M")
        self.amb_FAR = src.mask(A, "AFARc_M")
        self.Rt = self._uniform(A, "T_A_RtArray_M")          # gas constant table (all 0.0686)
        self.gam = self._uniform(A, "T_A_gammaArray_M")      # gamma table (all 1.4)
        for blk, ks in (("Engine/Bypass/BypNoz", ("T_N_RtArray_M", "T_N_MAP_gammaArray_M")),
                        ("Engine/BackEng/CoreNoz", ("T_N_RtArray_M", "T_N_MAP_gammaArray_M")),
                        ("Engine/FrontEng/duct2/StaticCalc", ("T_RtArray_M", "T_gammaArray_M"))):
            assert self._uniform(blk, ks[0]) == self.Rt and self._uniform(blk, ks[1]) == self.gam
        self.amb_g = src.lib["ambient_g"]
        # --- Inlet
        I = "Engine/FrontEng/Inlet"
        self.eRam, self.eRamVec, self.eRamTbl = src.mask(I, "eRambase_M"), src.vec(I, "X_eRamVec_M"), src.vec(I, "T_eRamtbl_M")
        # --- turbomachinery
        self.fan, self.lpc = _CompMap(src, "Engine/FrontEng/FAN"), _CompMap(src, "Engine/FrontEng/LPC")
        self.hpc = _CompMap(src, "Engine/Core/HPC")
        self.hpt, self.lpt = _TurbMap(src, "Engine/Core/HPT"), _TurbMap(src, "Engine/BackEng/LPT")
        self.GR = float(src.mask("Engine/FrontEng/GearBox", "R_M"))
        self.gear_eff = float(src.mask("Engine/FrontEng/GearBox", "Eff_M"))
        # --- VBV valve
        V = "Engine/FrontEng/Valve"
        self.vlv_open, self.vlv_dead, self.vlv_Ae = src.mask(V, "VFO_M"), src.mask(V, "VDZ_M"), src.mask(V, "VAe_M")
        self.vlv_PR, self.vlv_W = src.vec(V, "X_VPRVec_M"), src.vec(V, "T_VWArray_M")
        # --- ducts (masked subsystems of the engine model) and S36 static calc
        self.ducts = {}
        for name, path in (("2", "Engine/FrontEng/duct2"), ("25", "Engine/FrontEng/duct25"),
                           ("45", "Engine/BackEng/duct45"), ("5", "Engine/BackEng/duct5"),
                           ("17", "Engine/Bypass/duct17")):
            assert src.mask(path + "/StaticCalc", "ST_M") < 0.5           # SolveType 0: area known
            self.ducts[name] = (float(src.mask(path, "dP_M")), float(src.mask(path, "MNdes_M")),
                                float(src.mask(path, "Ath_M")), float(src.mask(path + "/StaticCalc", "MN_M")),
                                float(src.mask(path + "/Constant", "Value")))
        S = "Engine/Core/S36 Static Calc/StaticCalc"
        self.s36 = (float(src.mask(S, "Ath_M")), float(src.mask(S, "MN_M")))
        # --- burner
        B = "Engine/Core/Burner"
        assert src.mask(B, "LHVEn_M") > 0.5
        self.LHV, self.burn_eff, self.burn_dP = src.mask(B, "LHV_M"), src.mask(B, "Eff_M"), src.mask(B, "dPnormBurner_M")
        # --- nozzles
        self.noz = {}
        for name, path in (("byp", "Engine/Bypass/BypNoz"), ("core", "Engine/BackEng/CoreNoz")):
            assert src.mask(path, "SWType_M") < 1.5 and src.mask(path, "CfgEn_M") > 0.5
            assert src.mask(path, "iDesign_M") == 2.0
            self.noz[name] = dict(fl=src.mask(path, "N_FL_M"), PQ=src.vec(path, "X_N_PEQPaVec_M"),
                                  Cd=src.vec(path, "T_N_CdThArray_M"), Cfg=src.vec(path, "T_N_CfgArray_M"),
                                  TGx=src.vec(path, "X_N_TtVec2_M"), TG=src.vec(path, "T_N_TGArray_M"))
        self.core_Ath = float(src.mask("Engine/BackEng/CoreNoz", "N_TArea_M"))
        # --- shafts
        self.I_LP = float(src.mask("Engine/LP_Shaft", "Shaft_Inertia_M"))
        self.I_HP = float(src.mask("Engine/HP_Shaft", "Shaft_Inertia_M"))
        self.LPS_eff = float(src.mask("Engine/Constant", "Value"))        # multiplies the LPT torque
        self.pwr_LP = float(src.mask("Engine/Constant1", "Value"))        # hp taken from the LP shaft (0)
        self.pwr_HP = float(src.mask("Engine/PwrDebit [HP]", "Value"))    # hp taken from the HP shaft (350)
        self.shaft_gain = src.lib["shaft_gain"]
        self.secs2hrs = src.lib["secs2hrs"]
        # --- NASA schedules
        C = M["Cntrl"]
        self._vafn = (list(C["VAFN_MN"]), list(C["VAFN_Nc1"]), [list(r) for r in C["VAFN_sch"]])
        self._vbv = (list(C["VBV_MN"]), list(C["VBV_Nc1"]), [list(r) for r in C["VBV_sch"]])
        self.s_Nc_fan = float(M["FAN"]["s_Nc"])
        self._amb_cache: dict[tuple[float, float, float], tuple] = {}
        self.last_info: dict = {}

    def _uniform(self, block: str, key: str) -> float:
        """Value of a table whose entries are all equal (Rt and gamma tables of this model)."""
        v = self.src.vec(block, key)
        if max(v) != min(v):
            raise NotImplementedError(f"{block}:{key} is not constant; general lookup not implemented")
        return v[0]

    # ------------------------------------------------------------------ schedules
    def schedule_VAFN(self, Nc_fan: float, MN: float) -> float:
        """NASA baseline VAFN schedule MWS.Cntrl.VAFN_sch(Nc1, MN) [in^2], clamped bilinear."""
        return interp2Ac(self._vafn[0], self._vafn[1], self._vafn[2], MN, Nc_fan)

    def schedule_VBV(self, Nc_fan: float, MN: float) -> float:
        """NASA baseline VBV schedule MWS.Cntrl.VBV_sch(Nc1, MN) [0..1], clamped bilinear."""
        return interp2Ac(self._vbv[0], self._vbv[1], self._vbv[2], MN, Nc_fan)

    # ------------------------------------------------------------------ components
    def ambient(self, env: Env) -> tuple:
        """Ambient_TMATS_body.c -> (ht, Tt, Pt, FAR, Ps, Ts, Veng). Cached per flight condition."""
        key = (env.alt_ft, env.MN, env.dT)
        hit = self._amb_cache.get(key)
        if hit is not None:
            return hit
        FAR, MN = self.amb_FAR, env.MN
        Rs, gam, GJ = self.Rt, self.gam, self.G * self.J
        Ts = interp1Ac(self.amb_alt, self.amb_Ts, env.alt_ft) + env.dT
        Ps = interp1Ac(self.amb_alt, self.amb_Ps, env.alt_ft)
        Sout = self.pt2sc(Ps, Ts, FAR)
        hs = self.t2hc(Ts, FAR)
        Ttg = Ts * (1 + MN * MN * (self.GAMMA - 1) / 2)
        Ptg = Ps * divby(powT(Ts * divby(Ttg), self.GAMMA * divby(self.GAMMA - 1)))
        Ttg = self.sp2tc(Sout, Ptg, FAR)
        htg = self.t2hc(Ttg, FAR)
        Vg = sqrtT(2 * (htg - hs) * GJ)
        Vsg = sqrtT(gam * Rs * Ts * GJ)
        er = MN - Vg * divby(Vsg)
        Ptg_new = Ptg + 0.05
        it = 0
        while abs(er) > 0.001 and it < 15:
            er_old, Ptg_old = er, Ptg
            Ptg = Ptg + 0.05 if abs(Ptg - Ptg_new) < 0.03 else Ptg_new
            Ttg = self.sp2tc(Sout, Ptg, FAR)
            htg = self.t2hc(Ttg, FAR)
            Vg = sqrtT(2 * (htg - hs) * GJ)
            er = MN - Vg * divby(Vsg)
            if abs(er) > 0.001:
                Ptg_new = Ptg - er * (Ptg - Ptg_old) * divby(er - er_old)
            it += 1
        res = (htg, Ttg, Ptg, FAR, Ps, Ts, Vsg * MN)
        self._amb_cache[key] = res
        return res

    def _pcalc(self, Pt: float, Ps: float, Tt: float, ht: float, FAR: float) -> tuple:
        """PcalcStat_TMATS.c -> (Ts, hs, rhos, V)."""
        S = self.pt2sc(Pt, Tt, FAR)
        Ts = self.sp2tc(S, Ps, FAR)
        if Ts > Tt:
            Ts = Tt
        hs = self.t2hc(Ts, FAR)
        if hs > ht:
            hs = ht
        rhos = Ps * self.PSF * divby(self.Rt * Ts * self.J)
        return Ts, hs, rhos, sqrtT(2 * (ht - hs) * self.G * self.J)

    def static_area(self, W: float, Tt: float, Pt: float, FAR: float, Ath: float, MNin: float) -> tuple:
        """StaticCalc_TMATS_body.c, SolveType 0 (area known) -> (Ts, Ps, rhos, MN, Acalc)."""
        gam, Rs, GJ = self.gam, self.Rt, self.G * self.J
        htin = self.t2hc(Tt, FAR)
        Tsg = Tt * divby(1 + MNin * MNin * (1.4 - 1) / 2)
        Psg = Pt * powT(Tsg * divby(Tt), 1.4 * divby(1.4 - 1))
        Tsg, hsg, rhosg, Vg = self._pcalc(Pt, Psg, Tt, htin, FAR)
        Acalc = W * divby(Vg * rhosg / self.SFT)
        MNg = Vg * divby(sqrtT(gam * Rs * Tsg * GJ))
        erA = (Ath - Acalc) * divby(Ath)
        it = 0
        Psg_new = Psg + 0.05
        while abs(erA) > 0.0001 and it < 1000:
            erA_old, Psg_old = erA, Psg
            Psg = Psg + 0.0005 if abs(Psg - Psg_new) < 0.0003 else Psg_new
            Tsg, hsg, rhosg, Vg = self._pcalc(Pt, Psg, Tt, htin, FAR)
            MNg = Vg * divby(sqrtT(gam * Rs * Tsg * GJ))
            if Vg > 0.0001:
                Acalc = W * divby(Vg * rhosg / self.SFT)
                erA = (Ath - Acalc) * divby(Ath)
            else:
                erA, Psg, Tsg, Acalc = 0.0, Pt, Tt, 999.0
            if abs(erA) > 0.0001:
                Psg_new = Psg - erA * (Psg - Psg_old) * divby(erA - erA_old)
                if Psg_new > 1.001 * Psg:
                    Psg_new = 1.002 * Psg
                elif Psg_new < 0.999 * Psg:
                    Psg_new = 0.998 * Psg
            it += 1
        return Tsg, Psg, rhosg, MNg, Acalc

    def duct(self, name: str, flow: tuple) -> tuple:
        """Engine-model duct subsystem: Pt_out = Pt_in*(1 - (MN/MNdes)^2*dP); returns (flow, MN)."""
        dP, MNdes, Ath, MNguess, one = self.ducts[name]
        W, ht, Tt, Pt, FAR = flow
        MN = self.static_area(W, Tt, Pt, FAR, Ath, MNguess)[3]
        q = MN / MNdes
        return (W, ht, Tt, Pt * (one - q * q * dP), FAR), MN

    def compressor(self, m: _CompMap, flow: tuple, Nmech: float, Rline: float,
                   th_flow: float, th_eff: float) -> dict:
        """Compressor_TMATS_body.c (iDesign 2, alpha fixed, constant-Wc stall margin)."""
        WIn, _, TtIn, PtIn, FAR = flow
        t2hc, h2tc = self.t2hc, self.h2tc
        htin = t2hc(TtIn, FAR)
        Sin = self.pt2sc(PtIn, TtIn, FAR)
        delta, theta = PtIn / self.PSTD, TtIn / self.TSTD
        Wcin = WIn * sqrtT(theta) * divby(delta)
        Nc = Nmech * divby(sqrtT(theta))
        NcMap = Nc * divby(m.s_Nc)
        C_Wc, C_Eff = (1 + (th_flow - 1.0)) * m.s_Wc, (1 + (th_eff - 1.0)) * m.s_Eff   # scalarMod: (1+mod)*s
        WcMap = interp2Ac(m.R, m.Nc, m.Wc, Rline, NcMap)
        WcCalc = WcMap * C_Wc
        PRMap = interp2Ac(m.R, m.Nc, m.PR, Rline, NcMap)
        PR = m.s_PR * (PRMap - 1) + 1
        EffMap = interp2Ac(m.R, m.Nc, m.Eff, Rline, NcMap)
        Eff = EffMap * C_Eff
        PtOut = PtIn * PR
        htIdeal = t2hc(self.sp2tc(Sin, PtOut, FAR), FAR)
        htOut = ((htIdeal - htin) * divby(Eff)) + htin
        TtOut = h2tc(htOut, FAR)
        Wbleeds, PwrBld, bleeds = 0.0, 0.0, []
        for i, frac in enumerate(m.FBD):
            if frac <= 0:
                bleeds.append((0.0, 0.0, 0.0, 0.0, 0.0))
                continue
            Wb = frac * WIn
            Wbleeds += Wb
            Ptb = PtIn + m.BFPt[i] * (PtOut - PtIn)
            htb = htin + m.BFht[i] * (htOut - htin)
            bleeds.append((Wb, htb, h2tc(htb, FAR), Ptb, FAR))
            PwrBld += Wb * (htb - htOut) * self.HP
        Pwrout = WIn * (htin - htOut) * self.HP - PwrBld
        SPRMap = interp1Ac(m.SMWc, m.SMPR, WcMap)
        SPR = m.s_PR * (SPRMap - 1) + 1
        return dict(flow=(WIn - Wbleeds, htOut, TtOut, PtOut, FAR), bleeds=bleeds,
                    Trq=self.TRQ * Pwrout * divby(Nmech),
                    NErr=100.0 if WIn == 0 else (Wcin - WcCalc) * divby(Wcin),
                    SM=(SPR - PR) * divby(PR) * 100, Nc=Nc, PR=PR, Wcin=Wcin, NcMap=NcMap, Eff=Eff)

    def burner(self, Wf: float, flow: tuple) -> tuple:
        """Burner_TMATS_body.c (LHV mode)."""
        WIn, _, TtIn, PtIn, FAR = flow
        htin = self.t2hc(TtIn, FAR)
        WOut = WIn + Wf
        FARo = (WIn * FAR + Wf) * divby(WIn * (1 - FAR))
        htOut = (WIn * htin + Wf * self.LHV * self.burn_eff) * divby(WOut)
        return (WOut, htOut, self.h2tc(htOut, FARo), (1 - self.burn_dP) * PtIn, FARo)

    def turbine(self, m: _TurbMap, cool: list, flow: tuple, PRIn: float, Nmech: float,
                th_flow: float, th_eff: float) -> dict:
        """Turbine_TMATS_body.c with ConfigNPSS = 1 (the 'Turbine, NPSS' block), iDesign 2."""
        WIn, _, TtIn, PtIn, FAR = flow
        t2hc = self.t2hc
        Wcools1 = Wcoolout = Wfcools1 = Wfcoolout = dHcools1 = dHcoolout = 0.0
        for (Wc, _h, Ttc, _Ptc, FARc), pos in zip(cool, m.pos):
            Wcools1 += Wc * (1 - pos)
            Wcoolout += Wc
            Wfcools1 += FARc * Wc * (1 - pos) * divby(1 + FARc)
            Wfcoolout += FARc * Wc * divby(1 + FARc)
        for (Wc, _h, Ttc, _Ptc, FARc), pos in zip(cool, m.pos):
            htc = t2hc(Ttc, FARc)
            dHcools1 += htc * Wc * (1 - pos)
            dHcoolout += htc * Wc * pos
        Ws1in, WOut = WIn + Wcools1, WIn + Wcoolout
        Wa = WIn * divby(1 + FAR)
        FARs1 = (FAR * Wa + Wfcools1) * divby(Wa + Wcools1 - Wfcools1)
        FARo = (FAR * Wa + Wfcoolout) * divby(Wa + Wcoolout - Wfcoolout)
        htin = t2hc(TtIn, FAR)
        hts1 = (htin * WIn + dHcools1) * divby(Ws1in)
        Tts1 = self.h2tc(hts1, FARs1)
        Ss1 = self.pt2sc(PtIn, Tts1, FARs1)
        Nc = Nmech * divby(sqrtT(TtIn))                      # NPSS: no standard-day correction
        NcMap = Nc * divby(m.s_Nc)
        PRmap = ((PRIn - 1) * divby(m.s_PR)) + 1
        PtOut = PtIn * divby(PRIn)
        C_Wc, C_Eff = (1 + (th_flow - 1.0)) * m.s_Wc, (1 + (th_eff - 1.0)) * m.s_Eff
        WcMap = interp2Ac(m.PRv, m.Nc, m.Wc, PRmap, NcMap)
        WcCalc = WcMap * C_Wc
        Wcin = WIn * sqrtT(TtIn) * divby(PtIn)
        Eff = interp2Ac(m.PRv, m.Nc, m.Eff, PRmap, NcMap) * C_Eff
        htIdeal = t2hc(self.sp2tc(Ss1, PtOut, FARo), FARo)
        Pwrout = ((hts1 - htIdeal) * Eff) * Ws1in * self.HP
        htOut = ((((htIdeal - hts1) * Eff) + hts1) * Ws1in + dHcoolout) * divby(WOut)
        return dict(flow=(WOut, htOut, self.h2tc(htOut, FARo), PtOut, FARo),
                    Trq=self.TRQ * Pwrout * divby(Nmech),
                    NErr=100.0 if Ws1in == 0 else (Wcin - WcCalc) * divby(Wcin),
                    Nc=Nc, NcMap=NcMap, PRmap=PRmap, Eff=Eff, Tts1=Tts1)

    def nozzle(self, name: str, flow: tuple, Pamb: float, Ath: float) -> dict:
        """Nozzle_TMATS_body.c for a convergent nozzle with the Cfg thrust method, iDesign 2."""
        p = self.noz[name]
        WIn, _, TtIn, PtIn, FAR = flow
        gam, Rs, GJ = self.gam, self.Rt, self.G * self.J
        htin = self.t2hc(TtIn, FAR)
        Ptin = PtIn if PtIn > Pamb else Pamb + 0.1
        Ts_s, _hs, rhos_s, V_s = self._pcalc(Ptin, Pamb, TtIn, htin, FAR)
        MN_s = V_s * divby(sqrtT(gam * Rs * Ts_s * GJ))
        TsM = TtIn * divby(1 + 1 * 1 * (gam - 1) / 2)
        PsM = Ptin * powT(TsM * divby(TtIn), gam * divby(gam - 1))
        TsM, _h, rhoM, Vg = self._pcalc(Ptin, PsM, TtIn, htin, FAR)
        er = 1 - Vg * divby(sqrtT(gam * Rs * TsM * GJ))
        PsM_new = PsM + 0.05
        it = 0
        while abs(er) > 0.001 and it < 200:
            er_old, PsM_old = er, PsM
            PsM = PsM + 0.005 if abs(PsM - PsM_new) < 0.003 else PsM_new
            TsM, _h, rhoM, Vg = self._pcalc(Ptin, PsM, TtIn, htin, FAR)
            er = 1 - Vg * divby(sqrtT(gam * Rs * TsM * GJ))
            if abs(er) > 0.001:
                PsM_new = PsM - er * (PsM - PsM_old) * divby(er - er_old)
            it += 1
        choked = not (PsM < Pamb)
        if not choked:
            Psth, Vth, rhosth, MNth = Pamb, V_s, rhos_s, MN_s
        else:
            Psth, MNth, rhosth = PsM, 1.0, rhoM
            Vth = MNth * sqrtT(gam * Rs * TsM * GJ)
        PQPa = Ptin * divby(Pamb)
        CdTh = interp1Ac(p["PQ"], p["Cd"], PQPa)
        TG = interp1Ac(p["TGx"], p["TG"], TtIn)
        Wcalc = (1 - p["fl"] / 100) * (Ath * TG) * CdTh * rhosth * Vth / self.PSF
        Cfg = interp1Ac(p["PQ"], p["Cfg"], PQPa)
        return dict(Fg=(WIn / self.G) * V_s * Cfg,
                    NErr=100.0 if WIn == 0 else (WIn - Wcalc) * divby(WIn),
                    choked=choked, Wcalc=Wcalc, MNth=MNth, Psth=Psth)

    def valve(self, Pt_byp: float, pos: float, W: float, Tt: float, Pt: float) -> float:
        """Valve_TMATS_body.c: VBV bleed flow [lbm/s] from the LPC exit into the bypass."""
        frac = (pos - self.vlv_dead) * divby(self.vlv_open - self.vlv_dead)
        PR = Pt * divby(Pt_byp)
        if frac <= 0 or PR <= 1.0:
            return 0.0
        Wth = interp1Ac(self.vlv_PR, self.vlv_W, PR) * Pt * divby(sqrtT(Tt)) * (frac * self.vlv_Ae)
        return W if Wth > W else Wth

    # ------------------------------------------------------------------ the engine
    def residuals(self, z: Sequence[float], x: Sequence[float], ctrl: Controls, env: Env,
                  theta: dict | None = None) -> tuple[np.ndarray, dict]:
        """One pass through AGTF30_eng.mdl.

        z = [W, HPT_PR, HPC_Rline, FAN_Rline, BPR, LPC_Rline, LPT_PR], x = [N_LP, N_HP].
        r = [fan, LPC, HPC, HPT, LPT, core nozzle, bypass nozzle] normalized flow errors.
        """
        th = {k: 1.0 for k in THETA_KEYS}
        if theta:
            for k, v in theta.items():
                if k not in th:
                    raise KeyError("unknown health parameter " + k)
                th[k] = float(v)
        W, HPT_PR, HPC_R, FAN_R, BPR, LPC_R, LPT_PR = (float(v) for v in z)
        N_LP, N_HP = float(x[0]), float(x[1])
        ht0, Tt0, Pt0, FAR0, Pamb, Ts0, Veng = self.ambient(env)
        Fdrag = W * Veng / self.amb_g                                          # Ambient block
        Pt2 = Pt0 * self.eRam * lookup_clip(self.eRamVec, self.eRamTbl, Pt0 / Pamb)   # Inlet block
        s2 = (W, ht0, Tt0, Pt2, FAR0)
        N_fan = (1.0 / self.GR) * N_LP                                         # GearBox block
        fan = self.compressor(self.fan, s2, N_fan, FAN_R, th["fan_flow"], th["fan_eff"])
        Trq_fan_LP = fan["Trq"] * (1.0 / self.gear_eff) * (1.0 / self.GR)
        s21 = fan["flow"]
        bprm = max(BPR, 0.0)                                                   # Splitter block
        f = 1.0 / (bprm + 1.0)
        s13 = (bprm * f * s21[0],) + s21[1:]
        s22 = (f * s21[0],) + s21[1:]
        s23, MN2 = self.duct("2", s22)
        lpc = self.compressor(self.lpc, s23, N_LP, LPC_R, th["lpc_flow"], th["lpc_eff"])
        s24a = lpc["flow"]
        W_vbv = self.valve(s13[3], ctrl.VBV, s24a[0], s24a[2], s24a[3])        # Valve block
        byp = (s13[0] + W_vbv,) + s13[1:]              # bypass keeps its own ht, Tt, Pt, FAR
        s24 = (s24a[0] - W_vbv,) + s24a[1:]
        s25, MN25 = self.duct("25", s24)
        hpc = self.compressor(self.hpc, s25, N_HP, HPC_R, th["hpc_flow"], th["hpc_eff"])
        s36 = hpc["flow"]
        bl_lpt_exit, bl_hpt_exit, bl_hpt_inlet = hpc["bleeds"]
        s4 = self.burner(ctrl.Wf, s36)
        hpt = self.turbine(self.hpt, [bl_hpt_exit, bl_hpt_inlet], s4, HPT_PR, N_HP,
                           th["hpt_flow"], th["hpt_eff"])
        s45 = hpt["flow"]
        s48, MN45 = self.duct("45", s45)
        lpt = self.turbine(self.lpt, [bl_lpt_exit], s48, LPT_PR, N_LP, th["lpt_flow"], th["lpt_eff"])
        s5 = lpt["flow"]
        s7, MN5 = self.duct("5", s5)
        cnoz = self.nozzle("core", s7, Pamb, self.core_Ath)
        s17, MN17 = self.duct("17", byp)
        bnoz = self.nozzle("byp", s17, Pamb, ctrl.VAFN)
        r = np.array([fan["NErr"], lpc["NErr"], hpc["NErr"], hpt["NErr"], lpt["NErr"],
                      cnoz["NErr"], bnoz["NErr"]])
        # Shaft blocks: Ndot = (sum of torques - 5252.11*Pwr/N) * 60/(2 pi J)
        Trq_LPT_eff = lpt["Trq"] * self.LPS_eff
        Trq_LP = lpc["Trq"] + Trq_fan_LP + Trq_LPT_eff + self.shaft_gain * (self.pwr_LP / N_LP)
        Trq_HP = hpc["Trq"] + hpt["Trq"] + self.shaft_gain * (self.pwr_HP / N_HP)
        k = 60.0 / (2.0 * math.pi)
        Fg = bnoz["Fg"] + cnoz["Fg"]
        Fnet = Fg - Fdrag
        Ps3 = self.static_area(s36[0], s36[2], s36[3], s36[4], self.s36[0], self.s36[1])[1]
        out = dict(
            Trq_LP=Trq_LP, Trq_HP=Trq_HP, Ndot_LP=Trq_LP * k / self.I_LP, Ndot_HP=Trq_HP * k / self.I_HP,
            Fnet=Fnet, Fgross=Fg, Fdrag=Fdrag, Fg_byp=bnoz["Fg"], Fg_core=cnoz["Fg"], W_in=W,
            TSFC=ctrl.Wf * self.secs2hrs / Fnet if Fnet != 0 else float("inf"),
            BPR=BPR, OPR=s36[3] / Pt2, T4=s4[2], T45=s45[2], Ps3=Ps3,
            Pt2=Pt2, Tt2=Tt0, Pt21=s21[3], Tt21=s21[2], Pt25=s25[3], Tt25=s25[2], Pt3=s36[3], Tt3=s36[2],
            Pt4=s4[3], Pt45=s45[3], Pt5=s5[3], Tt5=s5[2], Pt7=s7[3], Pt17=s17[3], Pamb=Pamb, Tamb=Ts0,
            N_fan=N_fan, Nc_fan=fan["Nc"], N_LP=N_LP, N_HP=N_HP,
            SM_fan=fan["SM"], SM_lpc=lpc["SM"], SM_hpc=hpc["SM"],
            W25=s25[0], W_byp=byp[0], W13=s13[0], W_vbv=W_vbv, W3=s36[0], W4=s4[0], W5=s5[0],
            FAR4=s4[4], Wf=ctrl.Wf, VAFN=ctrl.VAFN, VBV=ctrl.VBV,
            Trq_fan=fan["Trq"], Trq_fan_LP=Trq_fan_LP, Trq_LPC=lpc["Trq"], Trq_LPT=lpt["Trq"],
            Trq_HPC=hpc["Trq"], Trq_HPT=hpt["Trq"], Tt41=hpt["Tts1"],
            MN_duct2=MN2, MN_duct25=MN25, MN_duct45=MN45, MN_duct5=MN5, MN_duct17=MN17,
            choked_core=float(cnoz["choked"]), choked_byp=float(bnoz["choked"]),
            PR_fan=fan["PR"], PR_lpc=lpc["PR"], PR_hpc=hpc["PR"],
        )
        return r, out

    # ------------------------------------------------------------------ inner solver
    def _fd_jac(self, z: np.ndarray, r0: np.ndarray, x, ctrl, env, theta) -> np.ndarray:
        """Forward-difference Jacobian dr/dz (7x7), relative step 1e-6."""
        Jm = np.empty((7, 7))
        for j in range(7):
            h = 1e-6 * max(abs(z[j]), 1e-3)
            zp = z.copy()
            zp[j] += h
            Jm[:, j] = (self.residuals(zp, x, ctrl, env, theta)[0] - r0) / h
        return Jm

    def _newton(self, z0, x, ctrl, env, theta, tol, max_iter, jac) -> tuple[np.ndarray, dict, dict]:
        z = np.array(z0, dtype=float)
        info = dict(converged=False, iters=0, resnorm=float("inf"), nfev=0)
        try:
            r, out = self.residuals(z, x, ctrl, env, theta)
        except _NUMERIC_ERRORS:
            return z, {}, info
        info["nfev"] = 1
        res = float(np.max(np.abs(r)))
        Jm = None if jac is None else np.array(jac, dtype=float)
        fresh = False
        for it in range(max_iter + 1):
            info.update(iters=it, resnorm=res)
            if res < tol:
                info["converged"] = True
                break
            if it == max_iter:
                break
            if Jm is None:
                try:
                    Jm = self._fd_jac(z, r, x, ctrl, env, theta)
                except _NUMERIC_ERRORS:
                    break
                info["nfev"] += 7
                fresh = True
            try:
                dz = np.linalg.solve(Jm, -r)
            except np.linalg.LinAlgError:
                if fresh:
                    break
                Jm = None
                continue
            lam, ok = 1.0, False
            while lam >= 1.0 / 1024:
                zn = z + lam * dz
                try:
                    rn, on = self.residuals(zn, x, ctrl, env, theta)
                    info["nfev"] += 1
                    resn = float(np.max(np.abs(rn)))
                    if math.isfinite(resn) and resn < res:
                        ok = True
                        break
                except _NUMERIC_ERRORS:
                    pass
                lam *= 0.5
            if not ok:
                if fresh:
                    break                       # even a fresh Jacobian gives no descent: give up
                Jm = None
                continue
            stalled = resn > 0.3 * res
            z, r, out, res = zn, rn, on, resn
            fresh = False
            if stalled:
                Jm = None                       # chord iteration stalls -> refresh the Jacobian
        return z, out, info

    def _anchor(self, env: Env) -> dict:
        return NASA_IC["cruise"] if env.alt_ft >= 20000.0 else NASA_IC["takeoff"]

    def solve_flow(self, x: Sequence[float], ctrl: Controls, env: Env, theta: dict | None = None,
                   tol: float = 1e-12, z0: Sequence[float] | None = None, max_iter: int = 60,
                   jac: np.ndarray | None = None) -> tuple[np.ndarray, dict, dict]:
        """Newton(-chord) solve of r(z) = 0 for fixed x. Only stopping rule: max|r| < tol.

        Without z0 the start is the nearest NASA initial condition; if that fails, a
        deterministic continuation walks from that NASA point to the requested one.
        """
        a = self._anchor(env)
        start = np.array(a["z"] if z0 is None else z0, dtype=float)
        z, out, info = self._newton(start, x, ctrl, env, theta, tol, max_iter, jac)
        if info["converged"]:
            return z, out, info
        nfev = info["nfev"]
        p0 = np.array([a["N2"], a["N3"], a["Wf"], a["VAFN"], a["VBV"], a["alt"], a["MN"], a["dT"]])
        p1 = np.array([x[0], x[1], ctrl.Wf, ctrl.VAFN, ctrl.VBV, env.alt_ft, env.MN, env.dT])

        def at(s: float):
            p = p0 + s * (p1 - p0)
            return p[:2], Controls(p[2], p[3], p[4]), Env(p[5], p[6], p[7])

        zc, _, inf0 = self._newton(np.array(a["z"]), *at(0.0), theta, 1e-9, max_iter, None)
        nfev += inf0["nfev"]
        s, ds = 0.0, 0.25
        while inf0["converged"] and s < 1.0 and ds > 1e-4:
            sn = min(1.0, s + ds)
            last = sn >= 1.0
            if last:
                zt, ot, it_ = self._newton(zc, x, ctrl, env, theta, tol, max_iter, None)
            else:
                zt, ot, it_ = self._newton(zc, *at(sn), theta, 1e-9, max_iter, None)
            nfev += it_["nfev"]
            if it_["converged"]:
                s, zc, ds = sn, zt, min(ds * 2, 0.5)
                if last:
                    it_["nfev"] = nfev
                    return zt, ot, it_
            else:
                ds *= 0.5
        info["nfev"] = nfev
        return z, out, info

    def F(self, x: Sequence[float], ctrl: Controls, env: Env, theta: dict | None = None,
          tol: float = 1e-12, z0: Sequence[float] | None = None,
          jac: np.ndarray | None = None) -> np.ndarray:
        """State derivative [dN_LP/dt, dN_HP/dt] in rpm/s after the inner flow solve."""
        z, out, info = self.solve_flow(x, ctrl, env, theta, tol=tol, z0=z0, jac=jac)
        self.last_info = dict(info, z=z)
        if not info["converged"]:
            warnings.warn(f"replica solve_flow not converged (max|r| = {info['resnorm']:.3e})", RuntimeWarning)
            if not out:
                return np.array([float("nan"), float("nan")])
        return np.array([out["Ndot_LP"], out["Ndot_HP"]])


@lru_cache(maxsize=4)
def _load(vendor: str) -> Engine:
    return Engine(Sources(vendor))


def load_engine(vendor_dir: str | None = None) -> Engine:
    """Parse the NASA sources (maps, .mdl masks/wiring, C gas tables) and build the engine.

    Default vendor folder: the 'vendor' directory next to the 'replica' package directory.
    """
    v = Path(vendor_dir) if vendor_dir is not None else Path(__file__).resolve().parent.parent / "vendor"
    return _load(str(v.resolve()))
