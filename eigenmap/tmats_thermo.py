"""Gas properties of T-MATS, ported from the C sources
`vendor/tmats/TMATS_Library/MEX/{t2hc,h2tc,pt2sc,sp2tc}_TMATS.c` (D. K. Frederick's routines).

The gas is a mixture of air and stoichiometric combustion products; its enthalpy
h(T) and entropy function phi(T) are piecewise cubic polynomials in temperature,
blended with the fuel-air ratio `fa`.

Two families of inverse functions are provided:

* `h2t`, `sp2t`      : FAITHFUL copies of the T-MATS secant loops, including their
                       loose stopping rules (|dh| < 1e-3 BTU/lbm, |dS| < 1e-4) and
                       iteration caps. These reproduce NASA's numbers bit-for-bit
                       in spirit, but they leave a small non-smooth error.
* `h2t_tight`, `sp2t_tight` : the same equations solved by Newton's method to
                       machine precision (DERIVED deviation, used when the study
                       needs F(x,u) to be a smooth function so that the ONLY
                       solver noise is the tunable tolerance `tol`).

Units: T in deg R, P in psia, h in BTU/lbm, S in BTU/(lbm R).
Course topic: Engine model (not a course topic)
"""
from __future__ import annotations

import math

from .agtf30_data import load_params

_G = load_params()["gas"]
_ITAB = [int(v) - 1 for v in _G["ITAB"]]           # 60 entries -> index of the polynomial piece
_TT = [float(v) for v in _G["TTTAB"]]              # 11 piece start temperatures
_AH, _BH, _CH, _DH = _G["AHAIR"], _G["BHAIR"], _G["CHAIR"], _G["DHAIR"]
_AHS, _BHS, _CHS, _DHS = _G["AHSTOC"], _G["BHSTOC"], _G["CHSTOC"], _G["DHSTOC"]
_AP, _BP, _CP, _DP = _G["APAIR"], _G["BPAIR"], _G["CPAIR"], _G["DPAIR"]
_APS, _BPS, _CPS, _DPS = _G["APSTOC"], _G["BPSTOC"], _G["CPSTOC"], _G["DPSTOC"]

# piece index for every 100 R band, exactly as the C code computes it
# (Tindex = clamp(int(T/100), 1, 60); it = ITAB[Tindex-1]-1)
_PIECE = [_ITAB[min(max(k, 1), 60) - 1] for k in range(0, 200)]


def _piece(T: float) -> int:
    """Index of the polynomial piece used for temperature T (same rule as the C code)."""
    k = int(0.01 * (T - math.fmod(T, 100.0)))
    if k < 1:
        return _PIECE[1]
    if k > 60:
        return _PIECE[60]
    return _PIECE[k]


def _mix(fa: float) -> tuple[float, float, float, float]:
    """Mixture mole numbers (zmea, zmsp, tmlsr, zmwtr) for fuel-air ratio fa (eqns 35-42 of the C source)."""
    if fa == 0.0:
        return 4.7642, 0.0, 0.2098988288, 0.0345194683
    zmea = 4.7642 - fa * 69.69056873
    zmsp = fa * 74.411931335
    tmls = 4.7642 + fa * 4.721362582
    zmwtr = tmls / (138.0148721 * (1.0 + fa))
    return zmea, zmsp, 1.0 / tmls, zmwtr


def t2h(T: float, fa: float) -> float:
    """Enthalpy h(T, fa) [BTU/lbm]. Port of t2hc()."""
    zmea, zmsp, tmlsr, zmwtr = _mix(fa)
    it = _piece(T)
    dl = T - _TT[it]
    hgea = ((_DH[it] * dl + _CH[it]) * dl + _BH[it]) * dl + _AH[it]
    hgsp = (((_DHS[it] * dl + _CHS[it]) * dl + _BHS[it]) * dl + _AHS[it]) if fa > 0.0 else 0.0
    return (hgsp * zmsp + hgea * zmea) * (tmlsr * zmwtr)


def _dh_dT(T: float, fa: float) -> float:
    """dh/dT (specific heat cp) from the same polynomials (used by the tight Newton solver)."""
    zmea, zmsp, tmlsr, zmwtr = _mix(fa)
    it = _piece(T)
    dl = T - _TT[it]
    d_ea = (3.0 * _DH[it] * dl + 2.0 * _CH[it]) * dl + _BH[it]
    d_sp = ((3.0 * _DHS[it] * dl + 2.0 * _CHS[it]) * dl + _BHS[it]) if fa > 0.0 else 0.0
    return (d_sp * zmsp + d_ea * zmea) * (tmlsr * zmwtr)


def h2t(H: float, fa: float) -> float:
    """Temperature from enthalpy: FAITHFUL port of h2tc() (secant, stops at |dh| <= 1e-3 or 11 iterations)."""
    zmea, zmsp, tmlsr, zmwtr = _mix(fa)
    tg = 3.55 * H
    tgo = 1.0
    hgo = 0.282
    zz = tmlsr * zmwtr
    ii = 0
    hh = 99.99
    while abs(hh) > 1e-3 and ii < 11:
        ii += 1
        it = _piece(tg)
        d1 = tg - _TT[it]
        hgea = ((_DH[it] * d1 + _CH[it]) * d1 + _BH[it]) * d1 + _AH[it]
        hgsp = (((_DHS[it] * d1 + _CHS[it]) * d1 + _BHS[it]) * d1 + _AHS[it]) if fa > 0.0 else 0.0
        hg = (hgsp * zmsp + hgea * zmea) * zz
        hh = H - hg
        d = (tg - tgo) / (hg - hgo)
        tgo = tg
        hgo = hg
        tg = tg + d * hh
    return tg


def h2t_tight(H: float, fa: float) -> float:
    """Temperature from enthalpy solved to machine precision (Newton's method on h(T) - H = 0)."""
    tg = 3.55 * H
    for _ in range(40):
        f = t2h(tg, fa) - H
        step = f / _dh_dT(tg, fa)
        tg -= step
        if abs(step) <= 2e-14 * abs(tg):
            break
    return tg


def _phi(T: float, fa: float, zmea: float, zmsp: float, tmlsr: float) -> float:
    it = _piece(T)
    dl = T - _TT[it]
    phiea = ((_DP[it] * dl + _CP[it]) * dl + _BP[it]) * dl + _AP[it]
    phisp = (((_DPS[it] * dl + _CPS[it]) * dl + _BPS[it]) * dl + _APS[it]) if fa > 0.0 else 0.0
    return (phisp * zmsp + phiea * zmea) * tmlsr


def _dphi_dT(T: float, fa: float, zmea: float, zmsp: float, tmlsr: float) -> float:
    it = _piece(T)
    dl = T - _TT[it]
    d_ea = (3.0 * _DP[it] * dl + 2.0 * _CP[it]) * dl + _BP[it]
    d_sp = ((3.0 * _DPS[it] * dl + 2.0 * _CPS[it]) * dl + _BPS[it]) if fa > 0.0 else 0.0
    return (d_sp * zmsp + d_ea * zmea) * tmlsr


def pt2s(P: float, T: float, fa: float) -> float:
    """Entropy S(P, T, fa). Port of pt2sc()."""
    zmea, zmsp, tmlsr, zmwtr = _mix(fa)
    phig = _phi(T, fa, zmea, zmsp, tmlsr)
    rcas = 1.98587 * zmwtr
    return (phig * 0.5035576347 - 23.0258509) * rcas - 0.1841304 - rcas * math.log(P / 14.696)


def sp2t(S: float, P: float, fa: float) -> float:
    """Temperature from entropy and pressure: FAITHFUL port of sp2tc().

    Secant iteration from 1000 R, stops when |S - Sg| < 1e-4 (but at least 2 and at most 10 iterations).
    """
    zmea, zmsp, tmlsr, zmwtr = _mix(fa)
    rcas = 1.98587 * zmwtr
    lnp = rcas * math.log(P / 14.696)
    Tg = 1000.0
    jj = 0
    Sg = 1e3
    Sg1 = Tg1 = 0.0
    while (abs(S - Sg) >= 1e-4 and jj < 10) or jj < 2:
        phig = _phi(Tg, fa, zmea, zmsp, tmlsr)
        Sg = (phig * 0.5035576347 - 23.0258509) * rcas - 0.1841304 - lnp
        if jj == 0:
            Sg1 = Sg
            Tg1 = Tg
            Tg = Tg1 + 50.0
        else:
            dTdS = (Tg - Tg1) / (Sg - Sg1)
            Tg1 = Tg
            Sg1 = Sg
            Tg = Tg1 + (S - Sg1) * dTdS
        jj += 1
    return Tg


def sp2t_tight(S: float, P: float, fa: float, T_guess: float = 1000.0) -> float:
    """Temperature from entropy and pressure solved to machine precision (Newton's method)."""
    zmea, zmsp, tmlsr, zmwtr = _mix(fa)
    rcas = 1.98587 * zmwtr
    const = -23.0258509 * rcas - 0.1841304 - rcas * math.log(P / 14.696)
    k = 0.5035576347 * rcas
    Tg = T_guess
    for _ in range(40):
        f = _phi(Tg, fa, zmea, zmsp, tmlsr) * k + const - S
        step = f / (_dphi_dT(Tg, fa, zmea, zmsp, tmlsr) * k)
        Tg -= step
        if Tg < 150.0:
            Tg = 150.0
        if abs(step) <= 2e-14 * abs(Tg):
            break
    return Tg


if __name__ == "__main__":
    # quick self-check: faithful and tight inverses agree to the T-MATS stopping tolerance
    import random
    rnd = random.Random(1)
    worst_h = worst_s = worst_rt = 0.0
    for _ in range(2000):
        T = rnd.uniform(380.0, 3400.0)
        fa = rnd.choice([0.0, rnd.uniform(0.0, 0.03)])
        P = rnd.uniform(3.0, 800.0)
        h = t2h(T, fa)
        s = pt2s(P, T, fa)
        worst_rt = max(worst_rt, abs(h2t_tight(h, fa) - T) / T, abs(sp2t_tight(s, P, fa) - T) / T)
        worst_h = max(worst_h, abs(h2t(h, fa) - T) / T)
        worst_s = max(worst_s, abs(sp2t(s, P, fa) - T) / T)
    print("tight round-trip max rel error      : %.2e" % worst_rt)
    print("faithful h2t   max rel error        : %.2e" % worst_h)
    print("faithful sp2t  max rel error        : %.2e" % worst_s)
