"""Linearization of the engine about a steady operating point.

    dx/dt = F(x, u; theta),   x = [N_LP, N_HP]^T (rpm),  u = W_f (lbm/s)
    F(x0, u0) = 0             (steady state)
    F(x0 + dx, u0 + du) ~ F(x0,u0) + (dF/dx) dx + (dF/du) du     (first-order Taylor series)
    =>  d(dx)/dt = A dx + B du,   A = dF/dx (2x2 Jacobian),  B = dF/du (2x1)

The partial derivatives are taken numerically (finite differences), because F is
only available as a computer program: every evaluation of F solves the engine's
7 flow-matching equations with an inner Newton solver whose tolerance is `tol`.

Held constant during the perturbation (documented in ASSUMPTIONS.md): the
variable-area fan nozzle (VAFN) and the variable bleed valve (VBV) stay frozen at
their steady scheduled values; ambient conditions and health parameters theta.

Owned by the orchestrator (written at Gate 1 on the frozen engine).
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from . import mathcore as mc
from .engine_model import Controls

# Base relative step of the reference Jacobian. The engine maps are piecewise
# (bi)linear tables (as in T-MATS), so F has small kinks at table breakpoints;
# the "study" engine uses C1 maps, whose second derivative still jumps at breakpoints, so a
# central difference that straddles a breakpoint has an O(h) error; 1e-5 keeps that below 1e-6
# while solver noise at tol = 1e-13 stays near 1e-8 relative.
H_REF = 1.0e-5


def inner_jacobian(engine, sp, rel_step: float = 1.0e-6) -> np.ndarray:
    """7x7 Jacobian dr/dz of the flow-matching residuals at the steady point.

    Computed once per operating point by central differences and re-used by the
    inner solver (chord iterations), exactly like T-MATS re-uses its Jacobian.
    Course topic: Jacobians / partial derivatives
    """
    z0 = np.asarray(sp.z0, dtype=float)
    J = np.zeros((7, 7))
    for j in range(7):
        dz = rel_step * max(abs(z0[j]), 1e-3)
        zp, zm = z0.copy(), z0.copy()
        zp[j] += dz
        zm[j] -= dz
        rp, _ = engine.residuals(zp, sp.x0, sp.ctrl, sp.env, sp.theta)
        rm, _ = engine.residuals(zm, sp.x0, sp.ctrl, sp.env, sp.theta)
        J[:, j] = (np.asarray(rp) - np.asarray(rm)) / (2.0 * dz)
    return J


def start_vector(sp, start: str = "converged") -> np.ndarray:
    """Starting guess for the inner solver.

    'converged': the operating point's own converged z0 (what a linearization program normally does);
    'offset'   : z0 moved by +-0.1 % (alternating signs) - a start that is NOT the solution, so the
                 solver always iterates and stops somewhere below tol (the classical 'noisy F' case).
    """
    z0 = np.asarray(sp.z0, dtype=float)
    if start == "converged":
        return z0
    if start == "offset":
        return z0 * (1.0 + 1e-3 * np.array([1, -1, 1, -1, 1, -1, 1], dtype=float))
    raise ValueError(start)


def make_F(engine, sp, tol: float, jac: np.ndarray | None = None,
           start: str = "converged") -> Callable[[np.ndarray], np.ndarray]:
    """Return f(v) = F(x, u) with v = [N_LP, N_HP, W_f]; VAFN and VBV frozen at sp.ctrl.

    Each call starts the inner flow solver from `start_vector(sp, start)` and stops it at
    max|r| < tol - so `tol` is exactly the largest solver residual that can be left in F.
    Course topic: Engine model (not a course topic)
    """
    z_start = start_vector(sp, start)

    def f(v: np.ndarray) -> np.ndarray:
        ctrl = Controls(Wf=float(v[2]), VAFN=sp.ctrl.VAFN, VBV=sp.ctrl.VBV)
        return np.asarray(engine.F([float(v[0]), float(v[1])], ctrl, sp.env, sp.theta,
                                   tol=tol, z0=z_start, jac=jac), dtype=float)
    return f


def v0_of(sp) -> np.ndarray:
    """Operating point as the 3-vector [N_LP0, N_HP0, W_f0]."""
    return np.array([sp.x0[0], sp.x0[1], sp.ctrl.Wf], dtype=float)


def split_AB(J: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Split the 2x3 Jacobian d F / d[N_LP, N_HP, W_f] into A (2x2) and B (2x1).

    Course topic: Jacobians / partial derivatives
    """
    return J[:, :2].copy(), J[:, 2:3].copy()


def normalized_A(A: np.ndarray, x0) -> np.ndarray:
    """A_n = D^-1 A D with D = diag(x0): the state matrix for fractional speeds dx_i/x0_i.

    A similarity transform: same eigenvalues, different eigenvectors and kappa(V).
    Course topic: Eigenvalues (similarity transformations)
    """
    D = np.diag(np.asarray(x0, dtype=float))
    return np.linalg.solve(D, A @ D)


def fd_AB(engine, sp, h: float, scheme: str, tol: float,
          jac: np.ndarray | None = None, start: str = "converged") -> tuple[np.ndarray, np.ndarray]:
    """A and B by plain forward or central differences with relative step h.

    forward:  dF/dv_j ~ [F(v + h v_j e_j) - F(v)] / (h v_j)          error O(h)
    central:  dF/dv_j ~ [F(v + h v_j e_j) - F(v - h v_j e_j)] / (2 h v_j)   error O(h^2)
    Course topic: Jacobians / partial derivatives; Taylor series (truncation error)
    """
    f = make_F(engine, sp, tol, jac, start)
    J = mc.fd_jacobian(f, v0_of(sp), h, scheme=scheme, rel=True)
    return split_AB(np.asarray(J, dtype=float))


def reference_AB(engine, sp, tol: float = 1.0e-13, jac: np.ndarray | None = None,
                 h: float = H_REF) -> dict:
    """Reference linear model: central differences + Richardson extrapolation at tight tol.

    Self-consistency: the same calculation is repeated with half the base step;
    `rel_diff` is the largest relative difference between the two A matrices
    (must be < 1e-6 for "6 consistent digits").
    Course topic: Jacobians / partial derivatives; Taylor series
    """
    if jac is None:
        jac = inner_jacobian(engine, sp)
    f = make_F(engine, sp, tol, jac)
    v0 = v0_of(sp)
    J1, err1 = mc.richardson_jacobian(f, v0, h, levels=2)
    J2, _ = mc.richardson_jacobian(f, v0, h / 2.0, levels=2)
    J1 = np.asarray(J1, dtype=float)
    J2 = np.asarray(J2, dtype=float)
    A, B = split_AB(J1)
    A2, B2 = split_AB(J2)
    scale = np.max(np.abs(A))
    rel_diff = float(np.max(np.abs(A - A2)) / scale)
    relB = float(np.max(np.abs(B - B2)) / np.max(np.abs(B)))
    return {
        "A": A, "B": B, "A_n": normalized_A(A, sp.x0),
        "rel_diff_A": rel_diff, "rel_diff_B": relB,
        "err_est": float(np.max(np.abs(np.asarray(err1)))),
        "h_ref": h, "tol": tol,
        "method": "central difference + Richardson (2 levels), relative step %g and %g" % (h, h / 2),
        "F0": f(v0),
    }
