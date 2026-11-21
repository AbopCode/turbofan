"""Math core of the eigenvalue study: OWN 2x2 linear algebra and finite differences.

Everything here is engine-independent pure mathematics. The routines are written
out from the course formulas on purpose (no numpy.linalg / scipy inside them);
numpy is used only to store arrays and for elementary functions (+ - * /, sqrt,
exp, sin, cos). Every routine is cross-checked against numpy / scipy in
`tests/test_mathcore.py`.

Notation follows CONTRACT.md section 2: A (2x2 state matrix), lam (eigenvalues),
tau (time constants), V (eigenvector matrix, unit 2-norm columns), kappaV
(2-norm condition number of V), E (perturbation matrix), h (finite-difference step).

Course topic: Eigenvalues / Diagonalization / Linear ODE systems / Gauss-Seidel /
Jacobians and partial derivatives / Taylor series / Matrix norms and conditioning
"""
from __future__ import annotations

from typing import Callable, Sequence

import numpy as np

# ---------------------------------------------------------------------------
# Numerical switches (documented, used by diagonalize / expm_diag)
# ---------------------------------------------------------------------------
#: |det V| (V has unit columns, so |det V| = |sin(angle between eigenvectors)|)
#: below which A is treated as numerically DEFECTIVE (no usable diagonalization).
DEFECTIVE_DET_TOL = 1e-10
#: expm_diag leaves the diagonalization formula and uses the closed-form
#: (Cayley-Hamilton) formula when kappa2(V) exceeds this value. The rounding
#: error of P e^{Dt} P^-1 grows like machine-eps * kappa2(V).
EXPM_KAPPA_SWITCH = 1e3


class DefectiveMatrixError(ValueError):
    """Raised by `diagonalize` when A has no (numerically usable) eigenvector basis.

    Course topic: Diagonalization
    """


# ---------------------------------------------------------------------------
# Small private helpers
# ---------------------------------------------------------------------------
def _real2x2(A) -> np.ndarray:
    """Return A as a real float 2x2 array (error if it is not real 2x2).

    Course topic: (none - input checking)
    """
    M = np.asarray(A)
    if M.shape != (2, 2):
        raise ValueError(f"expected a 2x2 matrix, got shape {M.shape}")
    if np.iscomplexobj(M):
        if np.any(M.imag != 0.0):
            raise ValueError("expected a REAL 2x2 matrix")
        M = M.real
    return M.astype(float)


def _abs2(z):
    """|z|^2 = Re(z)^2 + Im(z)^2 (works for real or complex numbers / arrays).

    Course topic: (none - complex arithmetic)
    """
    z = np.asarray(z)
    return z.real * z.real + z.imag * z.imag


def _input_vector(B, du) -> np.ndarray:
    """Return the constant forcing vector g = B du as a length-2 array.

    B may be (2,), (2,1) or (2,m); du a scalar or a length-m vector.
    g_i = sum_k B[i,k] * du[k]   (written out, no library matrix product).

    Course topic: Linear ODE systems
    """
    Bm = np.asarray(B, dtype=float)
    if Bm.ndim == 1:
        Bm = Bm.reshape(2, 1)
    duv = np.atleast_1d(np.asarray(du, dtype=float)).ravel()
    if Bm.shape[0] != 2 or Bm.shape[1] != duv.size:
        raise ValueError(f"B shape {Bm.shape} does not match du of length {duv.size}")
    g = np.zeros(2)
    for k in range(duv.size):
        g[0] += Bm[0, k] * duv[k]
        g[1] += Bm[1, k] * duv[k]
    return g


# ---------------------------------------------------------------------------
# Eigenvalues
# ---------------------------------------------------------------------------
def eig2x2(A) -> tuple[np.ndarray, dict]:
    """Eigenvalues of a real 2x2 matrix from its characteristic equation.

    Formula: det(A - lam I) = lam^2 - tr(A) lam + det(A) = 0, so
        lam = ( tr +/- sqrt(disc) ) / 2,   disc = tr^2 - 4 det.
    Two numerical safeguards (same mathematics, fewer rounding errors):
      * disc is evaluated as (a11 - a22)^2 + 4 a12 a21, which is algebraically
        equal to tr^2 - 4 det but does not subtract two nearly equal numbers;
      * for real roots the LARGER-magnitude root is q = (tr + sign(tr) sqrt(disc))/2
        (no cancellation) and the other one is det / q (Vieta: lam1 lam2 = det).
    disc < 0 gives the complex pair tr/2 +/- i sqrt(-disc)/2; disc == 0 the
    repeated root tr/2.

    Returns (lam, info): lam is a complex array of length 2, sorted "slow first":
    a real pair by |Re| ascending (ties: larger value first), a complex pair with
    the +imaginary member first. info = {'tr', 'det', 'disc'}.

    Course topic: Eigenvalues
    """
    M = _real2x2(A)
    a, b, c, d = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    tr = a + d
    det = a * d - b * c
    disc = (a - d) * (a - d) + 4.0 * b * c
    info = {"tr": float(tr), "det": float(det), "disc": float(disc)}
    if disc < 0.0:                                   # complex conjugate pair
        w = 0.5 * np.sqrt(-disc)
        lam = np.array([complex(0.5 * tr, w), complex(0.5 * tr, -w)])
        return lam, info
    if disc == 0.0:                                  # repeated real root
        return np.array([complex(0.5 * tr), complex(0.5 * tr)]), info
    s = np.sqrt(disc)
    q = 0.5 * (tr + s) if tr >= 0.0 else 0.5 * (tr - s)   # larger-magnitude root
    l_big = q
    l_small = det / q if q != 0.0 else 0.0
    pair = sorted([float(l_big), float(l_small)], key=lambda v: (abs(v), -v))
    return np.array([complex(pair[0]), complex(pair[1])]), info


def eig2x2_batch(As) -> np.ndarray:
    """Eigenvalues of MANY real 2x2 matrices at once (array shape (n, 2, 2)).

    Same formula and same ordering as `eig2x2`, evaluated with numpy arrays
    instead of a Python loop: lam = (tr +/- sqrt(disc))/2 with
    disc = (a11 - a22)^2 + 4 a12 a21 and the stable pair q, det/q.
    Returns a complex array of shape (n, 2).

    Course topic: Eigenvalues
    """
    M = np.asarray(As, dtype=float)
    if M.ndim != 3 or M.shape[1:] != (2, 2):
        raise ValueError(f"expected shape (n, 2, 2), got {M.shape}")
    a, b, c, d = M[:, 0, 0], M[:, 0, 1], M[:, 1, 0], M[:, 1, 1]
    tr = a + d
    det = a * d - b * c
    disc = (a - d) * (a - d) + 4.0 * b * c
    s = np.sqrt(np.abs(disc))
    # --- real branch (disc >= 0) ---
    q = np.where(tr >= 0.0, 0.5 * (tr + s), 0.5 * (tr - s))
    safe_q = np.where(q != 0.0, q, 1.0)
    other = np.where(q != 0.0, det / safe_q, 0.0)
    q = np.where(disc == 0.0, 0.5 * tr, q)
    other = np.where(disc == 0.0, 0.5 * tr, other)
    # slow first: smaller |value| first; on a tie the larger value first
    q_first = (np.abs(q) < np.abs(other)) | ((np.abs(q) == np.abs(other)) & (q >= other))
    first = np.where(q_first, q, other)
    second = np.where(q_first, other, q)
    # --- complex branch (disc < 0) ---
    is_cplx = disc < 0.0
    lam = np.empty((M.shape[0], 2), dtype=complex)
    lam[:, 0] = np.where(is_cplx, 0.5 * tr + 0.5j * s, first)
    lam[:, 1] = np.where(is_cplx, 0.5 * tr - 0.5j * s, second)
    return lam


def eigvec2x2(A, lam) -> np.ndarray:
    """Eigenvector matrix V of a real 2x2 matrix (columns = eigenvectors, unit 2-norm).

    Formula: (A - lam I) v = 0. Row 1 gives v = [a12, lam - a11]^T, row 2 gives
    v = [lam - a22, a21]^T. Both are valid; the one with the larger norm is used
    (the other may be the zero vector, e.g. for a triangular A). Each column is
    divided by its 2-norm and scaled so that its largest component is real and
    positive (an eigenvector is only defined up to a factor).
    Special cases: A = c I (every vector is an eigenvector) returns V = I;
    a defective A (repeated lam, only ONE eigenvector) returns two equal columns,
    so V is singular - see `diagonalize`.

    Course topic: Eigenvalues
    """
    M = _real2x2(A)
    lam = np.asarray(lam, dtype=complex).ravel()
    if lam.size != 2:
        raise ValueError("lam must hold the two eigenvalues")
    a, b, c, d = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    scale = max(abs(a), abs(b), abs(c), abs(d), 1e-300)
    V = np.zeros((2, 2), dtype=complex)
    for i in range(2):
        v_row1 = np.array([b, lam[i] - a], dtype=complex)
        v_row2 = np.array([lam[i] - d, c], dtype=complex)
        n1 = np.sqrt(_abs2(v_row1[0]) + _abs2(v_row1[1]))
        n2 = np.sqrt(_abs2(v_row2[0]) + _abs2(v_row2[1]))
        v, nv = (v_row1, n1) if n1 >= n2 else (v_row2, n2)
        if nv <= 1e-14 * scale:                      # A - lam I = 0: A is a multiple of I
            v = np.array([1.0, 0.0] if i == 0 else [0.0, 1.0], dtype=complex)
            nv = 1.0
        v = v / nv
        k = 0 if abs(v[0]) >= abs(v[1]) else 1       # largest component -> real, positive
        v = v * (np.conj(v[k]) / abs(v[k]))
        V[:, i] = v / np.sqrt(_abs2(v[0]) + _abs2(v[1]))
    if np.all(lam.imag == 0.0):
        return V.real.copy()
    return V


# ---------------------------------------------------------------------------
# 2x2 inverse, singular values, norms, condition number
# ---------------------------------------------------------------------------
def inv2x2(M) -> np.ndarray:
    """Inverse of a 2x2 matrix (real or complex) by the adjugate formula.

    Formula: [[a, b], [c, d]]^-1 = (1 / (a d - b c)) [[d, -b], [-c, a]].
    Raises ValueError when the determinant is exactly zero.

    Course topic: Diagonalization
    """
    M = np.asarray(M)
    if M.shape != (2, 2):
        raise ValueError(f"expected a 2x2 matrix, got shape {M.shape}")
    a, b, c, d = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    det = a * d - b * c
    if det == 0:
        raise ValueError("matrix is singular (determinant = 0)")
    return np.array([[d, -b], [-c, a]]) / det


def singular_values_2x2(M) -> tuple[float, float]:
    """The two singular values (sigma_max, sigma_min) of a 2x2 real or complex matrix.

    Formula: the singular values are the square roots of the eigenvalues of the
    Hermitian matrix G = M^H M = [[g11, g12], [conj(g12), g22]]. Its characteristic
    equation gives
        sigma_max^2 = ( g11 + g22 + sqrt( (g11 - g22)^2 + 4 |g12|^2 ) ) / 2,
        sigma_min   = |det M| / sigma_max        (because sigma_max sigma_min = |det M|).

    Course topic: Matrix norms and conditioning
    """
    M = np.asarray(M, dtype=complex)
    if M.shape != (2, 2):
        raise ValueError(f"expected a 2x2 matrix, got shape {M.shape}")
    a, b, c, d = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    g11 = _abs2(a) + _abs2(c)
    g22 = _abs2(b) + _abs2(d)
    g12 = np.conj(a) * b + np.conj(c) * d
    root = np.sqrt((g11 - g22) ** 2 + 4.0 * _abs2(g12))
    smax = float(np.sqrt(0.5 * (g11 + g22 + root)))
    if smax == 0.0:
        return 0.0, 0.0
    smin = float(abs(a * d - b * c)) / smax
    return smax, smin


def norm2(M) -> float:
    """Spectral norm ||M||_2 of a 2x2 matrix = its largest singular value.

    Formula: ||M||_2 = sigma_max(M), see `singular_values_2x2`.

    Course topic: Matrix norms and conditioning
    """
    return singular_values_2x2(M)[0]


def norm2_batch(Ms) -> np.ndarray:
    """Spectral norms of many REAL 2x2 matrices (array shape (n, 2, 2)).

    Formula (same as `singular_values_2x2`, vectorised):
        sigma_max^2 = ( g11 + g22 + sqrt((g11 - g22)^2 + 4 g12^2) ) / 2,  G = M^T M.

    Course topic: Matrix norms and conditioning
    """
    M = np.asarray(Ms, dtype=float)
    a, b, c, d = M[:, 0, 0], M[:, 0, 1], M[:, 1, 0], M[:, 1, 1]
    g11 = a * a + c * c
    g22 = b * b + d * d
    g12 = a * b + c * d
    return np.sqrt(0.5 * (g11 + g22 + np.sqrt((g11 - g22) ** 2 + 4.0 * g12 * g12)))


def kappa2(V) -> float:
    """2-norm condition number of a 2x2 (real or complex) matrix V.

    Formula: kappa_2(V) = sigma_max(V) / sigma_min(V), with the closed-form 2x2
    singular values of `singular_values_2x2`. Returns inf for a singular V.
    For an eigenvector matrix with unit columns kappa_2 = 1 means orthogonal
    eigenvectors; a large value means nearly parallel eigenvectors.

    Course topic: Matrix norms and conditioning
    """
    smax, smin = singular_values_2x2(V)
    if smin == 0.0:
        return float("inf")
    return smax / smin


# ---------------------------------------------------------------------------
# Diagonalization, stability, time constants
# ---------------------------------------------------------------------------
def diagonalize(A) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Diagonalize a real 2x2 matrix: A = P D P^-1.

    Formula: P = V (eigenvector columns, unit 2-norm), D = diag(lam1, lam2),
    P^-1 from the 2x2 adjugate formula. All three are returned as COMPLEX arrays
    so that a complex-conjugate eigenvalue pair needs no special treatment.

    Repeated eigenvalue: if A = c I then P = I, D = c I (fine). If A has a
    repeated eigenvalue but only one eigenvector (a "defective" matrix such as
    [[c, 1], [0, c]]) NO diagonalization exists. Then, and also when the two
    computed eigenvectors are numerically parallel (|det P| <= DEFECTIVE_DET_TOL),
    this function raises `DefectiveMatrixError` instead of returning a useless
    P^-1. Fallback for such matrices: `expm_closed_form`, which needs no
    eigenvectors (`expm_diag` switches to it automatically).

    Course topic: Diagonalization
    """
    lam, _ = eig2x2(A)
    P = eigvec2x2(A, lam).astype(complex)
    detP = P[0, 0] * P[1, 1] - P[0, 1] * P[1, 0]
    if abs(detP) <= DEFECTIVE_DET_TOL:
        raise DefectiveMatrixError(
            "matrix is (numerically) defective: its eigenvectors are parallel, "
            f"|det V| = {abs(detP):.3e}; use expm_closed_form instead")
    D = np.array([[lam[0], 0.0], [0.0, lam[1]]], dtype=complex)
    return P, D, inv2x2(P)


def is_stable2x2(A) -> bool:
    """True if both eigenvalues of the real 2x2 matrix A have negative real part.

    Formula (Routh-Hurwitz test for lam^2 - tr lam + det = 0):
        stable  <=>  tr(A) < 0  and  det(A) > 0.
    No eigenvalue has to be computed.

    Course topic: Stability
    """
    M = _real2x2(A)
    tr = M[0, 0] + M[1, 1]
    det = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
    return bool(tr < 0.0 and det > 0.0)


def time_constants(lam) -> tuple[float, float]:
    """Time constants of the two modes, returned as (tau_slow, tau_fast).

    Formula: a mode e^{lam t} decays like e^{-t/tau} with tau = -1 / Re(lam).
    tau_slow >= tau_fast. A mode that does not decay (Re(lam) >= 0) gets
    tau = inf. For a complex pair both time constants are equal.

    Course topic: Stability
    """
    lam = np.asarray(lam, dtype=complex).ravel()
    taus = [(-1.0 / z.real) if z.real < 0.0 else float("inf") for z in lam]
    return float(max(taus)), float(min(taus))


# ---------------------------------------------------------------------------
# Matrix exponential and step response
# ---------------------------------------------------------------------------
def expm_closed_form(A, t) -> np.ndarray:
    """e^{At} of a real 2x2 matrix WITHOUT eigenvectors (works for defective A).

    Formula (Cayley-Hamilton: every function of a 2x2 matrix is c I + s (A - m I)):
        m = tr/2,  delta^2 = disc/4 = (tr^2 - 4 det)/4,
        disc > 0 :  e^{At} = e^{mt} [ cosh(delta t) I + (sinh(delta t)/delta) (A - m I) ]
        disc < 0 :  e^{At} = e^{mt} [ cos(omega t) I + (sin(omega t)/omega) (A - m I) ],  omega^2 = -disc/4
        disc = 0 :  e^{At} = e^{mt} [ I + t (A - m I) ]          (repeated eigenvalue)
    The hyperbolic case is evaluated with exp/expm1 of the two real eigenvalues so
    that nothing overflows or cancels. t may be a number (returns 2x2) or an array
    of times (returns shape (len(t), 2, 2)).

    Course topic: Linear ODE systems
    """
    M = _real2x2(A)
    tt = np.atleast_1d(np.asarray(t, dtype=float)).ravel()
    _, info = eig2x2(M)
    m = 0.5 * info["tr"]
    disc = info["disc"]
    if disc > 0.0:
        delta = 0.5 * np.sqrt(disc)
        e_hi = np.exp((m + delta) * tt)
        e_lo = np.exp((m - delta) * tt)
        c = 0.5 * (e_hi + e_lo)                                  # e^{mt} cosh(delta t)
        s = e_hi * (-np.expm1(-2.0 * delta * tt)) / (2.0 * delta)  # e^{mt} sinh(delta t)/delta
    elif disc < 0.0:
        omega = 0.5 * np.sqrt(-disc)
        c = np.exp(m * tt) * np.cos(omega * tt)
        s = np.exp(m * tt) * np.sin(omega * tt) / omega
    else:
        c = np.exp(m * tt)
        s = tt * np.exp(m * tt)
    N = M - m * np.array([[1.0, 0.0], [0.0, 1.0]])               # A - m I
    out = c[:, None, None] * np.array([[1.0, 0.0], [0.0, 1.0]]) + s[:, None, None] * N
    return out[0] if np.ndim(t) == 0 else out


def _projectors(M: np.ndarray):
    """Eigenvalues and spectral projectors R_i = p_i q_i^T of a real 2x2 matrix.

    p_i = column i of P, q_i^T = row i of P^-1, so that A = lam1 R_1 + lam2 R_2 and
    f(A) = f(lam1) R_1 + f(lam2) R_2. Returns (lam, [R_1, R_2]), or None when A is
    defective or kappa2(P) > EXPM_KAPPA_SWITCH (callers then use the closed form).

    Course topic: Diagonalization
    """
    try:
        P, D, Pinv = diagonalize(M)
    except DefectiveMatrixError:
        return None
    if kappa2(P) > EXPM_KAPPA_SWITCH:
        return None
    Rs = []
    for i in range(2):
        R = np.empty((2, 2), dtype=complex)          # outer product p_i q_i^T, written out
        for r in range(2):
            for c in range(2):
                R[r, c] = P[r, i] * Pinv[i, c]
        Rs.append(R)
    return np.array([D[0, 0], D[1, 1]]), Rs


def expm_diag(A, t) -> np.ndarray:
    """Matrix exponential e^{At} of a real 2x2 matrix through its diagonalization.

    Formula: A = P D P^-1  =>  e^{At} = P diag(e^{lam1 t}, e^{lam2 t}) P^-1
                                      = e^{lam1 t} p1 q1^T + e^{lam2 t} p2 q2^T,
    where p_i is column i of P and q_i^T is row i of P^-1 (complex arithmetic; the
    imaginary parts of a conjugate pair cancel, the real part is returned).
    Fallback: if A is defective (`diagonalize` raises) or its eigenvector matrix
    is ill-conditioned (kappa2(P) > EXPM_KAPPA_SWITCH), the eigenvector-free
    formula `expm_closed_form` is used, because the rounding error of
    P e^{Dt} P^-1 grows in proportion to kappa2(P).
    t may be a number (returns 2x2) or an array of times (returns (len(t), 2, 2)).

    Course topic: Diagonalization
    """
    M = _real2x2(A)
    spec = _projectors(M)
    if spec is None:
        return expm_closed_form(M, t)
    lam, Rs = spec
    tt = np.atleast_1d(np.asarray(t, dtype=float)).ravel()
    out = np.zeros((tt.size, 2, 2), dtype=complex)
    for i in range(2):
        out += np.exp(lam[i] * tt)[:, None, None] * Rs[i]
    out = out.real
    return out[0] if np.ndim(t) == 0 else out


def step_response(A, B, du, t) -> np.ndarray:
    """Response dx(t) of dx/dt = A dx + B du to a constant step du applied at t = 0.

    Formula (solution of the linear ODE system with dx(0) = 0):
        dx(t) = A^-1 ( e^{At} - I ) B du.
    With the diagonalization A = P D P^-1 this is evaluated mode by mode,
        dx(t) = sum_i  (e^{lam_i t} - 1)/lam_i * R_i B du,   R_i = p_i q_i^T,
    which is the same formula (A^-1 R_i = R_i / lam_i) but gives dx(0) = 0 exactly.
    For a defective / ill-conditioned A (see `expm_diag`) the first form is used
    with `expm_closed_form` and the own 2x2 inverse. As t -> inf (stable A) the
    result tends to the steady gain -A^-1 B du.
    B may be (2,), (2,1) or (2,m); du a scalar or a length-m vector; t an array
    of times. Returns an array of shape (len(t), 2).

    Course topic: Linear ODE systems
    """
    M = _real2x2(A)
    Ainv = inv2x2(M)
    g = _input_vector(B, du)
    tt = np.atleast_1d(np.asarray(t, dtype=float)).ravel()
    spec = _projectors(M)
    if spec is not None:
        lam, Rs = spec
        dxc = np.zeros((tt.size, 2), dtype=complex)
        for i in range(2):
            phi = np.expm1(lam[i] * tt) / lam[i]     # (e^{lam t} - 1) / lam
            dxc[:, 0] += phi * (Rs[i][0, 0] * g[0] + Rs[i][0, 1] * g[1])
            dxc[:, 1] += phi * (Rs[i][1, 0] * g[0] + Rs[i][1, 1] * g[1])
        return dxc.real.copy()
    eAt = expm_closed_form(M, tt)                    # (len(t), 2, 2)
    # w = (e^{At} - I) g
    w0 = (eAt[:, 0, 0] - 1.0) * g[0] + eAt[:, 0, 1] * g[1]
    w1 = eAt[:, 1, 0] * g[0] + (eAt[:, 1, 1] - 1.0) * g[1]
    dx = np.empty((tt.size, 2))
    dx[:, 0] = Ainv[0, 0] * w0 + Ainv[0, 1] * w1
    dx[:, 1] = Ainv[1, 0] * w0 + Ainv[1, 1] * w1
    return dx


def steady_gain(A, B, du) -> np.ndarray:
    """Steady-state deviation after a constant input step du.

    Formula: 0 = A dx_ss + B du  =>  dx_ss = -A^-1 B du   (own 2x2 inverse).
    Returns a length-2 array.

    Course topic: Linear ODE systems
    """
    Ainv = inv2x2(_real2x2(A))
    g = _input_vector(B, du)
    return -np.array([Ainv[0, 0] * g[0] + Ainv[0, 1] * g[1],
                      Ainv[1, 0] * g[0] + Ainv[1, 1] * g[1]])


# ---------------------------------------------------------------------------
# Gauss-Seidel
# ---------------------------------------------------------------------------
def gauss_seidel(A, b, x0=None, tol: float = 1e-10, max_iter: int = 500
                 ) -> tuple[np.ndarray, list, dict]:
    """Solve A x = b by Gauss-Seidel iteration (own loop, any n x n system).

    Formula, for i = 1..n in order, always using the NEWEST available values:
        x_i^(k+1) = ( b_i - sum_{j<i} a_ij x_j^(k+1) - sum_{j>i} a_ij x_j^(k) ) / a_ii.
    Stops when the largest change of any component in one sweep is below `tol`
    (or after `max_iter` sweeps, or if the iterates stop being finite).

    Returns (x, history, info). history[0] is the start vector x0 (zeros if None)
    and history[k] the iterate after k sweeps. info = {'converged', 'iters',
    'last_change', 'residual_inf'} with residual_inf = max_i |b - A x|_i.

    Course topic: Gauss-Seidel
    """
    M = np.asarray(A, dtype=float)
    rhs = np.asarray(b, dtype=float).ravel()
    n = rhs.size
    if M.shape != (n, n):
        raise ValueError(f"A has shape {M.shape}, b has length {n}")
    for i in range(n):
        if M[i, i] == 0.0:
            raise ValueError(f"zero diagonal entry a[{i},{i}]: Gauss-Seidel cannot divide by it")
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float).ravel()
    history = [x.copy()]
    converged = False
    change = float("inf")
    iters = 0
    for k in range(1, max_iter + 1):
        change = 0.0
        for i in range(n):
            s = rhs[i]
            for j in range(n):
                if j != i:
                    s -= M[i, j] * x[j]              # x[j] is already updated for j < i
            new = s / M[i, i]
            change = max(change, abs(new - x[i]))
            x[i] = new
        history.append(x.copy())
        iters = k
        if not np.all(np.isfinite(x)):
            break
        if change < tol:
            converged = True
            break
    res = 0.0
    for i in range(n):
        s = rhs[i]
        for j in range(n):
            s -= M[i, j] * x[j]
        res = max(res, abs(s))
    info = {"converged": converged, "iters": iters,
            "last_change": float(change), "residual_inf": float(res)}
    return x, history, info


def gs_spectral_radius_2x2(A) -> float:
    """Spectral radius of the Gauss-Seidel iteration matrix of a 2x2 system.

    Formula: the iteration matrix -(D + L)^-1 U of a 2x2 matrix has eigenvalues
    0 and a12 a21 / (a11 a22), so
        rho_GS = | a12 a21 / (a11 a22) |.
    Gauss-Seidel converges if and only if rho_GS < 1; each sweep multiplies the
    error by about rho_GS.

    Course topic: Gauss-Seidel
    """
    M = _real2x2(A)
    if M[0, 0] == 0.0 or M[1, 1] == 0.0:
        raise ValueError("zero diagonal entry")
    return float(abs(M[0, 1] * M[1, 0] / (M[0, 0] * M[1, 1])))


def is_diag_dominant(A) -> bool:
    """True if the square matrix A is STRICTLY diagonally dominant by rows.

    Formula: |a_ii| > sum_{j != i} |a_ij| for every row i. This is a sufficient
    (not necessary) condition for Gauss-Seidel to converge.

    Course topic: Gauss-Seidel
    """
    M = np.asarray(A, dtype=float)
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError("square matrix expected")
    n = M.shape[0]
    for i in range(n):
        off = 0.0
        for j in range(n):
            if j != i:
                off += abs(M[i, j])
        if not abs(M[i, i]) > off:
            return False
    return True


# ---------------------------------------------------------------------------
# Finite-difference Jacobians
# ---------------------------------------------------------------------------
def fd_steps(x0, h: float, rel: bool = True, floor: float = 0.0) -> np.ndarray:
    """Perturbation sizes used by the finite-difference Jacobians.

    Formula: rel=True  -> step_i = h * max(|x0_i|, floor)   (h is a RELATIVE step);
             rel=False -> step_i = h                         (h is an absolute step).
    If x0_i = 0 and floor = 0 the relative rule would give a zero step, so
    step_i = h is used for that component.

    Course topic: Jacobians and partial derivatives
    """
    x0 = np.asarray(x0, dtype=float).ravel()
    if not rel:
        return np.full(x0.shape, float(h))
    scale = np.maximum(np.abs(x0), float(floor))
    scale = np.where(scale > 0.0, scale, 1.0)
    return float(h) * scale


def fd_jacobian(func: Callable, x0, h: float, scheme: str = "forward",
                rel: bool = True, floor: float = 0.0) -> np.ndarray:
    """Finite-difference Jacobian J[i, j] = d f_i / d x_j of func: R^n -> R^m.

    Formulas (one column per perturbed variable, e_j = unit vector, s_j = step):
        forward : J[:, j] = ( f(x0 + s_j e_j) - f(x0) ) / s_j              error O(h)
        central : J[:, j] = ( f(x0 + s_j e_j) - f(x0 - s_j e_j) ) / (2 s_j)  error O(h^2)
    Steps s_j come from `fd_steps` (relative by default). The divisor is the
    step that was ACTUALLY applied in floating point, (x0_j + s_j) - x0_j.
    Returns an array of shape (m, n).

    Course topic: Jacobians and partial derivatives
    """
    if scheme not in ("forward", "central"):
        raise ValueError("scheme must be 'forward' or 'central'")
    x0 = np.array(x0, dtype=float).ravel()
    steps = fd_steps(x0, h, rel, floor)
    f0 = np.asarray(func(x0.copy()), dtype=float).ravel() if scheme == "forward" else None
    cols = []
    for j in range(x0.size):
        xp = x0.copy()
        xp[j] = x0[j] + steps[j]
        fp = np.asarray(func(xp), dtype=float).ravel()
        if scheme == "forward":
            cols.append((fp - f0) / (xp[j] - x0[j]))
        else:
            xm = x0.copy()
            xm[j] = x0[j] - steps[j]
            fm = np.asarray(func(xm), dtype=float).ravel()
            cols.append((fp - fm) / (xp[j] - xm[j]))
    J = np.empty((cols[0].size, x0.size))
    for j, col in enumerate(cols):
        J[:, j] = col
    return J


def richardson_jacobian(func: Callable, x0, h: float, levels: int = 2,
                        rel: bool = True, floor: float = 0.0
                        ) -> tuple[np.ndarray, np.ndarray]:
    """Central-difference Jacobian improved by Richardson extrapolation.

    Formula: the central difference has the error expansion
        D(h) = J + c2 h^2 + c4 h^4 + ...     (Taylor series, even powers only),
    so combining two step sizes cancels the leading term:
        T[k][0] = D(h / 2^k),
        T[k][j] = T[k][j-1] + ( T[k][j-1] - T[k-1][j-1] ) / (4^j - 1).
    levels = 2 uses h and h/2:  J = (4 D(h/2) - D(h)) / 3, error O(h^4);
    levels = 3 adds h/4, error O(h^6).

    Returns (J, err_est). err_est has the same shape as J and equals
    |T[L-1][L-1] - T[L-1][L-2]|, the size of the last correction. While
    truncation dominates it is the error of the previous (lower-order) estimate,
    so it OVER-estimates the error of J; for very small h it shows the rounding
    noise instead. Use max(err_est) as a single number. (levels = 1: plain
    central difference, err_est = nan.)

    Course topic: Taylor series (truncation error)
    """
    if levels < 1:
        raise ValueError("levels must be >= 1")
    table = []
    for k in range(levels):
        row = [fd_jacobian(func, x0, h / 2.0 ** k, "central", rel, floor)]
        for j in range(1, k + 1):
            row.append(row[j - 1] + (row[j - 1] - table[k - 1][j - 1]) / (4.0 ** j - 1.0))
        table.append(row)
    J = table[-1][-1]
    if levels == 1:
        return J, np.full(J.shape, np.nan)
    return J, np.abs(J - table[-1][-2])


def complex_step_jacobian(func: Callable, x0, h: float = 1e-30) -> np.ndarray:
    """Jacobian by the complex-step method (for complex-safe analytic functions only).

    Formula (Taylor series in the imaginary direction):
        f(x + i h e_j) = f(x) + i h df/dx_j - h^2/2 d2f/dx_j^2 - ...
        =>  J[:, j] = Im( f(x0 + i h e_j) ) / h + O(h^2).
    There is NO subtraction of nearly equal numbers, so h can be tiny (1e-30)
    and the result is accurate to machine precision. `func` must accept complex
    arrays and use only analytic operations (no abs, no comparisons on values).
    h is an absolute step.

    Course topic: Jacobians and partial derivatives
    """
    x0 = np.asarray(x0, dtype=float).ravel()
    cols = []
    for j in range(x0.size):
        xc = x0.astype(complex)
        xc[j] = xc[j] + 1j * h
        cols.append(np.imag(np.asarray(func(xc)).ravel()) / h)
    J = np.empty((cols[0].size, x0.size))
    for j, col in enumerate(cols):
        J[:, j] = col
    return J


def optimal_h_forward(eps: float, M2: float, scale: float = 1.0) -> float:
    """Step that minimises the total error of a FORWARD difference.

    Formula: total error(h) ~ (M2 / 2) h  +  2 e_f / h
             (Taylor truncation)     (rounding / solver noise in f),
    minimised at  h* = 2 sqrt(e_f / M2),  with e_f = eps * scale.
    eps = relative noise of one function value (machine epsilon, or the solver
    tolerance), scale = typical |f| (scale = 1 means eps is already an absolute
    noise), M2 = bound on |f''|. h* is an ABSOLUTE step in x; divide by |x0| to
    get the relative step used by `fd_jacobian`. Returns inf if M2 = 0.

    Course topic: Taylor series (truncation error)
    """
    if eps < 0.0 or M2 < 0.0 or scale < 0.0:
        raise ValueError("eps, M2 and scale must be non-negative")
    if M2 == 0.0:
        return float("inf")
    return float(2.0 * np.sqrt(eps * scale / M2))


def optimal_h_central(eps: float, M3: float, scale: float = 1.0) -> float:
    """Step that minimises the total error of a CENTRAL difference.

    Formula: total error(h) ~ (M3 / 6) h^2  +  e_f / h,
    minimised at  h* = (3 e_f / M3)^(1/3),  with e_f = eps * scale.
    eps, scale as in `optimal_h_forward`; M3 = bound on |f'''|. h* is an
    ABSOLUTE step in x. Returns inf if M3 = 0.

    Course topic: Taylor series (truncation error)
    """
    if eps < 0.0 or M3 < 0.0 or scale < 0.0:
        raise ValueError("eps, M3 and scale must be non-negative")
    if M3 == 0.0:
        return float("inf")
    return float((3.0 * eps * scale / M3) ** (1.0 / 3.0))


# ---------------------------------------------------------------------------
# Bauer-Fike: eigenvalue sensitivity
# ---------------------------------------------------------------------------
def bauer_fike_bound(V, E) -> float:
    """Bauer-Fike upper bound on how far a perturbation E can move the eigenvalues.

    Theorem: if A = V D V^-1 is diagonalizable and mu is any eigenvalue of A + E,
    then   min_j |mu - lam_j|  <=  kappa_2(V) * ||E||_2.
    This function returns the right-hand side (own kappa2 and own 2-norm).

    Course topic: Matrix norms and conditioning
    """
    return float(kappa2(V) * norm2(E))


def random_perturbations(rng: np.random.Generator, n: int, norm2: float) -> np.ndarray:
    """n random real 2x2 matrices E, each rescaled so that ||E||_2 == norm2.

    Formula: draw the four entries from a standard normal distribution, then
    E <- E * norm2 / sigma_max(E) (own closed-form sigma_max). The rescaling is
    applied twice so that the norm is exact to the last floating-point digit.
    `rng` is a numpy Generator (e.g. numpy.random.default_rng(seed)), so the
    result is reproducible. Returns an array of shape (n, 2, 2).

    Course topic: Matrix norms and conditioning
    """
    if n < 0 or norm2 < 0.0:
        raise ValueError("n and norm2 must be non-negative")
    E = rng.standard_normal((int(n), 2, 2))
    for _ in range(2):
        s = norm2_batch(E)
        s = np.where(s > 0.0, s, 1.0)
        E = E * (float(norm2) / s)[:, None, None]
    return E


def eig_shift(A, E) -> float:
    """How far the eigenvalues of A move when A is replaced by A + E.

    Formula: shift = max_i min_j | mu_i - lam_j |,  lam = eig(A), mu = eig(A + E),
    i.e. the largest distance from a perturbed eigenvalue to the NEAREST original
    one. This is exactly the quantity bounded by Bauer-Fike
    (shift <= kappa_2(V) ||E||_2).

    Course topic: Eigenvalues
    """
    lam, _ = eig2x2(A)
    mu, _ = eig2x2(_real2x2(A) + _real2x2(E))
    shift = 0.0
    for i in range(2):
        nearest = min(abs(mu[i] - lam[0]), abs(mu[i] - lam[1]))
        shift = max(shift, nearest)
    return float(shift)


def eig_shift_batch(A, Es) -> np.ndarray:
    """`eig_shift(A, E)` for a whole batch of perturbations Es (shape (n, 2, 2)).

    Formula: shift_k = max_i min_j | mu_i(A + E_k) - lam_j(A) |, with the
    closed-form 2x2 eigenvalues evaluated for all k at once (`eig2x2_batch`).
    Returns an array of length n. Fast enough for Monte Carlo studies
    (about a million 2x2 eigenvalue problems per second).

    Course topic: Eigenvalues
    """
    M = _real2x2(A)
    lam, _ = eig2x2(M)
    mu = eig2x2_batch(M[None, :, :] + np.asarray(Es, dtype=float))     # (n, 2)
    d0 = np.abs(mu - lam[0])                                           # distance to lam_1
    d1 = np.abs(mu - lam[1])                                           # distance to lam_2
    nearest = np.minimum(d0, d1)                                       # (n, 2)
    return np.maximum(nearest[:, 0], nearest[:, 1])


__all__ = [
    "DEFECTIVE_DET_TOL", "EXPM_KAPPA_SWITCH", "DefectiveMatrixError",
    "eig2x2", "eig2x2_batch", "eigvec2x2", "inv2x2", "singular_values_2x2", "norm2",
    "norm2_batch", "kappa2", "diagonalize", "is_stable2x2", "time_constants",
    "expm_closed_form", "expm_diag", "step_response", "steady_gain",
    "gauss_seidel", "gs_spectral_radius_2x2", "is_diag_dominant",
    "fd_steps", "fd_jacobian", "richardson_jacobian", "complex_step_jacobian",
    "optimal_h_forward", "optimal_h_central",
    "bauer_fike_bound", "random_perturbations", "eig_shift", "eig_shift_batch",
]
