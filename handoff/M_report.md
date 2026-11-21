# M_report - math core and published-example anchor (agent M)

Status: all deliverables built; `python -m tests.test_mathcore` 68 PASS / 0 FAIL (about 4 s);
`python -m tests.test_anchors` 38 PASS / 0 FAIL; `python -m eigenmap.anchors` 64 PASS / 0 FAIL
(fast and full). `outputs/data/anchors.json` (full config) is byte-identical over two runs,
sha256 `bc90f8d1f75b5c37e0291cefac0f100651a987116ad1271991b9d569d3665c29`; the fast file is
also identical between an in-process run and a separate Python process.

## Files (only these were written)

- `eigenmap/mathcore.py` - own 2x2 linear algebra + finite differences (numpy only for storage,
  + - * /, sqrt, exp/expm1, sin, cos; no numpy.linalg / scipy inside).
- `eigenmap/anchors.py` - module `anchors`, `run(cfg, out_path=None) -> dict`, writes
  `outputs/data/anchors.json` with a `checks` list.
- `tests/test_mathcore.py`, `tests/test_anchors.py`, this report. (No `tests/__init__.py` was
  created; `python -m tests.<name>` works as a namespace package.)

## mathcore: what is in it

Every CONTRACT section 5 function exists with the listed name and signature. Extras that other
agents may use: `eig2x2_batch(As)`, `eig_shift_batch(A, Es)`, `inv2x2`, `singular_values_2x2`,
`norm2`, `norm2_batch`, `expm_closed_form`, `fd_steps`, `DefectiveMatrixError`.

Interface decisions where the contract left room (please confirm or tell me to change):

1. `fd_jacobian(..., rel=True, floor=0.0)`: step_i = h*max(|x0_i|, floor). Default floor is 0, so a
   fuel flow below 1 lbm/s still gets step h*|u0| as the notation table requires. If x0_i = 0 and
   floor = 0 the step is h (absolute). Divisor is the step actually applied in floating point.
2. `richardson_jacobian(func, x0, h, levels=2)` uses central differences at h, h/2 (, h/4 ...).
   `err_est` is an ARRAY of the shape of J (size of the last Richardson correction); take
   `err_est.max()` for one number. It over-estimates the error while truncation dominates and
   shows the rounding noise for very small h.
3. `optimal_h_forward(eps, M2, scale)` = 2 sqrt(eps*scale/M2), `optimal_h_central` =
   (3 eps*scale/M3)^(1/3): `scale` = typical |f| (eps*scale = absolute noise of one f value).
   The returned h is an ABSOLUTE step in x; divide by |x0| for the relative h of `fd_jacobian`.
4. `gauss_seidel` history: `history[0]` is the start vector, `history[k]` the iterate after k
   sweeps; stop rule = largest component change in a sweep < tol.
5. `eig2x2` accepts REAL matrices only. `info['disc']` is evaluated as (a11-a22)^2 + 4 a12 a21
   (equal to tr^2 - 4 det, without cancellation). Real ties in |Re| put the larger value first.
6. `time_constants` returns inf for a mode with Re(lambda) >= 0.
7. `step_response` evaluates A^-1 (e^{At} - I) B du mode by mode, sum_i (e^{lam_i t}-1)/lam_i R_i B du
   (same formula, dx(0) = 0 exactly); B may be (2,), (2,1) or (2,m).

Numerical caveats (near-defective matrices):

- `diagonalize` raises `DefectiveMatrixError` when |det V| <= 1e-10 (unit-norm columns), i.e. for a
  repeated eigenvalue with one eigenvector. A = c I returns P = I.
- `expm_diag` (and `step_response`) switch to the eigenvector-free Cayley-Hamilton formula
  `expm_closed_form` when A is defective or kappa2(V) > 1e3, because the rounding error of
  P e^{Dt} P^-1 grows like eps*kappa2(V). 6 of the 262 test matrices take that path.
- Eigenvalues of nearly defective matrices are ill-conditioned for ANY algorithm (error up to
  eps*kappa(V), sqrt(eps) at an exact double root). The test therefore uses 1e-10 for kappa(V) < 1e3
  (256 matrices, worst 3.6e-16), max(1e-10, 50 eps kappa(V)) for 1e3 <= kappa(V) < 1e7 (4 matrices)
  and 1e-7 for kappa(V) >= 1e7 (2 Jordan blocks, observed difference 0).
- `kappa2` / `inv2x2` comparisons are scaled by max(1, kappa/1e3) for the same reason.

## test_mathcore: 68 checks

240 seeded random stable 2x2 matrices (seed 20261114) + 22 edge cases (complex pairs, triangular,
diagonal, scalar, singular, zero, stiff, near-defective with gaps 2e-3 .. 2e-7, two Jordan blocks,
the four anchor matrices). Worst observed errors: eigenvalues 3.6e-16, eigenvectors 4.4e-16,
P D P^-1 1e-15, kappa 6e-14, expm 2.5e-12, step response vs exact augmented expm 8e-14, vs
solve_ivp 1.3e-11, Gauss-Seidel 2.7e-14 (120 systems, n = 2..6), complex step 2.2e-16.
FD orders measured: forward 0.995, central 2.000, Richardson 4.002 (levels 2), 5.989 (levels 3).
Optimal h on f = exp: forward predicted 2.98e-8 vs empirical 1.26e-8; central 8.73e-6 vs 5.01e-6
(both within half a decade; error at the predicted h within 3x of the sweep minimum).
Bauer-Fike: 0 violations in 5000 random cases (1570 complex pairs, 20 real<->complex flips),
largest shift/bound 0.985. `eig_shift_batch`: 20,000 samples x 20 points in 0.03 s.

Changes made to checks after the first run (all test-calibration, no library bug found):

- solve_ivp comparison: the 2e-8 difference was scipy's dense-output interpolation (`t_eval`), not
  the own routine (it agrees with the exact augmented matrix exponential to 8e-14). The test now
  steps the integrator exactly to each output time (rtol 1e-13) and passes at 1e-10. Systems with
  |lambda|_max * t_end > 2000 (2 stiff cases) are compared with the exact reference only.
- Richardson err_est: one Jacobian entry of the test function is exactly linear, so its error is
  pure rounding (1e-16) with err_est 0; a 1e-13 rounding floor was added and named in the check.
- Richardson vs analytic: at h = 1e-2 the O(h^4) truncation error itself is 1.8e-10; the check now
  uses h = 1e-3 (error 9.5e-14).
- "defective" test-name matching bug in the test; dx(0) now exactly 0 (code improvement, item 7).

## anchors: turboshaft JetCat SPT5 (illustration) - reproduced numbers

Every label, console line and check name carries "turboshaft JetCat SPT5 (illustration)"; the JSON
says it is NOT the AGTF30 turbofan. Matrices were transcribed from the brief (the paper was not
re-opened: no network).

| quantity | idle | mid | cruise | brief |
|---|---|---|---|---|
| lambda | -0.34614 / -0.37386 | -0.59917 / -0.89083 | -1.01270 / -1.78730 | -0.346/-0.374, -0.599/-0.891, -1.013/-1.787 |
| tau (model time units) | 2.88897 / 2.67482 | 1.66896 / 1.12255 | 0.98746 / 0.55950 | 2.89/2.67, 1.67/1.12, 0.99/0.56 |
| kappa(V) | 18.8746 | 2.3929 | 1.8357 | 18.9 / 2.39 / 1.84 |
| Bauer-Fike bound, norm 1e-3 | 0.018875 | 0.002393 | 0.001836 | - |
| MC worst shift, n = 20000 | 0.016171 / 0.016159 / 0.016156 | 0.001409 (all seeds) | 0.001191 (all seeds) | 0.016 / - / 0.0012 |

- MC idle/cruise ratio per seed: 13.58 / 13.57 / 13.57 (brief about 13x; allowed 10..16). Fast
  config (n = 500): idle 0.01606 / 0.01610 / 0.01613, cruise 0.00119, ratio 13.49 / 13.55 / 13.56.
  About 10 % of the perturbed idle matrices have a complex pair; none at mid or cruise.
- Gauss-Seidel (A_cruise x = [-1.2, -0.3], x0 = 0): 0.70588/0.65775, 0.74457/0.67886,
  0.74582/0.67954, 0.74586/0.67956; exact 0.74586/0.67956; spectral radius 0.03209; 8 sweeps to 1e-10.
- Cruise unit-fuel-step steady state [0.74586, 0.67956].
- Flip: a12 -0.0008 -> -0.0018 at idle gives -0.36000 +/- 0.00825i (disc 7.68e-4 -> -2.72e-4).

Nothing in the brief failed to reproduce. Two things to know:

1. THREE tau values meet the brief only to its printed rounding, not the contract's 2e-3 absolute:
   idle tau_fast 2.67482 vs 2.67 (diff 4.8e-3), mid tau_fast 1.12255 vs 1.12 (2.6e-3), cruise
   tau_slow 0.98746 vs 0.99 (2.5e-3). A 2-decimal brief value cannot be matched to 2e-3 in general.
   These three checks are printed as PASS with the text "matches brief ... ONLY to its printed
   rounding" and limit "0.005 (printed rounding; contract 2e-3 NOT met)", and are listed in the JSON
   key `matched_to_printed_rounding_only`. If you prefer FAIL, or want the brief quoted to 3 d.p.
   (2.675, 1.123, 0.987), say so; it is a one-line change.
2. Flip wording: the change of 1e-3 in a12 is below the rounding half-width (5e-3) of the three
   entries printed to 2 decimals, but it is LARGER than the rounding half-width (5e-5) of a12 itself
   as printed. Over the whole box of matrices that round to the printed one the discriminant goes
   down to about -1e-6 (extreme corner only), so the printed digits only just fix the idle pair as
   real. Stored in `flip.note`; the report text should not say the flip lies inside the printed
   precision of a12.

Other notes:

- The Monte Carlo uses `cfg.E_norm`; the comparison with the brief is made only when it equals 1e-3
  (otherwise one "skipped" check is recorded). The same perturbations are applied to all three
  points for a given seed.
- The step-response grid (0..20 model time units, 201 points) is a plotting grid defined in
  `anchors.py`, not in `config.py`.
- JSON keys: `points.<idle|mid|cruise>` (A, B, eig own/numpy/brief, tau, V, kappaV, steady_gain,
  step_fuel t/dx1/dx2), `gauss_seidel.table` (k, x1, x2, err_inf, err_ratio, brief), `flip`,
  `monte_carlo.per_seed` (worst/mean shift, quantiles, fraction_complex, worst_over_bound),
  `monte_carlo.bauer_fike_bound`, `checks`.
