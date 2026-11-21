# CODE_GUIDE.md — where each course topic lives in the code

Every function's docstring ends with a line `Course topic: ...`. Read the modules in this order.

| Step | Module | What it does | Course topic |
|---|---|---|---|
| 1 | `eigenmap/mathcore.py` | Our OWN implementations: 2x2 eigenvalues from the characteristic equation, eigenvectors, diagonalization A = P D P^-1, e^{At}, step response, steady gain, Gauss-Seidel, finite-difference / Richardson / complex-step Jacobians, kappa(V), Bauer-Fike tools. Each is cross-checked against numpy/scipy in `tests/test_mathcore.py`. | Eigenvalues, diagonalization, linear ODE systems, stability, Gauss-Seidel, Jacobians, Taylor series |
| 2 | `eigenmap/anchors.py` | The published turboshaft example (JetCat SPT5, illustration only): eigenvalues, tau, kappa(V), Gauss-Seidel table, Monte Carlo. | Eigenvalues, Gauss-Seidel, conditioning |
| 3 | `eigenmap/agtf30_data.py`, `tmats_thermo.py`, `engine_model.py`, `steady_state.py` | The NASA AGTF30 engine: reads NASA's map files, gas properties, components, the flow-matching residuals r(z) = 0, shaft equations I (pi/30) dN/dt = torque, Newton steady-state solver. | Engine model (not a course topic); Newton's method uses a Jacobian |
| 4 | `eigenmap/linearize.py` | Taylor expansion about F(x0, u0) = 0: A = dF/dx, B = dF/du by finite differences; reference = central + Richardson. | Partial derivatives, Jacobians, Taylor series |
| 5 | `eigenmap/opoints.py` | Chooses the steady operating points (idle to max, sea level and 35,000 ft). | — |
| 6 | `eigenmap/layer1_opmap.py` | Layer 1: eigenvalues, tau, kappa(V), tr/det/discriminant vs power; hypothesis H1. | Eigenvalues, stability |
| 7 | `eigenmap/steps.py` | Step response dx(t) = A^-1 (e^{At} - I) B du and the new steady state by Gauss-Seidel, at every point. | Linear ODE systems, diagonalization, Gauss-Seidel |
| 8 | `eigenmap/layer2_fd.py` | Layer 2: eigenvalue error vs finite-difference step h and solver tolerance tol; hypothesis H3. | Jacobians, Taylor series (truncation error) |
| 9 | `eigenmap/layer3_bauerfike.py` | Layer 3: Monte Carlo check of the Bauer-Fike theorem; hypothesis H2. | Eigenvalues, matrix norms, conditioning |
| 10 | `eigenmap/layer4_degradation.py` | Layer 4: eigenvalue / tau shifts when components wear; hypothesis H4. | Eigenvalue sensitivity |
| 11 | `eigenmap/validity.py` | Linear vs nonlinear response to fuel steps of 1-25 %. | Taylor series (range of validity), linear ODE systems |
| 12 | `eigenmap/validation.py` | Port vs NASA Table 1; Gate-1 summary. | Model validation |
| 13 | `eigenmap/evidence/build.py` | Figures, tables, RESULTS_FACTS.md, objectives_trace.md from the JSON. | — |
| 14 | `twin/` | Virtual lab extension (not part of the report objectives): `sim.py` live Euler simulation of the 2-state ODE with NASA's baseline controller; `experiments.py` the seven experiments and 16 presets; `app.py` the UI; `assets/engine_canvas.js` the animated engine; `export_viewer.py` offline replay. Tested by `tests/test_lab.py`. | ODE initial-value problems (Euler), equilibria, linearization, eigenvalues, iterative methods |

## 15 viva questions with short model answers

Numbers: take them from `RESULTS_FACTS.md` (full run), never from memory.

1. **Why can we linearize a nonlinear engine?** Near a steady point F(x0, u0) = 0 the first-order Taylor expansion F(x0+dx, u0+du) ≈ A dx + B du is accurate when dx, du are small; the neglected terms are second order. Our validity study measures how small "small" is (linear error ≤ 2 % for a 1 % fuel step).
2. **What are A and B physically?** a11 = how the LP-shaft acceleration changes with LP speed (negative: a faster fan needs more torque, which slows it), a12 = effect of core speed on the LP shaft, etc. B = how much each shaft accelerates per lbm/s of extra fuel.
3. **Why are the diagonal entries negative?** Each spool's load torque rises faster with its own speed than its turbine torque does, so a speed increase produces a decelerating torque: self-stabilizing. All our points confirm it.
4. **How do you get eigenvalues of a 2x2 matrix by hand?** Solve λ² − tr(A) λ + det(A) = 0: λ = (tr ± √(tr² − 4 det))/2. Our code does exactly this and checks it against numpy.
5. **When is the linear model stable?** When both eigenvalues have negative real part; for 2x2 this is equivalent to tr(A) < 0 and det(A) > 0.
6. **What is a time constant?** A mode e^{λt} decays to 37 % after τ = −1/Re(λ). The slow mode dominates how long the engine takes to settle.
7. **What does a complex eigenvalue pair mean?** The two shafts exchange energy and the response has a small oscillation (Im λ) on top of the decay. The port has one at max cruise; the published turboshaft has one only after a tiny change of a12.
8. **Why diagonalize?** With A = P D P^-1, e^{At} = P e^{Dt} P^-1 and the coupled system splits into two independent first-order equations (modes), each solved by an exponential.
9. **Why Gauss-Seidel for a 2x2 system?** It is the course's iterative method; it converges here because the spectral radius |a12 a21/(a11 a22)| is far below 1, and we show it matches the direct solution −A^-1 B du. (A in rpm units is diagonally dominant at every point; that property depends on the scaling of the states, the spectral radius does not.)
10. **Why does the finite-difference error curve look like a U?** Large h: Taylor truncation error (∝ h for forward, h² for central). Small h: textbook answer, the noise in F divided by h. In our engine that classical rise appears only at tight tolerance (round-off). At loose tolerance the left side is different: a flat plateau of ~15–30 % error (the solver does no iterations) that ends in a cliff. So the minimum exists, but not for the textbook reason.
11. **What surprised you in Layer 2?** The solver error is not random noise. With a warm start, a tiny perturbation needs zero solver iterations, so F is evaluated with the flow variables frozen and the derivative is about 15–30 % wrong. The best step sits just above that dead zone; it grows with tol with a log-log slope of about 0.6–0.8 (the classical noise law predicts 0.5), and from a start that is not the solution it does not shift consistently at all.
12. **What does the Bauer–Fike theorem say?** An error E in A moves each eigenvalue by at most κ(V)·‖E‖₂. A large κ(V) (nearly parallel eigenvectors) means the eigenvalues are fragile. Our Monte Carlo never exceeded the bound.
13. **Is the slowest point the most fragile here?** (And which idle?) At sea-level idle the model has two steady states for the same settings; we use the one reached from max power and report both. No. In the AGTF30 port the slowest point is cruise idle, while the most fragile is near sea-level idle; κ(V) stays modest everywhere. The published turboshaft example behaves differently (κ ≈ 19 at idle). We report the contradiction rather than hide it.
14. **How does wear change the dynamics?** It depends on the point and on the comparison rule: at sea level, combined degradation slows the engine under both rules; at max cruise the sign differs between rules, because that point sits next to the real/complex boundary. Hence "mixed".
15. **How do you know the engine port is right?** (a) An independent blind re-implementation agrees to 1e-12 on residuals, torques and steady states; (b) airflow, TSFC and thrust are within 0.32 % of NASA's published Table 1 and TIT within 0.67 %; BPR and OPR (printed rounded to whole numbers) are off by up to 1.48 % (take-off BPR 27.6 vs 28), within the printed rounding; no NASA A/B matrices exist to compare the linear model itself; (c) 500+ automated checks run with `./run.sh verify`.
