# RESULTS_FACTS (generated from outputs/data/*.json - facts only, no report prose)

Engine: NASA AGTF30 conceptual N+3 geared turbofan (NOT a production engine), Python port of the T-MATS model.
Run settings: full fidelity settings.

## Validation of the port

- No NASA A/B matrices are available (nasa_lin.mat not supplied; none in the repo). Validation is of the steady-state model.
- Gate 1 (independent blind replica): 140 checks passed, 0 failed, worst relative disagreement 8.45e-13.
- Worst |relative error| vs NASA Table 1 (faithful port): 1.48 % (take-off BPR; Table 1 prints BPR rounded to integers).
- cruise: BPR 24 vs 24 (-0.01 %); Fnet 6086 vs 6073 (+0.22 %); OPR 55.03 vs 55 (+0.06 %); T4 3163 vs 3150 (+0.40 %); TSFC 0.4622 vs 0.4637 (-0.32 %); W_in 815.6 vs 813 (+0.32 %)
- takeoff: BPR 27.59 vs 28 (-1.48 %); Fnet 2.863e+04 vs 2.862e+04 (+0.04 %); OPR 36.64 vs 37 (-0.98 %); T4 3191 vs 3170 (+0.67 %); TSFC 0.1746 vs 0.1751 (-0.27 %); W_in 1724 vs 1724 (-0.02 %)
- OPR is compared as the product of the fan, LPC and HPC pressure ratios (ASSUMED NASA definition). As Pt3/Pt2 (duct losses included) it would be 53.72 and 35.99, i.e. -2.3 % and -2.7 % low.
- The study uses the differentiable 'study' mode of the same engine (C1 map interpolation, inner loops converged); it moves the anchor quantities by at most 0.27 % relative to the faithful port.

## O1/O2 - linear models and eigen-analysis (18 operating points + 2 NASA anchors)

- SLS_00 SLS idle (Nc_fan 892 rpm): lam = -1.2438, -3.0362 1/s; tau_slow 0.804 s, tau_fast 0.329 s; kappa(V) 6.59 (normalized), 1.30 (rpm); tr -4.280, det 3.776
- SLS_01 SLS part power 1 (Nc_fan 1005 rpm): lam = -1.3438, -3.3487 1/s; tau_slow 0.744 s, tau_fast 0.299 s; kappa(V) 6.87 (normalized), 1.24 (rpm); tr -4.693, det 4.500
- SLS_02 SLS part power 2 (Nc_fan 1118 rpm): lam = -1.4172, -3.3239 1/s; tau_slow 0.706 s, tau_fast 0.301 s; kappa(V) 6.70 (normalized), 1.29 (rpm); tr -4.741, det 4.711
- SLS_03 SLS part power 3 (Nc_fan 1232 rpm): lam = -1.4915, -3.5348 1/s; tau_slow 0.670 s, tau_fast 0.283 s; kappa(V) 6.27 (normalized), 1.28 (rpm); tr -5.026, det 5.272
- SLS_04 SLS part power 4 (Nc_fan 1345 rpm): lam = -1.6359, -4.0132 1/s; tau_slow 0.611 s, tau_fast 0.249 s; kappa(V) 5.47 (normalized), 1.20 (rpm); tr -5.649, det 6.565
- SLS_05 SLS part power 5 (Nc_fan 1458 rpm): lam = -1.8362, -4.7622 1/s; tau_slow 0.545 s, tau_fast 0.210 s; kappa(V) 4.59 (normalized), 1.08 (rpm); tr -6.598, det 8.744
- SLS_06 SLS part power 6 (Nc_fan 1572 rpm): lam = -2.0255, -5.3553 1/s; tau_slow 0.494 s, tau_fast 0.187 s; kappa(V) 4.12 (normalized), 1.05 (rpm); tr -7.381, det 10.847
- SLS_07 SLS part power 7 (Nc_fan 1685 rpm): lam = -2.2531, -5.6135 1/s; tau_slow 0.444 s, tau_fast 0.178 s; kappa(V) 3.94 (normalized), 1.08 (rpm); tr -7.867, det 12.648
- SLS_08 SLS part power 8 (Nc_fan 1798 rpm): lam = -2.5118, -5.3485 1/s; tau_slow 0.398 s, tau_fast 0.187 s; kappa(V) 4.17 (normalized), 1.21 (rpm); tr -7.860, det 13.434
- SLS_09 SLS max power (Nc_fan 1912 rpm): lam = -2.7942, -5.8918 1/s; tau_slow 0.358 s, tau_fast 0.170 s; kappa(V) 3.98 (normalized), 1.22 (rpm); tr -8.686, det 16.463
- CRUISE_00 35k ft / M0.8 idle (Nc_fan 1012 rpm): lam = -0.3339, -1.1922 1/s; tau_slow 2.995 s, tau_fast 0.839 s; kappa(V) 4.14 (normalized), 1.23 (rpm); tr -1.526, det 0.398
- CRUISE_01 35k ft / M0.8 part power 1 (Nc_fan 1205 rpm): lam = -0.5570, -1.5478 1/s; tau_slow 1.795 s, tau_fast 0.646 s; kappa(V) 5.38 (normalized), 1.16 (rpm); tr -2.105, det 0.862
- CRUISE_02 35k ft / M0.8 part power 2 (Nc_fan 1397 rpm): lam = -0.6703, -1.7580 1/s; tau_slow 1.492 s, tau_fast 0.569 s; kappa(V) 4.90 (normalized), 1.14 (rpm); tr -2.428, det 1.178
- CRUISE_03 35k ft / M0.8 part power 3 (Nc_fan 1590 rpm): lam = -0.8184, -2.1870 1/s; tau_slow 1.222 s, tau_fast 0.457 s; kappa(V) 4.33 (normalized), 1.09 (rpm); tr -3.005, det 1.790
- CRUISE_04 35k ft / M0.8 part power 4 (Nc_fan 1782 rpm): lam = -1.0024, -2.8136 1/s; tau_slow 0.998 s, tau_fast 0.355 s; kappa(V) 3.54 (normalized), 1.03 (rpm); tr -3.816, det 2.820
- CRUISE_05 35k ft / M0.8 part power 5 (Nc_fan 1975 rpm): lam = -1.1918, -3.1788 1/s; tau_slow 0.839 s, tau_fast 0.315 s; kappa(V) 2.97 (normalized), 1.10 (rpm); tr -4.371, det 3.789
- CRUISE_06 35k ft / M0.8 part power 6 (Nc_fan 2167 rpm): lam = -1.4421, -2.5155 1/s; tau_slow 0.693 s, tau_fast 0.398 s; kappa(V) 3.70 (normalized), 1.20 (rpm); tr -3.958, det 3.628
- CRUISE_07 35k ft / M0.8 max cruise (Nc_fan 2360 rpm): lam = -1.3967+0.1627i, -1.3967-0.1627i 1/s; tau_slow 0.716 s, tau_fast 0.716 s; kappa(V) 4.59 (normalized), 4.73 (rpm); tr -2.793, det 1.977; COMPLEX PAIR
- ANCHOR_CRUISE NASA cruise anchor (35k ft, M0.8): lam = -1.4007+0.1683i, -1.4007-0.1683i 1/s; tau_slow 0.714 s, tau_fast 0.714 s; kappa(V) 4.45 (normalized), 4.51 (rpm); tr -2.801, det 1.990; COMPLEX PAIR
- ANCHOR_TAKEOFF NASA max-thrust anchor (SLS, +27 F): lam = -2.7327, -5.6890 1/s; tau_slow 0.366 s, tau_fast 0.176 s; kappa(V) 4.04 (normalized), 1.24 (rpm); tr -8.422, det 15.547
- All points stable, all diagonal entries of A negative (see checks).
- Complex eigenvalue pairs (real-to-complex flip exists in THIS model): CRUISE_07.

### Multiple steady states at idle (found by the audit, confirmed in the faithful port)

- At low power the AGTF30 model (faithful port too, checked) has more than one steady state for the same corrected fan speed, VAFN and VBV. The study uses the branch reached by continuation from max power along NASA's VAFN/VBV schedules (what a decelerating engine follows); 'other' is the steady state found from NASA's default initial guess. Each operating point is reproducible from its tabulated z0 (warm start). The flow-matching equations themselves can also have more than one root at a fixed (x, u): from NASA's default start the inner solver can converge to a spurious root (e.g. LPT pressure ratio ~150 at SLS idle) that is not the engine's state; all F evaluations of the study therefore start from the operating point's own z0.
- CRUISE_00: no second steady state found from NASA's default initial guess.
- SLS_00: study branch W_f 0.2329 lbm/s, N_HP 17230 rpm, T4 1937 R; OTHER steady state with the same Nc_fan, VAFN, VBV: W_f 0.1540, N_HP 16542, T4 1723, lam -0.4531, -2.8112, tau_slow 2.207 s, kappa(V) 15.83, stable True
- H1 with the other idle branch instead: SLS tau_slow(idle)/tau_slow(max) = 6.17

### H1 - slower at low power: **supported**

- Rule (applied by code): supported if, in BOTH the SLS and the cruise sweep, Spearman rank correlation(power, tau_slow) < -0.7 AND tau_slow(idle)/tau_slow(max) > 1.5; not supported if the ratio < 1 in both; otherwise mixed
- Numbers: `{"CRUISE": {"ratio_idle_over_max": 4.18326155034772, "spearman_power_vs_tau_fast": -0.4285714285714286, "spearman_power_vs_tau_slow": -0.9761904761904763, "tau_fast_idle": 0.8388001047694377, "tau_fast_max": 0.7159859194450939, "tau_slow_idle": 2.9951563674050212, "tau_slow_max": 0.7159859194450939}, "SLS": {"ratio_idle_over_max": 2.246539487008533, "spearman_power_vs_tau_fast": -0.9515151515151514, "spearman_power_vs_tau_slow": -0.9999999999999999, "tau_fast_idle": 0.3293603317016725, "tau_fast_max": 0.16972604283324116, "tau_slow_idle": 0.8039894227160611, "tau_slow_max": 0.3578790523671785}}`

## O3 - finite-difference accuracy

- [converged start] SLS_00 tol 1e-13: best forward h 3.2e-08 (error 2.3e-07), theory 2.4e-09; best central h 3.2e-06 (error 3.8e-09), theory 4.2e-07; eps_F 1.8e-15
- [converged start] SLS_00 tol 1e-09: best forward h 3.2e-08 (error 8.5e-07), theory 5.2e-08; best central h 3.2e-06 (error 4.2e-09), theory 3.3e-06; eps_F 8.5e-13
- [converged start] SLS_00 tol 1e-06: best forward h 3.2e-07 (error 9.2e-06), theory 4.5e-06; best central h 1.0e-06 (error 4.4e-09), theory 6.5e-05; eps_F 6.4e-09
- [converged start] SLS_00 tol 0.0001: best forward h 5.6e-05 (error 4.9e-04), theory 4.5e-05; best central h 5.6e-05 (error 5.4e-07), theory 3.0e-04; eps_F 6.4e-07
- [converged start] CRUISE_00 tol 1e-13: best forward h 5.6e-08 (error 5.3e-07), theory 1.6e-08; best central h 3.2e-06 (error 2.1e-09), theory 1.3e-06; eps_F 2.7e-15
- [converged start] CRUISE_00 tol 1e-09: best forward h 5.6e-08 (error 5.9e-07), theory 6.6e-07; best central h 5.6e-07 (error 4.9e-09), theory 1.6e-05; eps_F 4.7e-12
- [converged start] CRUISE_00 tol 1e-06: best forward h 5.6e-07 (error 8.1e-06), theory 3.6e-05; best central h 5.6e-07 (error 2.3e-09), theory 2.3e-04; eps_F 1.4e-08
- [converged start] CRUISE_00 tol 0.0001: best forward h 5.6e-05 (error 8.1e-04), theory 3.6e-04; best central h 5.6e-05 (error 9.0e-07), theory 1.1e-03; eps_F 1.4e-06
- [converged start] CRUISE_07 tol 1e-13: best forward h 1.0e-08 (error 6.5e-07), theory 1.4e-08; best central h 3.2e-06 (error 2.7e-09), theory 1.5e-06; eps_F 2.6e-15
- [converged start] CRUISE_07 tol 1e-09: best forward h 1.0e-08 (error 6.5e-07), theory 8.5e-07; best central h 3.2e-06 (error 3.5e-09), theory 2.3e-05; eps_F 9.6e-12
- [converged start] CRUISE_07 tol 1e-06: best forward h 1.0e-06 (error 3.4e-05), theory 5.5e-05; best central h 3.2e-06 (error 3.5e-09), theory 3.7e-04; eps_F 4.0e-08
- [converged start] CRUISE_07 tol 0.0001: best forward h 1.0e-04 (error 3.5e-03), theory 5.5e-04; best central h 1.0e-04 (error 4.2e-06), theory 1.7e-03; eps_F 4.0e-06
- [converged start] SLS_09 tol 1e-13: best forward h 1.0e-07 (error 8.0e-09), theory 1.7e-08; best central h 1.0e-06 (error 3.7e-09), theory 1.7e-06; eps_F 1.7e-15
- [converged start] SLS_09 tol 1e-09: best forward h 3.2e-08 (error 3.3e-07), theory 2.6e-07; best central h 3.2e-07 (error 1.6e-09), theory 1.1e-05; eps_F 4.1e-13
- [converged start] SLS_09 tol 1e-06: best forward h 1.0e-06 (error 1.2e-05), theory 2.6e-05; best central h 1.0e-06 (error 5.5e-09), theory 2.3e-04; eps_F 4.0e-09
- [converged start] SLS_09 tol 0.0001: best forward h 1.0e-04 (error 1.2e-03), theory 2.6e-04; best central h 1.0e-04 (error 1.3e-05), theory 1.1e-03; eps_F 4.0e-07
- [offset start] SLS_00 tol 1e-13: best forward h 3.2e-08 (error 1.2e-07), theory 2.0e-09; best central h 1.0e-06 (error 2.7e-09), theory 3.8e-07; eps_F 1.3e-15
- [offset start] SLS_00 tol 1e-09: best forward h 3.2e-08 (error 2.1e-07), theory 7.1e-09; best central h 1.0e-06 (error 3.0e-09), theory 8.8e-07; eps_F 1.6e-14
- [offset start] SLS_00 tol 1e-06: best forward h 1.0e-09 (error 7.0e-07), theory 6.7e-07; best central h 5.6e-05 (error 7.0e-07), theory 1.8e-05; eps_F 1.4e-10
- [offset start] SLS_00 tol 0.0001: best forward h 1.0e-09 (error 7.0e-07), theory 6.7e-07; best central h 3.2e-05 (error 9.9e-07), theory 1.8e-05; eps_F 1.4e-10
- [offset start] CRUISE_00 tol 1e-13: best forward h 5.6e-08 (error 5.6e-07), theory 3.4e-08; best central h 3.2e-06 (error 3.6e-09), theory 2.2e-06; eps_F 1.2e-14
- [offset start] CRUISE_00 tol 1e-09: best forward h 3.2e-08 (error 1.0e-06), theory 4.5e-07; best central h 5.6e-06 (error 2.5e-07), theory 1.2e-05; eps_F 2.2e-12
- [offset start] CRUISE_00 tol 1e-06: best forward h 3.2e-08 (error 1.0e-06), theory 8.8e-07; best central h 5.6e-08 (error 4.9e-07), theory 1.9e-05; eps_F 8.5e-12
- [offset start] CRUISE_00 tol 0.0001: best forward h 5.6e-05 (error 5.3e-04), theory 2.0e-04; best central h 1.8e-03 (error 4.7e-04), theory 7.1e-04; eps_F 4.2e-07
- [offset start] CRUISE_07 tol 1e-13: best forward h 1.8e-08 (error 4.5e-07), theory 2.2e-08; best central h 3.2e-06 (error 2.6e-09), theory 2.0e-06; eps_F 6.2e-15
- [offset start] CRUISE_07 tol 1e-09: best forward h 3.2e-07 (error 3.2e-06), theory 1.2e-06; best central h 3.2e-05 (error 4.4e-07), theory 2.9e-05; eps_F 1.8e-11
- [offset start] CRUISE_07 tol 1e-06: best forward h 3.2e-07 (error 3.2e-06), theory 2.2e-06; best central h 3.2e-10 (error 3.2e-06), theory 4.3e-05; eps_F 6.4e-11
- [offset start] CRUISE_07 tol 0.0001: best forward h 3.2e-12 (error 1.7e-03), theory 2.0e-04; best central h 3.2e-12 (error 1.5e-03), theory 8.7e-04; eps_F 5.2e-07
- [offset start] SLS_09 tol 1e-13: best forward h 3.2e-08 (error 2.6e-08), theory 2.0e-08; best central h 1.0e-06 (error 3.1e-09), theory 2.0e-06; eps_F 2.4e-15
- [offset start] SLS_09 tol 1e-09: best forward h 1.0e-08 (error 3.7e-06), theory 2.5e-06; best central h 1.8e-05 (error 8.6e-07), theory 4.8e-05; eps_F 3.6e-11
- [offset start] SLS_09 tol 1e-06: best forward h 1.0e-08 (error 3.7e-06), theory 2.5e-06; best central h 3.2e-05 (error 7.9e-07), theory 4.8e-05; eps_F 3.6e-11
- [offset start] SLS_09 tol 0.0001: best forward h 3.2e-05 (error 1.0e-04), theory 1.5e-04; best central h 3.2e-11 (error 4.0e-04), theory 7.5e-04; eps_F 1.4e-07
- Caveat (measured, not assumed): with the warm start from the converged point, a perturbation whose initial residual is already below tol produces ZERO solver iterations, so F is evaluated with the flow-matching variables frozen; the resulting partial derivative is wrong by ~15-30 % (plateau at small h). The best step therefore sits just above that dead zone, not at the noise-model optimum. With the 0.1 % offset start the solver error is a smooth function of the state and largely cancels in the difference, so the classical noise/h rise does not appear.

- Log-log slope of the best forward step vs tol (tol 1e-9..1e-4; 0.5 would be the classical sqrt(noise) law): converged/CRUISE_00 0.58, converged/CRUISE_07 0.79, converged/SLS_00 0.63, converged/SLS_09 0.68, offset/CRUISE_00 0.60, offset/CRUISE_07 -0.92, offset/SLS_00 -0.32, offset/SLS_09 0.64

- **The 'about 1e-3' part of H3 fails at every point and for both solver starts**: at tol = 1e-6 the best forward step is 1.0e-09 to 1.0e-06, i.e. 1000x to 1000000x smaller than 1e-3. With the warm start the best step does move with tol (>= 10x from 1e-13 to 1e-6); with the offset start it does not move consistently.
- The 'U-shape' at loose tolerances (1e-6, 1e-4) is NOT the classical noise/h rise: its left side is the zero-iteration plateau with a cliff (see the left_shape column of tables/h_sweep_optima.csv). The classical gradual left rise appears only at tight tolerances, where floating-point round-off is the noise.

### H3 - best h shifts with solver tolerance: **not supported**

- Rule (applied by code): the hypothesis has two parts, (a) the best forward step moves with tol (>= 10x from tol=1e-13 to tol=1e-6 at every Layer-2 point) and (b) at tol=1e-6 it lies near 1e-3 (within a factor 3.5). supported = (a) and (b) hold for both inner-solver starts; not supported = (b) fails at EVERY point for both starts, or (a) fails everywhere; otherwise mixed
- Numbers: `{"converged/CRUISE_00": {"err_min_fwd_tol1e-6": 8.05031940316486e-06, "h_opt_fwd_tol1e-13": 5.6234132519034905e-08, "h_opt_fwd_tol1e-6": 5.62341325190349e-07, "h_theory_fwd_tol1e-6": 3.6233141399487954e-05, "shift_factor": 10.0}, "converged/CRUISE_07": {"err_min_fwd_tol1e-6": 3.44353023149175e-05, "h_opt_fwd_tol1e-13": 1e-08, "h_opt_fwd_tol1e-6": 1e-06, "h_theory_fwd_tol1e-6": 5.493739946849818e-05, "shift_factor": 100.0}, "converged/SLS_00": {"err_min_fwd_tol1e-6": 9.19689028947555e-06, "h_opt_fwd_tol1e-13": 3.162277660168379e-08, "h_opt_fwd_tol1e-6": 3.162277660168379e-07, "h_theory_fwd_tol1e-6": 4.4787297767273536e-06, "shift_factor": 10.0}, "converged/SLS_09": {"err_min_fwd_tol1e-6": 1.1691601411972596e-05, "h_opt_fwd_tol1e-13": 1e-07, "h_opt_fwd_tol1e-6": 1e-06, "h_theory_fwd_tol1e-6": 2.606989019234659e-05, "shift_factor": 10.0}, "offset/CRUISE_00": {"err_min_fwd_tol1e-6": 9.99197146143944e-07, "h_opt_fwd_tol1e-13": 5.6234132519034905e-08, "h_opt_fwd_tol1e-6": 3.162277660168379e-08, "h_theory_fwd_tol1e-6": 8.836252996885425e-07, "shift_factor": 0.5623413251903491}, "offset/CRUISE_07": {"err_min_fwd_tol1e-6": 3.239481783741398e-06, "h_opt_fwd_tol1e-13": 1.7782794100389228e-08, "h_opt_fwd_tol1e-6": 3.162277660168379e-07, "h_theory_fwd_tol1e-6": 2.1953911269996687e-06, "shift_factor": 17.78279410038923}, "offset/SLS_00": {"err_min_fwd_tol1e-6": 7.002371655362772e-07, "h_opt_fwd_tol1e-13": 3.162277660168379e-08, "h_opt_fwd_tol1e-6": 1e-09, "h_theory_fwd_tol1e-6": 6.658396658109643e-07, "shift_factor": 0.0316227766016838}, "offset/SLS_09": {"err_min_fwd_tol1e-6": 3.6948317250989355e-06, "h_opt_fwd_tol1e-13": 3.162277660168379e-08, "h_opt_fwd_tol1e-6": 1e-08, "h_theory_fwd_tol1e-6": 2.4612686695109745e-06, "shift_factor": 0.31622776601683794}}`

## Layer 3 - Bauer-Fike fragility

- idle (SLS) / max-cruise worst-shift ratio per seed: 1.407, 1.415, 1.405
- Spearman(kappa(V), tau_slow) over operating points: 0.131
- Monte Carlo: n = 20000 per point per seed, seeds ['20261114', '20261115', '20261116']; max relative seed spread of the worst shift: 7.30e-03

- Within each sweep (seed 20261114): CRUISE slowest CRUISE_00, most fragile CRUISE_01; SLS slowest SLS_00, most fragile SLS_01
- Contrast with the published turboshaft example (kappa(V) 18.9 at idle, ~13x idle/cruise ratio): on the study branch the AGTF30 port's idle eigenvalues are well separated, kappa(V) stays between 2.97 and 6.87, so the idle/cruise ratio is only ~1.4x (the other idle steady state, below, has kappa(V) 15.8).

- If the OTHER SLS_00 steady state were used: tau_slow 2.207 s, kappa(V) 15.83, worst shift 0.00797 (bound 0.01583); it would be the most fragile point overall: True; the slowest overall: False; so H2's rule would then read: not supported.

### H2 - slowest point is the most fragile: **not supported**

- Rule (applied by code): over all SLS and cruise operating points (normalized A_n, ||E||_2 = 1e-3): supported if the point with the largest tau_slow is also the point with the largest worst-case eigenvalue shift for every Monte Carlo seed; not supported if it is not for any seed; otherwise mixed
- Numbers: `{"20261114": {"fragility_rank_of_slowest": 12, "most_fragile_point": "SLS_01", "of_points": 18, "same": false, "slowest_point": "CRUISE_00", "worst_shift_max": 0.0035130114171337823, "worst_shift_slowest": 0.002194331667041194}, "20261115": {"fragility_rank_of_slowest": 12, "most_fragile_point": "SLS_01", "of_points": 18, "same": false, "slowest_point": "CRUISE_00", "worst_shift_max": 0.0035127751509782534, "worst_shift_slowest": 0.002194394863783755}, "20261116": {"fragility_rank_of_slowest": 12, "most_fragile_point": "SLS_01", "of_points": 18, "same": false, "slowest_point": "CRUISE_00", "worst_shift_max": 0.0035131049727250563, "worst_shift_slowest": 0.002194399942881564}}`

## O4 - degradation

- CRUISE_07/same_fan, largest |dtau_slow| first: lpc_flow +79.94 %, hpt_flow +69.22 %, fan_eff +62.49 %, lpt_eff +57.38 %, hpc_flow +54.27 %
- CRUISE_07/same_fuel, largest |dtau_slow| first: lpc_flow +87.02 %, lpt_eff +80.77 %, fan_eff +76.52 %, hpt_flow +70.25 %, hpc_flow +52.39 %
- SLS_00/same_fan, largest |dtau_slow| first: combined +11.80 %, lpt_flow +9.53 %, hpt_eff +9.25 %, fan_flow +7.51 %, lpt_eff -4.00 %
- SLS_00/same_fuel, largest |dtau_slow| first: combined +25.76 %, hpt_eff +15.21 %, lpt_flow +9.63 %, fan_flow +6.45 %, hpc_eff +5.28 %
- SLS_09/same_fan, largest |dtau_slow| first: lpt_eff -10.69 %, fan_flow +9.78 %, combined +9.17 %, fan_eff -7.96 %, lpc_flow -5.74 %
- SLS_09/same_fuel, largest |dtau_slow| first: combined +17.90 %, hpt_eff +7.58 %, lpc_flow -5.38 %, hpc_eff +5.21 %, fan_flow +4.22 %
- Caveat: at CRUISE_07 the healthy engine has (or the degraded engine acquires) a COMPLEX eigenvalue pair; there tau_slow = -1/Re(lambda) of the pair, and a degradation that splits the pair into two real eigenvalues makes tau_slow jump (large % changes are a real/complex transition, not a gradual slow-down). 22 of 66 cases affected.

### H4 - degradation makes the engine slower: **mixed**

- Rule (applied by code): combined degradation (all 10 health parameters reduced together) at the largest converged level: supported if tau_slow increases at every Layer-4 point under BOTH comparison rules; not supported if it decreases everywhere; otherwise mixed
- Numbers: `{"CRUISE_07/same_fan": {"dtau_fast_pct": -56.85537712953661, "dtau_slow_pct": 1.3149643975205016, "level": 0.05}, "CRUISE_07/same_fuel": {"dtau_fast_pct": -64.43867828745921, "dtau_slow_pct": -2.4322217978399285, "level": 0.05}, "SLS_00/same_fan": {"dtau_fast_pct": 5.841014696878388, "dtau_slow_pct": 11.804446320161954, "level": 0.05}, "SLS_00/same_fuel": {"dtau_fast_pct": -13.1620515285489, "dtau_slow_pct": 25.75525504059739, "level": 0.05}, "SLS_09/same_fan": {"dtau_fast_pct": -23.49814108051842, "dtau_slow_pct": 9.174881687331116, "level": 0.05}, "SLS_09/same_fuel": {"dtau_fast_pct": -20.317860764157093, "dtau_slow_pct": 17.903913087881172, "level": 0.05}}`

## O5 - new steady state after a 1 % fuel step (Gauss-Seidel)

- SLS_00: dx_ss = [9.597, 16.939] rpm (0.347 %, 0.098 % of speed); Gauss-Seidel 15 sweeps, rho 0.1038, matches direct to 5.8e-15; diagonally dominant (rpm form) True
- SLS_01: dx_ss = [13.556, 18.276] rpm (0.435 %, 0.103 % of speed); Gauss-Seidel 19 sweeps, rho 0.1701, matches direct to 6.4e-15; diagonally dominant (rpm form) True
- SLS_02: dx_ss = [16.153, 20.611] rpm (0.466 %, 0.113 % of speed); Gauss-Seidel 18 sweeps, rho 0.1557, matches direct to 8.8e-15; diagonally dominant (rpm form) True
- SLS_03: dx_ss = [19.060, 21.963] rpm (0.499 %, 0.117 % of speed); Gauss-Seidel 19 sweeps, rho 0.1652, matches direct to 4.2e-15; diagonally dominant (rpm form) True
- SLS_04: dx_ss = [20.297, 22.323] rpm (0.487 %, 0.117 % of speed); Gauss-Seidel 20 sweeps, rho 0.1755, matches direct to 2.4e-15; diagonally dominant (rpm form) True
- SLS_05: dx_ss = [20.396, 22.282] rpm (0.451 %, 0.114 % of speed); Gauss-Seidel 21 sweeps, rho 0.1908, matches direct to 1.9e-15; diagonally dominant (rpm form) True
- SLS_06: dx_ss = [21.359, 22.223] rpm (0.438 %, 0.111 % of speed); Gauss-Seidel 20 sweeps, rho 0.1886, matches direct to 9.0e-15; diagonally dominant (rpm form) True
- SLS_07: dx_ss = [22.161, 22.474] rpm (0.424 %, 0.111 % of speed); Gauss-Seidel 19 sweeps, rho 0.1613, matches direct to 2.8e-15; diagonally dominant (rpm form) True
- SLS_08: dx_ss = [23.112, 23.111] rpm (0.415 %, 0.112 % of speed); Gauss-Seidel 16 sweeps, rho 0.1089, matches direct to 1.8e-15; diagonally dominant (rpm form) True
- SLS_09: dx_ss = [23.643, 23.133] rpm (0.399 %, 0.110 % of speed); Gauss-Seidel 16 sweeps, rho 0.1006, matches direct to 4.5e-16; diagonally dominant (rpm form) True
- CRUISE_00: dx_ss = [12.548, 29.065] rpm (0.432 %, 0.184 % of speed); Gauss-Seidel 15 sweeps, rho 0.0947, matches direct to 1.3e-15; diagonally dominant (rpm form) True
- CRUISE_01: dx_ss = [15.229, 21.351] rpm (0.440 %, 0.128 % of speed); Gauss-Seidel 21 sweeps, rho 0.2031, matches direct to 6.6e-15; diagonally dominant (rpm form) True
- CRUISE_02: dx_ss = [16.715, 22.012] rpm (0.417 %, 0.126 % of speed); Gauss-Seidel 21 sweeps, rho 0.1939, matches direct to 2.4e-15; diagonally dominant (rpm form) True
- CRUISE_03: dx_ss = [17.761, 22.229] rpm (0.389 %, 0.122 % of speed); Gauss-Seidel 21 sweeps, rho 0.2068, matches direct to 9.9e-15; diagonally dominant (rpm form) True
- CRUISE_04: dx_ss = [18.259, 21.580] rpm (0.357 %, 0.114 % of speed); Gauss-Seidel 22 sweeps, rho 0.2211, matches direct to 8.4e-15; diagonally dominant (rpm form) True
- CRUISE_05: dx_ss = [18.882, 21.359] rpm (0.333 %, 0.110 % of speed); Gauss-Seidel 20 sweeps, rho 0.1877, matches direct to 7.2e-15; diagonally dominant (rpm form) True
- CRUISE_06: dx_ss = [19.343, 23.532] rpm (0.311 %, 0.117 % of speed); Gauss-Seidel 13 sweeps, rho 0.0658, matches direct to 2.1e-15; diagonally dominant (rpm form) True
- CRUISE_07: dx_ss = [18.548, 37.815] rpm (0.274 %, 0.181 % of speed); Gauss-Seidel 13 sweeps, rho 0.0625, matches direct to 4.7e-16; diagonally dominant (rpm form) True
- ANCHOR_CRUISE: dx_ss = [18.562, 37.652] rpm (0.274 %, 0.181 % of speed); Gauss-Seidel 13 sweeps, rho 0.0618, matches direct to 3.8e-16; diagonally dominant (rpm form) True
- ANCHOR_TAKEOFF: dx_ss = [24.221, 23.575] rpm (0.398 %, 0.109 % of speed); Gauss-Seidel 15 sweeps, rho 0.0970, matches direct to 3.1e-15; diagonally dominant (rpm form) True

## Linearization validity

- SLS_00: +1 % fuel -> 0.55 % error; +3 % fuel -> 1.49 % error; +10 % fuel -> 4.63 % error; +25 % fuel -> 9.85 % error; valid (<= 2 %) up to 3 %, (<= 5 %) up to 10 % steps
- CRUISE_07: +1 % fuel -> 1.11 % error; +3 % fuel -> 3.21 % error; +10 % fuel -> 9.42 % error; +25 % fuel -> 13.21 % error; valid (<= 2 %) up to 1 %, (<= 5 %) up to 3 % steps
- SLS_09: +1 % fuel -> 0.46 % error; +3 % fuel -> 1.33 % error; +10 % fuel -> 3.99 % error; +25 % fuel -> 10.50 % error; valid (<= 2 %) up to 3 %, (<= 5 %) up to 10 % steps

## Published turboshaft example (JetCat SPT5, illustration only, not the AGTF30)

- idle: lam -0.3461, -0.3739; tau 2.889 / 2.675; kappa(V) 18.875
- mid: lam -0.5992, -0.8908; tau 1.669 / 1.123; kappa(V) 2.393
- cruise: lam -1.0127, -1.7873; tau 0.987 / 0.560; kappa(V) 1.836
- Monte Carlo worst shift idle 0.0162, cruise 0.00119, ratio 13.58 (brief: ~0.016, ~0.0012, ~13x)
- tau values matched to the brief only to its printed rounding: idle tau_fast 2.6748 vs 2.67, mid tau_fast 1.1226 vs 1.12, cruise tau_slow 0.9875 vs 0.99

