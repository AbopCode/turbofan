# ASSUMPTIONS.md

Every modelling choice of the study, tagged:

- **PUBLIC**: taken directly from a NASA source file or publication (traceable, not our choice).
- **DERIVED**: our own choice, made by a stated rule from public data.
- **ASSUMED**: our own choice where the sources are silent or ambiguous. Change it and results may change.

## Engine model (the port)

| # | Item | Tag | Detail |
|---|---|---|---|
| E1 | Engine | PUBLIC | NASA AGTF30, a conceptual N+3 30,000 lbf geared turbofan (github.com/nasa/AGTF30, commit a1521a7). Not a production engine. |
| E2 | Component maps, scalars, bleeds, duct losses, inertias, gear ratio 3.1, controller schedules | PUBLIC | Parsed programmatically from `AGTF30/SimSetup/setup_*.m` (no hand transcription); see `eigenmap/agtf30_data.py`. NASA's code comment notes the NPSS data gives gear ratio 3.0; the repo value 3.1 is used. |
| E3 | Component equations, gas properties, interpolation, residual definitions | PUBLIC | Ported from the T-MATS C sources (`TMATS_Library/MEX/*.c`, commit ad6e4d5) and the library block diagrams. Wiring from `AGTF30_eng.mdl`. |
| E4 | Shaft dynamics | PUBLIC | I (pi/30) dN/dt = net torque per shaft (T-MATS Shaft block). Fan torque reaches the LP shaft through the gearbox (efficiency `MWS.GearBox.Eff`); I_LP includes the fan inertia reflected through the gear. HP shaft carries NASA's constant 350 hp power off-take. |
| E5 | Two model modes | DERIVED | `mode="nasa"`: faithful port (T-MATS inner loops with NASA's thresholds; bilinear map interpolation). This mode is validated against NASA and against the blind replica. `mode="study"` (used for all eigenvalue results): the inner loops (gas-property inversions, duct static pressure, nozzle Mach-1 search) are converged to machine precision, and the maps use C1 bicubic Hermite interpolation through the same map points. Reason: the faithful model has derivative jumps at every map breakpoint (NASA's design point sits on one) and about 1e-6 internal loop noise, so A = dF/dx would not be well defined. Effect on the anchors: at most 0.27 % (checked in `validation.json`). |
| E6 | OPR definition for the Table 1 comparison | ASSUMED | Product of the fan, LPC and HPC pressure ratios. With Pt3/Pt2 (duct losses included) the port reads 2.3 % (cruise) and 2.7 % (take-off) low. Both values are reported. |
| E7 | TIT in Table 1 | ASSUMED | Compared with the burner-exit total temperature Tt4. |
| E8 | No volume (gas-dynamic) states | PUBLIC | T-MATS AGTF30 solves flow matching algebraically at every instant, so the only dynamic states are the two shaft speeds. |
| E9 | Faithful nozzle sub-ambient branch | PUBLIC | When the nozzle inlet pressure is <= ambient, Pt is replaced by Pamb + 0.1 psia and the entropy is recomputed from it (`PcalcStat_TMATS.c`). This was mis-ported at first and caught at Gate 1. |

## Linearization

| # | Item | Tag | Detail |
|---|---|---|---|
| L1 | States, input | PUBLIC/DERIVED | x = [N_LP, N_HP] (rpm), u = W_f (lbm/s). N_LP is the LPC/LPT shaft (AGTF30 code "N2", paper "N1"); N_fan = N_LP / 3.1. N_HP is the code's "N3" (paper "N2"). |
| L2 | Held constant in A and B | DERIVED | VAFN and VBV stay at their steady scheduled values; ambient conditions and theta are fixed. No actuator, sensor or controller dynamics. The full system has no other states (E8), so nothing else had to be held. |
| L3 | Reference Jacobian | DERIVED | Central difference + Richardson extrapolation (2 levels), relative steps 1e-5 and 5e-6, inner tol 1e-13; self-consistency between the two steps < 1e-6 relative at every point (checked). Complex-step is not possible for this engine (table look-ups, clamps, branches are not complex-analytic). |
| L4 | Normalized form | DERIVED | A_n = D^-1 A D, D = diag(x0) (states as fractions of the operating speed). The eigenvalues are identical to A's; kappa(V) and the Bauer-Fike study use A_n so that ||E||_2 = 1e-3 means the same relative error at every point. kappa(V) of the rpm form is reported as well. |
| L5 | Inner-solver start for Layer 2 | DERIVED | Primary: warm start from the converged z0 (what a linearization program does). Extra experiment: start 0.1 % away (alternating signs). |

## Operating points

| # | Item | Tag | Detail |
|---|---|---|---|
| P1 | Power setting | DERIVED | Set by corrected fan speed Nc_fan, NASA's thrust surrogate; VAFN/VBV on NASA's baseline schedules. |
| P2 | Idle | DERIVED | The point where Ps3 equals NASA's minimum-Ps3 limiter: 120 psia at sea level, 50 psia at 35,000 ft (linear in altitude, setup_Controller.m). |
| P3 | Max power | PUBLIC | SLS: Nc_fan of NASA's max-thrust point (0.81 s_Nc). Cruise sweep: Nc_fan of NASA's cruise point (1.0 s_Nc). |
| P4 | Points between | DERIVED | Equally spaced in Nc_fan: 10 at SLS standard day, 8 at 35,000 ft / M0.8 (full run). |
| P5 | Anchors | PUBLIC | Cruise 35,000 ft M0.8 with VAFN 4775.22 in^2, VBV 0; max thrust SLS +27 F with VAFN 6314.82, VBV 0 (setup_InputsSS.m). |

## Study layers

| # | Item | Tag | Detail |
|---|---|---|---|
| S1 | h | DERIVED | Relative step: variable j moves by h * abs(v0_j). Grid 1e-12 .. 1e-1, 4 per decade (full run). |
| S2 | Eigenvalue error | DERIVED | max over the two eigenvalues of abs(lam_FD - lam_ref) / abs(lam_ref), nearest-neighbour matched. |
| S3 | Noise level eps_F | DERIVED | Largest abs(F_tol - F_1e-14) at 6 nearby states, divided by the size of F for a 100 % change (abs(A) x0). M2, M3 from 2nd/3rd differences at h = 1e-2. The theoretical h* formulas assume random noise, which the measured solver error is NOT (see RESULTS_FACTS). |
| S4 | Monte Carlo E | DERIVED | Standard-normal entries rescaled to ||E||_2 = 1e-3 exactly; 20,000 per point per seed; seeds 20261114/15/16. |
| S5 | Health parameters | DERIVED | Multiplicative scalars on map efficiency and map corrected flow after NASA's map scalars; "x % degradation" = theta = 1 - x/100; "combined" = all ten reduced together. |
| S6 | Degradation rules | DERIVED | same_fuel: W_f fixed at the healthy value. same_fan: Nc_fan fixed (W_f re-trimmed); the fuel penalty is the W_f increase. |
| S7 | Linear-model validity | DERIVED | Fuel steps +1/3/10/25 %, VAFN/VBV held, simulated for 5 tau_slow; error = max_t abs(dx_nl - dx_lin) / max_t abs(dx_nl), worst shaft. |
| S8 | Hypothesis verdicts | DERIVED | Each verdict is computed by code from an explicit rule stored next to it in the JSON. |
| S9 | Published turboshaft example | PUBLIC | Pakmehr and Yucelen, arXiv:1403.3738, eqs. (57)-(59): JetCat SPT5 turboshaft, illustration only; states unlabeled, time unit unstated, 1-2 significant figures. |

## Virtual lab (extension, `twin/`)

| # | Item | Tag | Detail |
|---|---|---|---|
| L1 | "Max power" away from the anchors | ASSUMED | Corrected fan-speed demand at 100 % lever = 2359.98 rpm x interp(altitude, [0, 35,000 ft], [0.81, 1.00]); the controller limiters may stop the engine below it. "Idle" = the Ps3-min limit. Power % is linear in Nc_fan between the two. |
| L2 | Controller | DERIVED | NASA baseline structure (setup_Controller.m): gain-scheduled PI on Nc_fan, N1max, N3max, Ps3max, T45max, Ps3min and Wf/Ps3 ratio-unit limits, min/max select, written in velocity (incremental) form so switching loops is bumpless. Not a block-by-block port of the Simulink diagram. Lab switches: Kp x, Ki x, each limiter on/off. |
| L3 | Live integration | DERIVED | Explicit Euler, dt = 0.02 s, flow matching solved to 1e-9 every step; advanced by wall time x speed. Flight condition changes are slewed (2,500 ft/s, Mach 0.05/s, 5 F/s) and wear 2 %/s so the solver stays on its branch. |
| L4 | Actuators | ASSUMED | Rate limits: fuel valve 1 lbm/s per s (manual mode), VAFN 3,000 in^2/s, VBV 1 per s. No sensor dynamics. |
| L5 | Valid range | DERIVED | No start-up or sub-idle model: the live run pauses if N_LP < 25 % or N_HP < 55 % of design, or if flow matching fails (outside the component maps). |
| L6 | Engine drawing | ASSUMED | Schematic proportions; colours from the model's station total temperatures; blade rotation shown ~30x slower than real; particle density proportional to the model's mass flows. |

