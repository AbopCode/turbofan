# GATE 1 — port (A) vs blind replica (A2)

Run: `python -m tests.gate1_compare` (port in its faithful `mode="nasa"`: T-MATS inner loops with NASA's thresholds, bilinear maps).
Result: **PASS, 140/140 checks; worst relative disagreement 8.5e-13 (limit 1e-6).** Full output: `handoff/gate1_output.txt`.

Compared: raw flow-matching residuals, shaft torques, dN/dt, thrust, T4, Ps3, VBV flow, TSFC, OPR, T45, stall margins at
5 random states (+-3 % speeds, +-2 % flow-matching variables, +-5 % fuel, +-3 % VAFN, random VBV 0-0.3, one with 5 health
parameters degraded); and 6 steady points from each port's OWN Newton solver (cruise/takeoff with fixed Wf, cruise with
Nc_fan target, SLS part power and low power on NASA's VAFN/VBV schedules, SLS with hpc_eff = 0.97).

## Disagreement found and resolved
| # | Where | Cause (traced to T-MATS source) | Who was wrong | Fix |
|---|---|---|---|---|
| 1 | random state 3 (take-off, VBV 0.012): core-nozzle residual 1.0 (A) vs 0.775 (A2); Fnet differs 1.3 % | When the nozzle inlet total pressure is <= ambient, `Nozzle_TMATS_body.c` substitutes `Ptin = Pamb + 0.1`, and `PcalcStat_TMATS.c` RECOMPUTES the entropy `S = pt2sc(Pt, ...)` from the pressure it is given. The port re-used the entropy of the original (sub-ambient) pressure, which made the static temperature exceed Tt, clamp, and give zero velocity. | Port (A) | `engine_model._nozzle` now computes S from the substituted `Ptin`. Only reachable at off-design states inside solver iterations (never at a converged operating point), but it is a genuine porting error. |

## Independence notes (honest)
- A2 worked in a separate sandbox and never saw the port; its package imports nothing from `eigenmap` (checked by grep) and has its own MATLAB parser.
- Residuals agree to the last bit at most states because both ports transcribe the same C code operation by operation; the steady points (different Newton solvers) agree to 1e-13..1e-14.
- A2 was stopped (account usage limit, then stopped by the user) before writing its own anchor test and report. Its code was used exactly as it left it; `replica/` is that code.

## Engine frozen
`eigenmap/engine_model.py`, `steady_state.py`, `tmats_thermo.py`, `agtf30_data.py` frozen after this gate (2026-11-17).
