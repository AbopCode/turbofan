# AUDIT REPORT: "The Eigenvalue Map of a Turbofan Engine" (turbofan_eigenmap)

Auditor work area: `/home/claude/audit` (checkout, `delivered_outputs/` = backup of the builders' outputs taken before the
clean rebuild, `out2/` = second run, `scripts/` = auditor's own code, `sheets/` = figure contact sheets, `run1.log`, `run2.log`,
`h1.txt`, `h2.txt`).

## Summary verdict

The software works: it rebuilds cleanly, all 521 checks pass, the two runs give identical hashes, and the published-example numbers
check out. My own recomputation of the table eigen-quantities matches, and so does an independent re-linearization at max cruise.
**The study is not yet safe to hand in, for two HIGH findings:**

1. **The SLS idle operating point cannot be reproduced from its printed inputs, and the engine has more than one steady state for
   the same controls.** The package's own `solve_steady` gives a different equilibrium with a different A. At the table's own
   (x0, u0), the inner solver from its default start finds a root where the LP shaft is far from equilibrium. None of the
   documents mention this.
2. **H3's "about 1e-3" part is contradicted at every point, but the rule and the wording make the result read as "mixed".**
   The JSON records `near_1e-3_somewhere = false` for both starts. No document says in words that the location claim fails.

There are also several MEDIUM problems with labels and evidence: a y-axis that is wrong by a factor of 100, "within 1%" check
labels that use wider tolerances, a "U-shape" check that passes because of a solver artifact, and the missing mode-shape output.

## Clean-room rebuild

| item | run 1 (`outputs/`) | run 2 (`EIGENMAP_OUT=/home/claude/audit/out2`) |
|---|---|---|
| command | `rm -rf .venv outputs; ./run.sh full` (uv at /root/.local/bin, CPython 3.13.18) | `EIGENMAP_OUT=... ./run.sh full` (reuses .venv) |
| wall time | **98.8 s** (includes venv creation and package install); stages 85.7 s | **86.0 s**; stages 82.0 s |
| exit code | 0 | 0 |
| stages | 9 ok, 0 failed checks, 0 crashed, 0 skipped | same |
| PASS/FAIL table | **521 passed, 0 failed** | 521 passed, 0 failed |
| unit tests | test_mathcore 68 PASS, test_anchors 38 PASS, test_engine 23 PASS | same |

The run stays well inside the "< 6 min on 2 CPUs" target (machine: 2 CPUs).
`EIGENMAP_OUT` works as documented: it is read in `eigenmap/config.py:14` and `verify.py:54`.

## Determinism

- `python -m eigenmap.cli hashes` for run 1 and run 2: **identical**. All 16 files (8 `data/*.json` + 8 `tables/*.csv`) match,
  COMBINED = `447a0300...6df699`.
- The builders' delivered `outputs/` also hash identically to my rebuild.
- All 18 PNGs and the .md/.tex documents are byte-identical between the runs.
- **Differs:** `tables/tables.xlsx` (run 1 vs run 2). It is outside the hash set; the cause is probably zip/xlsx timestamps. LOW.

## Independent recomputations

### (a) Three random rows of operating_points.csv (`numpy.random.default_rng(7).choice(20, 3)`): CRUISE_01, CRUISE_03, CRUISE_07

My code: closed-form 2x2 eigenvalues, `A_n = D^-1 A D`, unit-column eigenvectors, kappa = 2-norm condition number
(`scripts/recompute.py`).

| row | quantity | mine | table | rel diff |
|---|---|---|---|---|
| CRUISE_01 | lambda1 / lambda2 | -0.557028 / -1.54779 | -0.557 / -1.5478 | 5.0e-5 / 6.4e-6 |
| | tau_slow / tau_fast | 1.79524 / 0.646082 | 1.79524 / 0.646082 | 1.1e-6 / 7.2e-7 |
| | kappaV_norm / kappaV_rpm | 5.37627 / 1.16045 | 5.37625 / 1.16045 | 4.2e-6 / 2.1e-6 |
| | tr / det / disc | -2.10482 / 0.862162 / 0.981609 | -2.10482 / 0.862166 / 0.981611 | 9.5e-7 / 4.3e-6 / 1.5e-6 |
| CRUISE_03 | lambda1 / lambda2 | -0.81837 / -2.18701 | -0.8184 / -2.187 | 3.7e-5 / 4.8e-6 |
| | tau_slow / tau_fast | 1.22194 / 0.457245 | 1.22194 / 0.457245 | 1.4e-6 / 4.4e-7 |
| | kappaV_norm / kappaV_rpm | 4.33469 / 1.08949 | 4.3347 / 1.08949 | 2.0e-6 / 1.1e-6 |
| CRUISE_07 | lambda | -1.39668 ± 0.16266i | -1.3967 ± 0.1627i | 3.3e-5 |
| | tau_slow = tau_fast | 0.715986 | 0.715986 | 2.5e-7 |
| | kappaV_norm / kappaV_rpm | 4.59484 / 4.73014 | 4.59496 / 4.73028 | 2.6e-5 / 2.9e-5 |
| | disc | -0.105834 | -0.105828 | 5.3e-5 |

**Max rel diff = 5.3e-5.** This is fully explained by the CSV printing a_ij to 6 significant figures and lambdas to 4 decimals.
AGREE.

### (b) Independent re-linearization (`scripts/relin.py`, my own central differences, inner tol 1e-13)

**CRUISE_07** (35k ft, M0.8, Nc_fan 2360.x, VAFN/VBV held at the table values). `solve_steady` re-solve gives
x0 = [6773.536, 20854.749], Wf = 0.782109 (table: 6773.54, 20854.8, 0.782113). F(x0) = [1.4e-11, -1.1e-10] rpm/s.

| step (relative) | A (mine) | max elementwise rel diff vs table A | B rel diff | eigenvalues (mine) |
|---|---|---|---|---|
| 1e-3 | [[-1.69650, 0.18962], [-0.61447, -1.09730]] | 9.6e-4 (truncation) | 1.6e-5 | -1.39690 ± 0.16358i |
| 1e-4 | [[-1.69656, 0.18945], [-0.61432, -1.09683]] | **5.2e-5** | 3.0e-6 | -1.39670 ± 0.16269i |

The table A is [[-1.69655, 0.18944], [-0.614355, -1.0968]]. AGREE to printing precision.

**SLS_00 (idle): DOES NOT REPRODUCE.** See FAILURE 1.

- `solve_steady(target=("Nc_fan", 891.586), VAFN=8000, VBV=0.578703, tol=1e-13)` converges (resnorm 7.6e-14) to
  x0 = [2763.92, **16542.4**], Wf = **0.15397**, Fnet 4177, Ps3 102.8. The table has [2763.92, 17230.0], Wf 0.232879,
  Fnet 4242, Ps3 133.53.
- My A there is [[-2.436, 3.112], [0.239, -0.828]]; the table A is [[-2.745, 0.937], [0.467, -1.535]]. a12 differs by 232%, and the
  eigenvalues are -0.453 / -2.811 against the table's -1.244 / -3.036.
- `target=("Wf", 0.232879)` with the same VAFN/VBV gives a **second equilibrium with identical controls**: x0 = [4399.3, 17750.6],
  Fnet 14444, F/N ~ 1e-14.
- `target=("Ps3", 133.528)` gives a third point: x0 = [3856.0, 17315.9], Wf 0.196.
- At the **table's own x0 and u0**, `solve_flow` from its default start (and also from the two z0 above) converges to an inner root
  with W_in = 768.889 and Ps3 = 133.53, matching the table. But on that root **dN_LP/dt / N_LP = 1.397 1/s** and Fnet = 4177.3
  (table 4241.8). So the flow-matching system has more than one root at the same (x, u).
- From the default inner start, `solve_flow` at the tabulated (x0, u0) **fails** ("Singular matrix" / "math domain error") at 13 of
  the 20 tabulated points. It "converges" to a non-equilibrium root at SLS_00 and CRUISE_01.
- So F(x, u) is only well defined when the builders' warm start (z0, carried by continuation) is supplied.

## Published-example numbers (SPEC item 7), own numpy

| item | SPEC | mine | delivered |
|---|---|---|---|
| eig idle / mid / cruise | -0.346/-0.374, -0.599/-0.891, -1.013/-1.787 | -0.3461/-0.3739, -0.5992/-0.8908, -1.0127/-1.7873 | same |
| tau | 2.89/2.67, 1.67/1.12, 0.99/0.56 | 2.889/2.6748, 1.669/1.1226, 0.9875/0.5595 | same (docs note the 2.67/1.12/0.99 rounding) |
| kappa(V) | 18.9 / 2.39 / 1.84 | 18.8746 / 2.3929 / 1.8357 | same |
| GS cruise, b = [-1.2, -0.3] | 0.70588/0.65775, 0.74457/0.67886, 0.74582/0.67954, 0.74586/0.67956 | identical to 5 d.p.; direct 0.745856/0.679558 | same |
| rho cruise | ~0.032 | 0.032086 | 0.0320856 |
| a12 flip to -0.0018 | ~ -0.36 ± 0.008i | -0.36 ± 0.008246i | plotted, correct |
| MC worst shift idle / cruise, ratio | ~0.016 / ~0.0012, ~13x | 0.01621 / 0.001191, 13.6x (my seed 12345, Gaussian E normalized to ‖E‖2 = 1e-3) | 0.01617 / 0.001191, 13.58x |

ALL AGREE.

## Figures (all 18 PNGs viewed; all are 300 dpi with an SVG twin)

| file | findings |
|---|---|
| bauer_fike_scale_sweep.png | Axes and units OK. The three dotted Bauer-Fike bound lines differ only by colour, so they cannot be told apart in B&W (LOW). |
| bauer_fike_shift_vs_bound.png | OK: distinct markers, labelled axes with units. Consistent with the bauer_fike numbers (max worst shift 0.00351 at kappa\*‖E‖ ≈ 0.0069). |
| degradation_fuel_penalty_cruise.png | Shows only the 5 **efficiency** cases plus combined. The **flow-capacity** cases (in degradation.csv) are missing. Same-fan rule only (LOW). The combined 5.9% penalty is below lpt_eff alone (9.4%) and fan_eff alone (7.2%), and the figure does not explain this (LOW, needs a sentence). |
| degradation_tau_cruise.png | Consistent with the H4 numbers (combined +1.3% same-fan, -2.4% same-fuel). Individual curves reach +60 to +87% and are non-monotonic (fan_eff same-fan: 76% at 4%, then 62% at 5%), because CRUISE_07 is a merged complex pair (tau_slow = tau_fast) that splits under degradation. No caveat in the figure or in RESULTS_FACTS (LOW/MEDIUM interpretation risk, see FAILURE 8). Legend fonts are small. |
| degradation_tau_idle.png / degradation_tau_max.png | Consistent with the H4 numbers (SLS_00 +25.8/+11.8, SLS_09 +17.9/+9.2; lpt_eff -10.7). OK. |
| eigen_map_complex_plane.png | Correct values. The SLS slow mode (blue) and the cruise fast mode (orange) overlap on the real axis between -1.2 and -2.8, which is hard to read; the different markers keep it B&W-readable. Conjugates are shown as slow = +Im, fast = -Im, as the title explains. OK. |
| fd_error_vs_h_*_converged.png | At tol 1e-9/1e-6/1e-4 the left branch is a **flat plateau of ~0.1 to 0.3** (the zero-iteration "dead zone"), not an eps/h rise. Visual evidence for FAILURE 3. The legend does not explain the open "best h" markers (only the title does). |
| fd_error_vs_h_*_offset.png | At tol 1e-4 the curve is flat (~1e-3) down to h = 1e-12, with no U. The "best h" markers land at 1e-9 (SLS_00) and ~3e-12 (CRUISE_07), which are arbitrary points on a plateau (FAILURE 7). At SLS_00, tol 1e-6 and 1e-4 give identical curves, so the loose-tolerance case is not actually exercised there. |
| gauss_seidel_convergence.png | OK. The SLS_00 and SLS_09 lines overlap (rho 0.104 vs 0.101); the markers differ. |
| graphical_summary.png | Very small fonts (about 5 to 6 pt at print size). The H3 panel does not say which point, start or scheme it shows. In the H4 panel "cruise" means CRUISE_07 (max cruise), not cruise idle, which is ambiguous (LOW). |
| kappa_vs_power.png | Values agree with the table (SLS peak 6.87 at SLS_01; cruise min 2.97; CRUISE_07 rpm 4.73). OK. |
| **linear_vs_nonlinear_steps.png** | **The y-axis "δN/N0 per unit δWf/Wf0 [-]" is wrong by a factor of 100.** The plotted plateaus are ≈35 (SLS idle N_LP), ≈9.8 (N_HP) and ≈27 (cruise N_LP). The O5 data give 0.347% / 0.098% / 0.274% per 1% fuel, i.e. 0.347 / 0.098 / 0.274 per unit. **The four nonlinear curves are all solid lines that differ only by colour, so they cannot be told apart in B&W.** The legend appears in the first panel only. (MEDIUM) |
| port_vs_nasa_validation.png | The takeoff BPR bar (-1.48%) lies **outside** the plotted ±1% band. The legend box covers the +1% dashed line. Hatching makes it B&W-safe. |
| published_example_turboshaft.png | Correct values. Labelled "turboshaft JetCat SPT5 (published example, illustration)", with units "1/model time unit". OK. |
| tau_vs_power.png | Consistent with the table; the complex-pair merge at 100% is visible. OK. |

## Verdict checks (does the verdict follow from the stored rule and numbers?)

| H | delivered | auditor | reason |
|---|---|---|---|
| H1 slower at low power | supported | **AGREE** | SLS: Spearman -1.00, ratio 2.25. Cruise: Spearman -0.976, ratio 4.18. Both are below -0.7 with ratios above 1.5. Caveat: the SLS idle tau (0.804 s) comes from an operating point that is not uniquely defined (FAILURE 1). Cruise tau_slow is non-monotonic at the top (0.693 at CRUISE_06, then 0.716 at CRUISE_07). |
| H2 slowest = most fragile | not supported | **AGREE** | The slowest point is CRUISE_00 (tau 2.995); the most fragile is SLS_01 (worst shift 0.00351); the slowest ranks 12 of 18 for all 3 seeds. The verdict holds even in rpm states, where the most ill-conditioned point is CRUISE_07 (kappa_rpm 4.73), which is still not the slowest. The contrast with the turboshaft (13x vs 1.4x) is reported openly. |
| H3 best h moves to ~1e-3 at tol ~1e-6 | mixed | **AGREE with the letter of the rule, DISAGREE with how it is presented** | Shift ≥ 10x: holds at every point for the converged start (factors 10, 100, 10, 10; three of these sit exactly on the half-decade grid limit) and fails at 3 of 4 points for the offset start. Location within 3.5x of 1e-3: **fails at every point for both starts** (best forward h at tol 1e-6 is 3.2e-7 to 1e-6, about 3 orders of magnitude too small; even the theoretical optimum is 4.5e-6 to 5.5e-5). The rule's "not supported" branch ignores location, so a decisive refutation of "about 1e-3" becomes "mixed". The converged-start shift is itself caused by the dead-zone artifact, as the docs say. FAILURE 2. |
| H4 degradation slows the engine | mixed | **AGREE** | Combined at 5%: +1.31, -2.43, +11.8, +25.8, +9.17, +17.9%. One decrease (CRUISE_07 same-fuel), so mixed follows. The sign is not assumed. All 396 cases converged. |

On hiding: H2 is reported openly (CODE_GUIDE Q13: "We report the contradiction"). H3 is softened (FAILURE 2). For the NASA anchors,
CODE_GUIDE Q15 quotes only the favourable quantities (FAILURE 5). The O4 rankings at CRUISE_07 lack the complex-pair caveat
(FAILURE 8).

Claims in the documents that I checked against outputs and found backed:
- O5 GS vs direct 9.9e-15
- replica 8.45e-13 ("1e-12")
- "≤ 2% for a 1% fuel step" (max 1.11%)
- slope 0.6 to 0.8 (0.58 to 0.79)
- MC seed spread 7.3e-3
- MC never exceeded the Bauer-Fike bound

Claims that are not quite right:
- "kappa(V) stays between 3.0 and 6.9": the minimum is 2.97 (rounding, trivial).
- CODE_GUIDE Q9 "A is diagonally dominant": true only in rpm form. In A_n at SLS_00, |a12_n| = 5.84 > |a11| = 2.74.
- CODE_GUIDE Q10 explains the U by "solver leftover error divided by h", which contradicts Q11 and the figures.

## Objectives trace (objectives_trace.md)

| O | evidence exists? | supports? |
|---|---|---|
| O1 | operating_points.csv yes; "docs/equations.tex" is listed under "Evidence (outputs/)" but the file there is `outputs/equations.tex` | Yes, with the SLS idle reproducibility caveat (FAILURE 1). Wrong path (LOW). |
| O2 | figures exist | Partly. It promises "eigenvectors", but no eigenvector or mode-shape output exists in the tables or figures (FAILURE 6). |
| O3 | the fd figures and h_sweep_optima.csv exist | Error vs h, yes. "Optimal h": the best-h values on the flat plateaus are not meaningful (FAILURE 7). "✓ met" is shown next to "H3: mixed" with no word on the 1e-3 refutation. |
| O4 | the figures and degradation.csv exist (incl. dRe(lam_slow)/dtheta) | Yes. The fuel-penalty figure omits the flow cases (LOW). |
| O5 | the GS figure and table exist | Yes (max 9.9e-15, rho ≤ 0.221). |

The "✓ met" statuses come from PASS checks only, so they cannot reflect FAILURES 1, 2, 3 and 6.

## Honesty findings

- SOURCES.md lists 4 "Opened" sources. RUN_LOG.md records the arXiv:1403.3738 WebFetch. Unopened references are clearly
  separated, with an instruction to open them before citing. No problem found. I could not independently check that the NTRS PDF
  was opened; the Table 1 values used match the SPEC.
- AGTF30 as a production engine: never presented that way. "Not a production engine" appears in README, RESULTS_FACTS,
  ASSUMPTIONS E1 and equations.tex.
- Turboshaft example: labelled "illustration only / NOT the AGTF30" in every place I grepped (tables .tex captions, GS table rows,
  figure title, RESULTS_FACTS heading, CODE_GUIDE, README). Units are labelled "1/model time unit". OK.
- Softening: H3 (FAILURE 2), CODE_GUIDE Q15 anchor selection (FAILURE 5), and the "within 1%" check labels (FAILURE 5).
- Undisclosed: the multiple steady states / inner roots at SLS idle (FAILURE 1). A grep for unique/multiple/branch/second root finds
  nothing relevant.

## FAILURES

1. **HIGH: SLS idle steady state is not unique and not reproducible from the printed inputs.**
   - Files: `outputs/tables/operating_points.csv` row SLS_00, `outputs/data/layer1_opmap.json`.
   - Evidence: re-solving with the table's Nc_fan, VAFN and VBV converges (7.6e-14) to N_HP 16542 / Wf 0.154, against the table's
     17230 / 0.2329. A there differs (a12 +232%; lambda -0.453/-2.811 vs -1.244/-3.036).
   - With identical Wf/VAFN/VBV there is a second equilibrium at N_LP 4399.
   - At the table's own (x0, u0) the inner solver finds a root with dN_LP/dt/N = 1.397 1/s (Fnet 4177 vs 4242).
   - The default-start `solve_flow` fails at 13 of the 20 tabulated points.
   - Nothing in the documents mentions this. The SLS_00 numbers feed H1, H2 and H4.
   - Possible cause (unverified): the bypass-nozzle sub-ambient clamp (ASSUMPTIONS E9) at VAFN = 8000 in^2. The 7th flow variable
     is 3.0, 40.7 and 150.1 on the three roots.
2. **HIGH: H3 contradiction softened.**
   - Files: `outputs/data/layer2_fd.json` verdicts, `RESULTS_FACTS.md` H3 section, `CODE_GUIDE.md` Q10/Q11, `objectives_trace.md` O3.
   - Evidence: `near_1e-3_somewhere: false` for both starts. The best forward h at tol 1e-6 is 1e-9 to 1e-6, about 1000x from 1e-3.
   - The rule's "not supported" condition ignores location, so the verdict is "mixed". No document states that the 1e-3 claim is
     contradicted.
3. **MEDIUM: the "U-shaped error curve at loose tolerances" acceptance check passes on an artifact.**
   - Files: `verify_report.txt` lines "U-shaped forward-difference error curve" (ratios 97 to 37,490 against a threshold of 10),
     `fd_error_vs_h_*`.
   - Evidence: the left branch is the zero-iteration dead-zone plateau (15 to 30% error), not eps/h. With the offset start the
     curves are flat and show no rise (RESULTS_FACTS admits this).
   - CODE_GUIDE Q10 still teaches the classical U.
4. **MEDIUM: linear_vs_nonlinear_steps.png has a y-axis wrong by 100x and is not B&W-safe.**
   - Evidence: plateaus of 35 / 9.8 / 27 against 0.347 / 0.098 / 0.274 from the O5 dx_ss.
   - The four nonlinear step sizes are all solid lines distinguished only by colour.
5. **MEDIUM: NASA anchor checks are labelled "within 1%" but use widened thresholds.**
   - Evidence: the takeoff BPR threshold is 2.79% and the measured error is 1.48%. OPR takeoff (threshold 2.35%) and cruise BPR
     (3.08%) are also widened. The 1.48% bar sits outside the 1% band in `port_vs_nasa_validation.png`.
   - CODE_GUIDE Q15 quotes only the quantities within 0.32% / 0.67%. The SPEC says "port within 1% of NASA anchors". The
     rounding argument is reasonable but has to be stated in the check label and in the viva answer.
6. **MEDIUM: the Layer-1 "mode shapes vs power" output is missing.**
   - Files: `operating_points.csv` (no eigenvector columns), `figures/` (no mode-shape figure).
   - V_n exists only inside `layer1_opmap.json`, yet O2 promises eigenvectors.
7. **LOW: "best h" is reported where the error curve is flat.**
   - Examples: offset CRUISE_07 tol 1e-4 best h 3.2e-12 (forward and central), SLS_09 central 3.2e-11. Log-log slopes of -0.92 and
     -0.32 are reported without a warning (`RESULTS_FACTS.md`, `h_sweep_optima.csv`).
8. **LOW: the CRUISE_07 degradation rankings (+80% tau_slow) come from complex-pair splitting at a disc ≈ -0.106 point.**
   - They are listed in `RESULTS_FACTS.md` O4 without that caveat.
9. **LOW: `degradation_fuel_penalty_cruise.png` omits the flow-capacity cases and the same-fuel rule.**
   - The combined penalty being below single-component penalties is not explained.
10. **LOW: `objectives_trace.md` O1 cites `docs/equations.tex` under "Evidence (outputs/)".** The outputs file is `outputs/equations.tex`.
11. **LOW: `tables/tables.xlsx` is not byte-reproducible between runs.** It is outside the hash set.
12. **LOW: small text and B&W issues.**
    - Typo "~15-30 %%" in RESULTS_FACTS.
    - graphical_summary.png has tiny fonts and an unlabelled H3 panel.
    - bauer_fike_scale_sweep.png bound lines differ only by colour.
    - CODE_GUIDE Q9 calls A diagonally dominant, which holds only in rpm form.
    - "kappa between 3.0 and 6.9" while the minimum is 2.97.

## Re-check after fix cycle 1

Checked: `/home/claude/audit/checkout2/turbofan_eigenmap` (outputs only; no rebuild). Python from the first checkout's venv
(`requirements.txt` identical). Auditor scripts are in `/home/claude/audit/scripts2/`; figure sheets are
`sheets/fix1_sheet0.png` and `sheets/fix1_sheet1.png`.

| # | status | evidence |
|---|---|---|
| 1 | **FIXED** | RESULTS_FACTS "Multiple steady states at idle" reports the other SLS_00 branch (Wf 0.1540, N_HP 16542, lam -0.4531/-2.8112, tau_slow 2.207, kappa 15.83). This matches my earlier re-solve exactly. `idle_branches` is in layer1_opmap.json and z0_* columns are in operating_points.csv. **Reproduction:** with the printed x0, u0 and z0 as warm start, F/N ≤ 1.7e-5 at every one of the 20 rows; this residual is consistent with x0 printed to 6 s.f. My central-difference A agrees within 7.6e-5 at SLS_00 and within 3.6e-4 at worst over all rows. **H1 sensitivity:** 2.207/0.3579 = 6.17, which matches. **H2 sensitivity:** my own MC on A_n of the other branch gives worst shift 0.00797 (delivered 0.00797, bound 0.01583). That exceeds the global max 0.00351, so it would be the most fragile point. Its tau_slow 2.207 < 2.995 (CRUISE_00), so it would not be the slowest, and H2 stays "not supported". Both statements follow from the numbers. Residual LOW: the spurious inner root at the study x0 (LPT_PR ≈ 150, dN_LP/dt/N = 1.40 1/s from the default z start) is covered only implicitly by "reproducible from its tabulated z0". |
| 2 | **FIXED** | The rule now has parts (a) and (b), with "not supported = (b) fails at every point for both starts". The verdict is **not supported**, and RESULTS_FACTS l.91 and README l.161 state that the 1e-3 part fails everywhere (about 1000x to 1e6x smaller). objectives_trace O3 shows "H3: not supported". |
| 3 | **FIXED** | The check label now reads "has an interior minimum, rising >=10x on both sides (left side shape: ...)". h_sweep_optima.csv has `fwd_left_branch`, and CODE_GUIDE Q10 now says the classical rise appears only at tight tol. Residual LOW: the jump heuristic (any >10x step) labels the tol=1e-13 warm-start curves of SLS_00 and SLS_09 "dead-zone plateau + cliff", although their left branch is a gradual round-off rise in the figure. |
| 4 | **FIXED** | Plateaus are now 0.35 / 0.098 (SLS idle) and 0.27 (cruise N_LP), matching the O5 dx_ss. Curves now use distinct markers and line styles, and the linear model is a black dotted line, so the figure is B&W-safe. The legend in the first panel slightly overlaps the curves (cosmetic). |
| 5 | **FIXED** | The check labels now read "\|rel err\| <= 1% + half a unit of Table 1's last printed digit". CODE_GUIDE Q15 now states BPR/OPR up to 1.48% (takeoff BPR 27.6 vs 28). Cosmetic: the legend in port_vs_nasa_validation.png still covers the +1% line. |
| 6 | **FIXED** | `tables/mode_shapes.csv` (+ .tex) and `figures/mode_shapes_vs_power.png` now exist, and objectives_trace O2 cites them. I spot-checked SLS_00 mode 1: atan(0.2488/0.9686) = 14.4°, which matches the plot. |
| 7 | **FIXED** (one questionable cell) | h_sweep_optima.csv now has `fwd/cen_min_well_defined`, flagging flat-plateau minima (e.g. offset CRUISE_07 tol 1e-4: "NO (flat)"). The one questionable cell: offset CRUISE_07 tol 1e-6 central best h 3.2e-10 is marked "yes", but that curve is flat in the figure. |
| 8 | **FIXED** | RESULTS_FACTS l.123 now carries the caveat that the CRUISE_07 % changes are a real/complex transition. |
| 9 | **FIXED** | `degradation_performance_{cruise,idle,max}.png` now show efficiency and flow-capacity cases under both rules (extra fuel and thrust change). fan_flow (-13% fuel) explains why combined < individual. |
| 10 | **FIXED** | objectives_trace O1 now cites `equations.tex`. |
| 11 | **PARTLY (not re-run)** | `evidence/build.py:130` now rewrites the xlsx zip entries with a fixed date. I did not do a second run to confirm byte identity. Current `cli hashes` COMBINED = 9406bd93... (17 files, now including mode_shapes.csv). RUN_LOG.md records no combined hash to compare against. |
| 12 | **PARTLY** | Fixed: the "%%" typo is gone; CODE_GUIDE Q9 now says diagonal dominance holds in rpm form only; the Bauer-Fike bounds have open markers (B&W-safe); the graphical-summary H3 panel is labelled. Still open: graphical_summary.png fonts remain about 5 pt; "kappa(V) stays between 3.0 and 6.9" (min 2.97) remains (l.106). That sentence ("idle eigenvalues are well separated ... ratio only ~1.4x") is now also branch-dependent, because the other idle branch has kappa 15.8, and it should say "on the study branch" (LOW). |

**New HIGH findings: none.**
