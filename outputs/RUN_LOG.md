# RUN_LOG.md

All times are sandbox wall-clock notes written by the orchestrator. Sandbox: Linux, 2 CPUs, Python 3.13.18, no MATLAB/Octave.

## Wave 0 (orchestrator) — 2026-11-14
- Working folder and uploads checked: no `agtf30/`, `T-MATS/` or `nasa_truth/*.mat` supplied. **No NASA A/B matrices were available**; validation uses Table 1 of NASA/TM-2017-219569 and the repo's own steady-state initial conditions.
- `git clone` of github.com/nasa/AGTF30 (commit a1521a7, 2020-06-01) and github.com/nasa/T-MATS (commit ad6e4d5, 2025-01-15) succeeded; vendored under `vendor/`.
- Opened with WebFetch: NTRS 20170011210 (Chapman & Litt, NASA/TM-2017-219569): Table 1 confirmed — cruise 35,000 ft/M0.8: 813 pps, TSFC 0.4637, BPR 24, OPR 55, TIT 3150 R, 6073 lbf; max thrust hot-day SLS: 1724 pps, TSFC 0.1751, BPR 28, OPR 37, TIT 3170 R, 28,622 lbf. The paper names the LP shaft N1, HP shaft N2, fan shaft Nf, gear 3.1:1 (the repo code names them N2, N3, N1).
- Opened with WebFetch: arXiv:1403.3738 (Pakmehr & Yucelen), Sec. 3.3 eqs (57)-(59): the three turboshaft A (and B) matrices confirmed as in the brief.
- CONTRACT.md and eigenmap/config.py written. Wave 1 launched: A (port), A2 (blind replica, isolated sandbox), M (math core + published-example anchors), D0 (docs/runners).

## Interruption — 12:49 IST
- Sub-agents A (port) and A2 (blind replica) were stopped by the account's 5-hour usage limit (resets 17:00 IST) about 12 minutes into Wave 1. A had written no files; A2 had only scratch dumps of the .mdl/.slx wiring in /home/claude/a2_work/scratch (its context is saved and it has NOT seen any port code).
- On disk at that moment: M -> eigenmap/mathcore.py, tests/test_mathcore.py (anchors.py not yet); D0 -> run.sh, run.bat, requirements.txt, verify.py, eigenmap/cli.py, docs/equations.tex (README not yet). None of these has been gate-checked by the orchestrator yet.
- Orchestrator decision: per the brief's fallback ("if you have no sub-agent tool, run the same waves yourself sequentially"), the orchestrator takes over the Port Engineer role (A) itself while agents are unavailable. The blind replica stays with agent A2 (to be resumed after the reset) so the Gate-1 independence is preserved. A wake-up is scheduled for 17:06 IST.

## Resumed 2026-11-17 (after the usage-limit pause)
- Agent A2 (blind replica) was stopped by the user; its code on disk (engine + steady solver, written blind) was used as-is for Gate 1. No new agents were started for Wave 2: the orchestrator built Layers 1-4, validity and the evidence pack itself (the brief's sequential fallback), with the same file ownership, contract and PASS/FAIL discipline.
- GATE 1: PASS, 140/140 checks, worst port-vs-replica disagreement 8.5e-13. One real porting bug in the port (nozzle entropy when Pt <= Pamb) was found by the replica and fixed; see handoff/gate1.md. Engine frozen.
- Port vs NASA Table 1 (faithful mode): airflow, TSFC, thrust within 0.32 %; TIT within 0.67 %; BPR/OPR within printed rounding (worst 1.48 %, take-off BPR printed "28"). OPR compared as product of compressor PRs (ASSUMED definition).

## Run 1 — SKELETON (`./run.sh full --fast`) — 2026-11-17
- Wall time 70 s on 2 CPUs. 353 checks passed, 0 failed.
- Fast verdicts (coarse grids, NOT for the report): H1 supported, H2 not supported, H3 mixed, H4 mixed.
- Fixes during the run: idle nonlinear simulation failed on solve_ivp's first trial step (inner-solver Jacobian singular far from the point) -> integrator step bounded by tau_fast (validity.py). Layer 2: found the zero-iteration "dead zone" of the warm-started solver; added the 0.1 % offset-start experiment; the U-shape acceptance check applies to the warm-start configuration only (the offset-start ratio is reported as information).
- Delivered v1.zip.

## Run 2 — FIDELITY (`./run.sh full`) — 2026-11-17
- Wall time 87 s on 2 CPUs (stages: layer2_fd 50 s, layer4_degradation 47 s, validity 14 s, evidence 12 s; two stages in parallel). Confirmed full settings in the JSON: 10 SLS + 8 cruise points + 2 anchors, 45-step h grid, 4 tolerances, 2 solver starts, 20,000 Monte Carlo samples x 3 seeds per point, 0-5 % degradation in 6 levels x 11 cases x 2 rules x 3 points, fuel steps 1/3/10/25 %.
- 521 checks passed, 0 failed (pipeline + test modules).
- Verdicts: H1 supported; H2 not supported (slowest point = cruise idle, fragility rank 12 of 18; most fragile = SLS part power 1, next to SLS idle; same for all 3 seeds); H3 mixed (warm start: best h moves 10x at tol 1e-6 and 500-10,000x at 1e-4, but sits at 3e-7..1e-6, not 1e-3; offset start: no consistent shift); H4 mixed (sea level: slower under both rules; max cruise: -2.4 % same fuel, +1.3 % same fan).
- Added after the run (re-ran layer2_fd, layer3_bauerfike, evidence only): h_opt-vs-tol log-log slope, within-sweep fragility ranking; hand-written docs copied into outputs/. Delivered v2.zip.

## GATE 2 — independent audit (clean-room) — 2026-11-17
- Auditor = separate agent given only the client spec and a clean copy of v2; never saw the builders' reasoning.
- Clean-room `./run.sh full` (fresh venv via uv): exit 0, 98.8 s incl. venv; 521 passed / 0 failed. Second run: hashes identical (16 files); only tables.xlsx differed.
- Independent recomputation of 3 seeded table rows: agree to the CSV print precision (<= 5.3e-5); independent re-linearization of CRUISE_07: A within 5.2e-5, B within 3e-6. Published turboshaft numbers all reproduced.
- 12 findings (2 HIGH, 4 MEDIUM, 6 LOW). The two HIGH ones:
  1. SLS idle has MORE THAN ONE steady state for the same Nc_fan, VAFN, VBV (W_f 0.233 vs 0.154; different A); undocumented, table not reproducible from its printed inputs. Orchestrator confirmed it is a genuine property of NASA's model (present in the faithful bilinear port too), not an artifact of the smoothed maps.
  2. H3 read "mixed" although its "about 1e-3" part fails at every point.
- Fix cycle 1 (orchestrator): both idle branches computed and reported with their eigenvalues; H1 and H2 re-evaluated on the other branch (H1 still supported, ratio 6.2; H2 still not supported); z0 added to operating_points.csv so every point is reproducible; H3 rule now says "not supported" when the location claim fails everywhere (verdict changed mixed -> not supported) and RESULTS_FACTS states it in words; U-shape check renamed truthfully with a left-branch shape classification; validity figure rescaled (was 100x) and made B&W-safe; NASA check labels state the widened tolerance; mode-shape table + figure added; flat-minimum flag; complex-pair caveat for degradation; complete degradation performance figures; xlsx made byte-reproducible; CODE_GUIDE Q9/Q10/Q13/Q15 corrected. Full report: handoff/AUDIT_REPORT.md.

## Run 3 — ADVERSARIAL — 2026-11-17
- After fix cycle 1: `./run.sh full` 518 pipeline checks + 129 test checks pass (one stage crashed on a format-string slip in new text, caught by the pipeline, fixed, stage re-run).
- Determinism: a second complete run into a separate folder gives identical sha256 for all 17 JSON/CSV files (combined 9406bd93...), identical RESULTS_FACTS.md, objectives_trace.md and (after the zip-timestamp fix) tables.xlsx.
- Monte Carlo spread over 3 seeds: worst-shift relative spread <= 0.73 % at every point; H2 verdict identical for all 3 seeds.
- Verdicts vs Run 2: H1 supported (same); H2 not supported (same); H3 mixed -> NOT SUPPORTED (rule corrected after audit; numbers unchanged); H4 mixed (same).

## Run 4 — 3D twin (extension) — 2026-11-17
- twin/: Dash app (offline, serve_locally) with schematic 3D cutaway (surfaces of revolution, rotors spinning in proportion to shaft speeds, sections coloured by model total temperature), gauges with NASA limits, throttle lever, scripted scenarios, NASA baseline min/max-select PI fuel control (gains/limits from setup_Controller.m), live s-plane pole map with kappa(V) and tau, linear-vs-nonlinear overlay, degradation slider, optional closed-loop poles (labelled extension).
- outputs/engine_viewer.html exported (plotly.js embedded, works offline). App launched headlessly, served HTTP 200; screenshots outputs/figures/twin_control_centre.png and twin_viewer.png taken and inspected. Found and fixed during inspection: a format-string bug that blanked two panels, an over-precise axis.
- Idle-to-max burst at SLS under the baseline controller: Ps3-min -> Wf/Ps3 acceleration schedule -> fan-speed loop; ~95 % thrust after ~5 s.
- Auditor re-check after fix cycle 1: 10/12 FIXED, 2 PARTLY, no new HIGH. Fix cycle 2 (LOW items): left-branch classifier now requires a flat plateau before the cliff (round-off curves no longer mislabelled); spurious inner flow-matching root documented; "on the study branch" wording; summary-figure fonts enlarged; legends moved off the data. xlsx byte-identity confirmed by the orchestrator's second run (identical md5).

## GATE 3 — integration test — 2026-11-17
- `rm -rf .venv outputs && ./run.sh full`: fresh venv (uv, CPython 3.13.18), 10/10 stages OK (incl. twin export + headless screenshots), 521 checks passed, 0 failed, wall time 2 min 49 s including package install (study stages ~95 s, twin export ~75 s).
- Verdicts after all fix cycles: H1 supported, H2 not supported, H3 not supported, H4 mixed. objectives_trace: O1-O5 all met.
- `./run.sh verify` (quick): 4 s, exit 0. Delivered final.zip.
