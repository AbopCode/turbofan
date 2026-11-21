# The Eigenvalue Map of a Turbofan Engine

A study for the course *Matrix Algebra and Calculus*: how the two eigenvalues of a
linearized turbofan-engine model move across the operating envelope, how accurately
they can be computed, and how sensitive they are.

## What this project is (and is not)

| Part | What it is | Status word |
|---|---|---|
| **Port** | A Python port of NASA's open-source **AGTF30** geared-turbofan simulation and of the T-MATS component library it is built on (`eigenmap/engine_model.py`, maps read from `vendor/`). | PORTED |
| **Study** | Linearization of the port about steady operating points, `A = dF/dx`, `B = dF/du`, and a four-layer eigenvalue study (`eigenmap/layer*.py`, `eigenmap/steps.py`, `eigenmap/validity.py`). | DERIVED |
| **Virtual lab** | A local browser app (`twin/`, Dash, fully offline): a Simulator tab with 16 one-click experiments, a Live engine tab (animated cutaway, Start/Stop, every input adjustable while it runs) and seven experiment tabs that expose every parameter. It is an *extension* for teaching; the study's numbers do not depend on it. | EXTENSION |

**AGTF30 is a conceptual engine**: NASA's "Advanced Geared Turbofan 30,000 lbf", an
N+3 technology concept model published for research. It is **not** a production engine
and no result here describes any engine in service.

No NASA linear-model (A/B) matrices are available to compare against. Validation is
against (i) Table 1 of Chapman and Litt, NASA/TM-2017-219569 and (ii) the steady-state
initial conditions in the NASA repository itself. See *Validation status* below.

## Quick start

macOS (Apple silicon or Intel) and Linux:

```bash
./run.sh
```

Windows:

```bat
run.bat
```

The first run creates a private Python environment in `.venv/` (Python 3.13 preferred,
3.12 supported) and installs the packages listed in `requirements.txt`; this needs an
internet connection once. After that nothing is downloaded and the study runs offline.

`./run.sh` with no command does a quick verification and then starts the virtual lab in your
browser. On a fresh copy there are no results yet, so run the study first:

```bash
./run.sh full --fast     # coarse "skeleton" run, a short check that everything works
./run.sh full            # the complete study (sized for < 6 min on 2 CPU cores)
```

## The virtual lab (`./run.sh app`)

Opens at http://127.0.0.1:8050 in Safari/Chrome. Nothing is downloaded; it runs offline.

| Tab | What you do there |
|---|---|
| **Simulator** | 16 recommended experiments. Press **Run** on a card for a one-shot result (figures, numbers, what to look for). **Edit** copies the settings into the matching experiment tab. |
| **Live engine** | The engine in action: an animated cutaway (fan, booster, HPC, combustor flame, HPT, LPT, shafts, gas temperature, airflow, bleed valve, nozzle) driven by the live model. **Start / Stop** - it runs until you stop it. While it runs you can move the power lever or switch to manual fuel (with ±5 % fuel steps), change altitude / Mach / day temperature, set the nozzle and bleed valve by hand, wear any of the ten components, scale the controller gains, switch each NASA limiter on or off, change the simulation speed, play scripted scenarios, and press **Linearize now** to see the eigenvalues of A at the current state. |
| **Eigenvalue map, Jacobian accuracy, Fragility, Degradation, Gauss-Seidel, Linear vs nonlinear, Throttle transients** | One tab per experiment, every parameter editable, **Run** and **Reset to recommended**. |
| **Help** | What each tab does and how the model works, in plain words. |

Useful links: `?preset=burst` runs a Simulator preset, `?tab=live&scenario=climb` plays a live scenario,
`?tab=fragility&run=1` runs an experiment tab with its recommended settings.

## Commands

| Command | What it does |
|---|---|
| `./run.sh` | Quick verify, then start the virtual lab (`python -m twin.app --port 8050`) and open the browser. If verification fails the lab is not started (use `./run.sh app`). |
| `./run.sh verify` | Runs the test modules in `tests/` and prints the PASS/FAIL table from `outputs/data/*.json`. Target: under 60 s. |
| `./run.sh full [--fast] [--jobs N]` | Runs the pipeline `anchors`, `layer1_opmap`, then `steps`, `layer2_fd`, `layer3_bauerfike`, `layer4_degradation`, `validity` (N at a time, default 2), then `evidence.build`, then the PASS/FAIL table. `--fast` = coarse grids. |
| `./run.sh app` | Starts the virtual lab only (http://127.0.0.1:8050). `PORT=8051 ./run.sh app` changes the port. Stop it with Ctrl-C. |
| `./run.sh export` | Writes `outputs/engine_viewer.html` (offline replay of the animated engine during a throttle burst) and screenshots of the lab. |
| `./run.sh fetch-data [--yes]` | **Optional.** Tries to download NASA's N-CMAPSS data set into `data/` (several GB). The study does not need it; a failed download is not an error. |
| `./run.sh clean` | Deletes generated files (`outputs/` contents, `__pycache__`). Never touches source, `vendor/` or `.venv/`. |

The same commands work without the script, from the project folder:

```bash
python -m eigenmap.cli verify | full [--fast] | export | clean | hashes
python verify.py [--no-tests] [--strict]
```

Exit codes: `0` = everything passed, `1` = at least one stage or check failed (or no
results exist yet), `2` = unknown command. A pipeline module that has not been written
yet is reported as `SKIP` and is not a failure.

**Determinism check.** Two runs must give byte-identical data files:

```bash
./run.sh full && python -m eigenmap.cli hashes > run1.txt
EIGENMAP_OUT=/tmp/second_run ./run.sh full && EIGENMAP_OUT=/tmp/second_run python -m eigenmap.cli hashes > run2.txt
diff run1.txt run2.txt        # no output = identical
```

## Folder map

```
CONTRACT.md            binding notation, interfaces and tolerances of the project
run.sh / run.bat       one-command runners (macOS-Linux / Windows)
requirements.txt       Python packages
verify.py              gathers every check into one PASS/FAIL table
eigenmap/              the study (Python package)
  config.py              every tunable number: seeds, tolerances, grids, paths
  engine_model.py        PORT of the AGTF30 engine (two shaft states, 7 flow-matching residuals)
  steady_state.py        Newton-Raphson solver for steady operating points
  mathcore.py            own 2x2 linear algebra + finite differences (cross-checked with numpy/scipy)
  anchors.py             published turboshaft example (illustration only)
  linearize.py, opoints.py   reference Jacobians and the standard operating points
  layer1_opmap.py        Layer 1: eigenvalues and time constants across the operating map
  steps.py               step response and Gauss-Seidel steady state
  layer2_fd.py           Layer 2: finite-difference step h versus solver tolerance tol
  layer3_bauerfike.py    Layer 3: kappa(V) and the Bauer-Fike bound
  layer4_degradation.py  Layer 4: health parameters theta and eigenvalue movement
  validity.py            linear model versus nonlinear model for growing step sizes
  evidence/              builds every figure and table from outputs/data/*.json only
  cli.py                 the runner behind run.sh
replica/               independent second implementation of the engine (cross-check of the port)
twin/                  virtual lab extension: app.py (UI), sim.py (live simulation + controller),
                       experiments.py (7 experiments, 16 presets), assets/engine_canvas.js (animated engine)
tests/                 test modules (run by `verify`)
docs/equations.tex     every key equation as editable LaTeX
vendor/                NASA sources, unmodified, read-only (AGTF30 + T-MATS)
handoff/               working notes of the build
outputs/               GENERATED - never edit by hand
  data/*.json            one file per study module: numbers, checks, verdicts
  figures/*.png, *.svg   report figures (300 dpi, black-and-white safe)
  tables/*.csv, *.tex, tables.xlsx
  logs/*.log             console text of each pipeline stage
  verify_report.txt      the PASS/FAIL table as plain text
data/                  optional downloads (`fetch-data`); not used by the study
```

## How the outputs map to the report

Every figure and table is built from `outputs/data/*.json`, so each number in the report
can be traced to one file. *(Section names below are provisional; align them with the
report template.)*

| Report section | Take it from |
|---|---|
| Mathematical approach (equations) | `docs/equations.tex` (copy the LaTeX blocks, tags `E01` ... `E59` = the printed equation numbers) |
| Worked example by hand (turboshaft, illustration only) | `outputs/data/anchors.json` |
| Engine model, operating points, steady states | `outputs/data/layer1_opmap.json`, tables in `outputs/tables/` |
| Result 1: eigenvalue map, stability, time constants | `outputs/data/layer1_opmap.json`, `outputs/figures/` |
| Result 2: step response, Gauss-Seidel iteration | `outputs/data/steps.json` |
| Result 3: finite-difference step `h` versus solver tolerance `tol` | `outputs/data/layer2_fd.json` |
| Result 4: conditioning `kappa(V)`, Bauer-Fike bound | `outputs/data/layer3_bauerfike.json` |
| Result 5: degradation (health parameters `theta`) | `outputs/data/layer4_degradation.json` |
| Limits of the linear model | `outputs/data/validity.json` |
| Verification and validation table | `outputs/verify_report.txt` |
| Conclusions (hypothesis verdicts H1 to H4) | the `"verdicts"` entries of the JSON files, listed at the end of `outputs/verify_report.txt` |

Hypothesis verdicts (`supported`, `not supported`, `mixed`) are written by code from the
numbers, with the decision rule stored next to them. They are never typed by hand.

## Notation

`N_LP`, `N_HP` shaft-speed states (rpm; `N_LP = 3.1 N_fan`), `W_f` fuel flow (the input),
`theta` health parameters, `h` finite-difference step only, `delta` deviations from the
operating point, `lambda` eigenvalues, `tau = -1/Re(lambda)` time constants, `kappa(V)`
eigenvector condition number, `tol` inner-solver tolerance. Full table: `CONTRACT.md`,
section 2. Units are the English engineering units of the NASA source.

## Validation status

Numbers below are from Run 2 (full settings); `outputs/RESULTS_FACTS.md` is regenerated on every run.

- **Port vs NASA/TM-2017-219569 Table 1** (faithful `mode="nasa"`): airflow, TSFC and thrust within 0.32 %, TIT within 0.67 %
  at both published points; BPR and OPR within the printed rounding (worst 1.48 %: take-off BPR 27.59 vs printed "28").
  OPR is compared as the product of the compressor pressure ratios (ASSUMED NASA definition; as Pt3/Pt2 it is 2.3-2.7 % low).
- **No NASA A/B matrices exist** (no `nasa_lin.mat` was supplied and the repository has none), so the *linear* model is
  validated indirectly: against the nonlinear model (linear vs nonlinear fuel steps) and by self-consistency checks.
- **Port vs independent blind replica** (`replica/`, Gate 1): 140/140 checks, worst relative disagreement 8.5e-13. The replica
  caught one real porting bug (nozzle entropy when the inlet pressure is below ambient), fixed before the engine was frozen.
- **Own linear algebra vs numpy/scipy**: eigenvalues, inverse, Gauss-Seidel, diagonalization, e^{At}, step responses agree to
  better than 1e-8 at every operating point (`steps.json`, `layer1_opmap.json`, `tests/test_mathcore.py`).
- **H3**: the 'best step near 1e-3 at tol 1e-6' expectation fails at every point (best step ~1000x smaller); reported as such.
- **Linear vs nonlinear**: a 1 % fuel step agrees within 0.5-1.1 %; the linear model stays within 2 % up to 1-3 % fuel steps
  and within 5 % up to 3-10 % steps, depending on the operating point.
- **Determinism**: see `RUN_LOG.md` (Run 3) - two full runs into separate folders, hashes of every JSON and CSV compared.
- **Independent audit** (Gate 2): see `RUN_LOG.md` and `handoff/AUDIT_REPORT.md`.

## Known limitations

- AGTF30 is a conceptual engine; nothing here describes an engine in service.
- **More than one steady state at idle.** At sea-level idle the AGTF30 model (the faithful port as well) has a second steady
  state with the same corrected fan speed, VAFN and VBV but ~34 % less fuel and a different A matrix. The study uses the
  branch reached by continuation from max power along NASA's schedules; RESULTS_FACTS reports the other branch's eigenvalues
  and whether H1 changes. Reproduce any operating point by starting the solver from its tabulated z0.
- The study uses a differentiable version of the port (C1 bicubic map interpolation through NASA's map points, inner
  T-MATS loops converged to machine precision). It moves the published anchor quantities by at most 0.27 % relative to
  the faithful port, but it is a deliberate deviation (DERIVED, see ASSUMPTIONS.md E5).
- The linear model holds VAFN and VBV at their scheduled steady values; there are no actuator, sensor or controller
  dynamics in A and B.
- Layer 2's theoretical optimal-step formulas assume random noise in F. The measured solver error is not random (zero-iteration
  dead zone with a warm start; smooth, cancelling error from an offset start), so theory and experiment differ by 10-100x.
- Seed-to-seed spread of the Monte Carlo worst shift is below 1 %; a different perturbation distribution (e.g. structured
  finite-difference errors instead of Gaussian E) could give different worst cases.
- The published turboshaft example is printed to 1-2 significant figures; three of its time constants match the brief only
  to that printed rounding.
- The virtual lab is an extension: the engine drawing is schematic (not an engineering drawing; rotation is shown about
  30x slower than real), and the controller is our reading of NASA's baseline min/max-select PI structure in velocity form
  (gains and limits from NASA's files), not a block-by-block port of NASA's Simulink controller; actuators have simple rate
  limits and there are no sensor dynamics. "Max power" away from the two NASA anchor points is an assumed interpolation
  (ASSUMPTIONS L1). There is no start-up or sub-idle model: the live run pauses with a message if speeds leave the valid range.
- `run.bat` (Windows) and the macOS paths of `run.sh` were written for, but not executed on, those systems in the build sandbox.
- `./run.sh fetch-data` (N-CMAPSS) uses a URL that could not be verified from the build sandbox; the study does not need it.

## Licence and attribution

The folders `vendor/agtf30/` and `vendor/tmats/` contain unmodified copies of NASA's
open-source **AGTF30** (https://github.com/nasa/AGTF30) and **T-MATS**
(https://github.com/nasa/T-MATS) repositories; the exact commits are recorded in
`vendor/VENDOR_INFO.txt`. The licence files shipped with those copies,
`vendor/agtf30/LICENSE` and `vendor/tmats/LICENSE`, are the **Apache License,
Version 2.0**; they apply to the vendored files and to the engine data derived from them
(the parsed maps cached under `eigenmap/data/`). Read those two files for the binding
terms. This project is a student study. It is not affiliated with, or endorsed by, NASA.

Reference for the engine: J. W. Chapman and J. S. Litt, NASA/TM-2017-219569 (the source of
the Table 1 anchor values quoted in `CONTRACT.md`).
