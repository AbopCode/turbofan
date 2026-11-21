# CONTRACT.md — "The Eigenvalue Map of a Turbofan Engine"

Owner: the ORCHESTRATOR. Sub-agents may NOT edit this file. If something here is
wrong or blocks you, write a note to `handoff/<your-id>_note.md` and carry on with
the most reasonable reading; the orchestrator decides.

Project root: `/home/claude/turbofan_eigenmap` (called `$P`). Python 3.13 (must also
run on 3.12). Allowed packages only: numpy, scipy, pandas, matplotlib, openpyxl,
sympy, plotly, dash, (h5py optional). No network at run time. No CDNs.

## 0. What this project is (and is not)

A Python port of NASA's open-source AGTF30 geared-turbofan simulation (a
*conceptual* N+3 30,000 lbf engine; never call it a production engine), a
linearization/eigenvalue study of it, and a separate 3D "digital twin" extension.
Honesty rules: say what is PORTED, what is DERIVED, what is ASSUMED; never write a
hypothesis verdict before its numbers exist; report contradictions.

No NASA A/B matrices are available (`nasa_truth/nasa_lin.mat` was not supplied and
the repo has none). Validation is against (i) Table 1 of Chapman & Litt,
NASA/TM-2017-219569 and (ii) the repo's own steady-state initial conditions.

Published anchors (opened and confirmed by the orchestrator on 2026-11-14):

| Point | Alt ft | MN | dT | Airflow lbm/s | TSFC lbm/hr/lbf | BPR | OPR | TIT R | Fn lbf |
|---|---|---|---|---|---|---|---|---|---|
| Cruise | 35,000 | 0.8 | 0 | 813 | 0.4637 | 24 | 55 | 3150 | 6073 |
| Max thrust, hot-day SLS | 0 | 0 | +27 F | 1724 | 0.1751 | 28 | 37 | 3170 | 28622 |

Repo steady-state ICs (`vendor/agtf30/AGTF30/SimSetup/setup_InputsSS.m`), near-truth:
cruise: Nc_fan = 1.0*MWS.FAN.s_Nc, Wf 0.7837, N2 6772.0, N3 20796, VAFN 4775.22, VBV 0,
ICcom = [815.48, 4.0955, 1.994, 2.0007, 23.19937, 2.19956, 11.0178];
takeoff: Nc_fan = 0.810*s_Nc, dT 27, Wf 1.3896, N2 6079.0, N3 21583.2, VAFN 6314.82, VBV 0,
ICcom = [1723.81, 4.17, 2.0325, 1.75, 27.5078, 1.8740, 7.15].

## 1. Directory layout and file ownership (disjoint)

```
$P/CONTRACT.md                     orchestrator
$P/eigenmap/config.py              orchestrator  (StudyConfig, seeds, paths)
$P/eigenmap/engine_model.py        A   (+ eigenmap/tmats_thermo.py, eigenmap/agtf30_data.py, eigenmap/data/*)
$P/eigenmap/steady_state.py        A
$P/tests/test_engine.py            A
$P/replica/**                      A2  (works in /home/claude/a2_work, copied in at Gate 1)
$P/eigenmap/mathcore.py            M
$P/eigenmap/anchors.py             M   (published turboshaft example)
$P/tests/test_mathcore.py, tests/test_anchors.py   M
$P/docs/equations.tex, README.md, run.sh, run.bat, requirements.txt, eigenmap/cli.py, verify.py   D0
$P/eigenmap/linearize.py, eigenmap/opoints.py      orchestrator (written at Gate 1)
$P/eigenmap/layer1_opmap.py, eigenmap/steps.py     L  (Layer 1 + items 6: step response, Gauss-Seidel)
$P/eigenmap/layer2_fd.py, eigenmap/layer3_bauerfike.py   N
$P/eigenmap/layer4_degradation.py, eigenmap/validity.py  G
$P/eigenmap/evidence/**, outputs/*.md              E
$P/twin/**                         T
$P/handoff/**                      everyone (notes only; file name starts with your id)
$P/vendor/**                       read-only NASA sources (AGTF30 + T-MATS library)
$P/outputs/**                      generated; never hand-edited
```

Agents talk to each other only through `handoff/`. An agent edits only what it owns.

## 2. Notation table (use these symbols EVERYWHERE: code comments, LaTeX, tables)

| Symbol | Code name | Meaning | Unit |
|---|---|---|---|
| N_LP | `x[0]` | LP-shaft speed STATE = LPC/LPT shaft. AGTF30 *code* calls it `N2`; the Chapman-Litt *paper* calls it N1. N_LP = 3.1 x N_fan | rpm |
| N_fan | `N_fan` | fan-shaft speed (code "N1", paper "Nf") = N_LP / GR, GR = 3.1 | rpm |
| N_HP | `x[1]` | HP-shaft (core) speed STATE. Code "N3", paper "N2" | rpm |
| x | `x` | state vector [N_LP, N_HP]^T | rpm |
| u | `Wf` | fuel flow W_f (the only input of the linear model) | lbm/s |
| theta | `theta` | health parameters: multiplicative scalars on map efficiency and map flow capacity, 1.0 = healthy | - |
| h | `h` | finite-difference step ONLY (relative: perturbation = h * |x0_i| or h * |u0|) | - |
| delta | `d...` | deviation from the operating point, dx = x - x0, du = u - u0 | |
| F | `F` | F(x,u;theta) = dx/dt | rpm/s |
| A, B | `A`, `B` | A = dF/dx (2x2), B = dF/du (2x1) at F(x0,u0)=0 | 1/s, rpm/(lbm/s) |
| lambda | `lam` | eigenvalues of A | 1/s |
| tau | `tau` | time constant -1/Re(lambda); tau_slow >= tau_fast | s |
| V, kappa(V) | `V`, `kappaV` | eigenvector matrix (unit 2-norm columns), its 2-norm condition number | - |
| tol | `tol` | inner (flow-matching) Newton solver tolerance: stop when max|r_i| < tol | - |
| z | `z` | 7 flow-matching independents of the inner solver (NASA order, sec. 4) | mixed |
| r | `r` | 7 normalized flow-matching residuals (as T-MATS defines them) | - |
| I_LP, I_HP | `I_LP`, `I_HP` | shaft inertias; I_LP includes the fan reflected through the gear | slug ft^2 |
| E | `E` | perturbation matrix added to A (Bauer-Fike study) | 1/s |

Shaft equations (T-MATS Shaft block; verify against the source):
I_LP (pi/30) dN_LP/dt = net LP torque;  I_HP (pi/30) dN_HP/dt = net HP torque.
I_LP = (450102/3.1^2 + 17796 + 16237)/32.2/144 slug ft^2, I_HP = 8627/32.2/144.

Units are English engineering units as in T-MATS: lbm/s, psia, deg R, BTU/lbm,
lbf, ft lbf, in^2, rpm, hp where the source uses hp. Never convert silently.

Normalized-state form (used for kappa(V) and Layer 3): with D = diag(x0),
A_n = D^-1 A D (states in fraction-of-operating-speed). Eigenvalues are identical;
V and kappa(V) are NOT. Always say which form a kappa(V) is in. Primary = normalized.

## 3. Health parameters theta

`theta: dict[str, float]`, keys exactly
`fan_eff, fan_flow, lpc_eff, lpc_flow, hpc_eff, hpc_flow, hpt_eff, hpt_flow, lpt_eff, lpt_flow`
(lpc = booster). Missing keys = 1.0. `X_eff` multiplies the component's map
efficiency after the map scalar (Eff = EffMap * s_Eff * theta_eff); `X_flow`
multiplies the map corrected flow after the map scalar (Wc = WcMap * s_Wc * theta_flow).
"3% degradation" means theta = 0.97.

## 4. Engine interface (A and A2 implement the SAME interface)

A: `from eigenmap.engine_model import load_engine, Env, Controls`
A2: `from replica.engine_replica import load_engine, Env, Controls`

```python
@dataclass(frozen=True)
class Env:      alt_ft: float; MN: float; dT: float = 0.0     # dT = delta ambient temp, deg F (= deg R), as MWS.In.dTambIC
@dataclass(frozen=True)
class Controls: Wf: float; VAFN: float; VBV: float            # lbm/s ; bypass-nozzle throat area in^2 ; 0 = closed .. 1 = open

load_engine(vendor_dir: str | None = None) -> Engine          # parses the REAL maps from vendor/agtf30 (CR line endings!)

Engine.GR, Engine.I_LP, Engine.I_HP                           # floats
Engine.residuals(z, x, ctrl, env, theta=None) -> (r, out)
    # z: sequence of 7 floats in NASA "ICcom" order:
    #   [W (inlet airflow lbm/s), HPT_PR, HPC_Rline, FAN_Rline, BPR, LPC_Rline, LPT_PR]
    # x: [N_LP, N_HP] rpm.   r: numpy array (7,), normalized residuals in the order of
    #   the dynamic model's dependents in AGTF30SysDyn.mdl / AGTF30_eng.mdl (document the order in your notes).
    # out: dict of floats with AT LEAST these keys:
    #   'Trq_LP' net torque on the LP shaft, ft lbf (LPT + LPC + fan-through-gearbox, sign so that >0 accelerates)
    #   'Trq_HP' net torque on the HP shaft, ft lbf
    #   'Ndot_LP', 'Ndot_HP'  rpm/s
    #   'Fnet' lbf, 'Fgross' lbf, 'W_in' lbm/s, 'TSFC' lbm/hr/lbf, 'BPR', 'OPR' (= Pt3/Pt2),
    #   'T4' deg R (burner exit total), 'T45' deg R, 'Ps3' psia, 'Pt2','Tt2','Pt25','Tt25','Pt3','Tt3',
    #   'N_fan' rpm, 'Nc_fan' rpm (fan corrected speed = N_fan/sqrt(Tt2/518.67)),
    #   'SM_fan','SM_lpc','SM_hpc' (stall margins, %, as T-MATS computes them),
    #   'W25' HPC inlet flow, 'W_byp' bypass flow, 'W_vbv' VBV bleed flow.
Engine.solve_flow(x, ctrl, env, theta=None, tol=1e-12, z0=None, max_iter=60, jac=None) -> (z, out, info)
    # Inner Newton(-chord) solve of r(z)=0 for fixed x. Stops when max|r| < tol (the ONLY stopping rule;
    # no extra polishing, so tol really is the residual left behind). info = {'converged': bool, 'iters': int,
    # 'resnorm': float, 'nfev': int}. `jac` = optional precomputed 7x7 dr/dz to reuse (chord iterations);
    # default: finite-difference Jacobian at the start, refreshed if convergence stalls. Deterministic.
Engine.F(x, ctrl, env, theta=None, tol=1e-12, z0=None, jac=None) -> numpy array (2,)  # [dN_LP/dt, dN_HP/dt] rpm/s
Engine.schedule_VAFN(Nc_fan, MN) -> float ; Engine.schedule_VBV(Nc_fan, MN) -> float
    # NASA baseline schedules MWS.Cntrl.VAFN_sch / VBV_sch (2-D linear interpolation, clamped at table edges)
```

steady state (A: `eigenmap.steady_state`, A2: `replica.steady_replica`):

```python
@dataclass
class SteadyPoint: env: Env; ctrl: Controls; x0: np.ndarray; z0: np.ndarray; out: dict; theta: dict; info: dict

solve_steady(engine, env, *, target: tuple[str, float], theta=None, VAFN=None, VBV=None,
             tol=1e-12, guess: SteadyPoint | None = None) -> SteadyPoint
    # Newton-Raphson on the 9 unknowns [z(7), N_LP, N_HP] with 7 flow residuals + 2 torque balances,
    # for target = ('Wf', value). For any other target key (('Nc_fan', v), ('N_fan', v), ('Ps3', v), ('Fnet', v),
    # ... any key of `out`) Wf becomes a 10th unknown. VAFN/VBV = None -> NASA schedules evaluated at the
    # solution's own Nc_fan and MN (solved consistently); a float -> held fixed.
    # Converged when max|r| < tol AND |Ndot_i|/N_i < tol (1/s). info has 'converged', 'iters', 'resnorm'.
```

Performance targets (the whole study must run in < 6 min on 2 CPUs): one `residuals`
call <= ~1 ms (scalar `math`, bisect table lookups, no pandas, no per-call numpy allocations in
hot paths where avoidable), one warm-started `F` call <= ~20 ms. Maps and gas tables are parsed
once (cache parsed data as .npz/.json under `eigenmap/data/`, regenerated from vendor if missing).

Fidelity: interpolation, gas-property tables, component equations, residual definitions and
the solver structure come from the T-MATS / AGTF30 SOURCE (C files in `vendor/tmats/TMATS_Library/MEX`,
`.slx` = zip of XML, `.mdl` = text), not from memory. Record every place where you could not
follow the source in `handoff/<id>_deviations.md` with a PUBLIC / DERIVED / ASSUMED tag.

## 5. Math core interface (M) — `eigenmap/mathcore.py`, engine-independent

```python
eig2x2(A) -> (lam: complex ndarray (2,), info: dict)      # OWN characteristic-equation formula; lam sorted: slow first
     # info: {'tr','det','disc'}; real pair sorted by |Re| ascending; complex pair: +imag first
eigvec2x2(A, lam) -> V (2x2, unit 2-norm columns)           # OWN formula
diagonalize(A) -> (P, D, Pinv)                              # A = P D P^-1, own code, complex-safe
kappa2(V) -> float                                          # 2-norm condition number (own 2x2 SVD formula, cross-checked)
is_stable2x2(A) -> bool                                     # tr<0 and det>0
time_constants(lam) -> (tau_slow, tau_fast)
expm_diag(A, t) -> ndarray                                  # e^{At} via own diagonalization (fallback for defective A documented)
step_response(A, B, du, t) -> dx(t) array (len(t), 2)       # dx(t) = A^-1 (e^{At} - I) B du
steady_gain(A, B, du) -> dx_ss                              # -A^-1 B du (own 2x2 inverse)
gauss_seidel(A, b, x0=None, tol=1e-10, max_iter=500) -> (x, history: list[ndarray], info)   # OWN loop
gs_spectral_radius_2x2(A) -> float                          # |a12 a21 / (a11 a22)|
is_diag_dominant(A) -> bool
fd_jacobian(func, x0, h, scheme='forward'|'central', rel=True) -> J    # func: R^n -> R^m ; step_i = h*max(|x0_i|, floor)
richardson_jacobian(func, x0, h, levels=2) -> (J, err_est)  # central + Richardson extrapolation
complex_step_jacobian(func, x0, h=1e-30) -> J               # for complex-safe analytic funcs only
optimal_h_forward(eps, M2, scale=1.0), optimal_h_central(eps, M3, scale=1.0)
bauer_fike_bound(V, E) -> float                             # kappa2(V) * ||E||_2
random_perturbations(rng, n, norm2) -> ndarray (n,2,2)      # each with ||E||_2 == norm2 exactly
eig_shift(A, E) -> float                                    # max_i min_j |mu_i - lambda_j|  (mu = eig(A+E))
```
Every own routine is cross-checked against numpy/scipy in `tests/test_mathcore.py` to 1e-10.

## 6. Linearization and operating points (orchestrator, written at Gate 1)

```python
eigenmap.linearize.reference_AB(engine, sp, tol=cfg.tol_tight) -> dict   # 'A','B','A_n','err_est', method notes
     # central difference + Richardson at tight tol, steps relative to x0/u0, VAFN & VBV FROZEN at sp.ctrl.
eigenmap.linearize.fd_AB(engine, sp, h, scheme, tol) -> (A, B)
eigenmap.opoints.standard_points(engine, cfg) -> list[OpPoint]            # cached; id, label, group ('SLS'|'CRUISE'|'ANCHOR'), sp
```
Linear model held-constant list (documented in ASSUMPTIONS.md): VAFN, VBV at their steady
scheduled values; no actuator, sensor, heat-soak or volume dynamics (T-MATS AGTF30 has no volume
dynamics: flow matching is algebraic, solved every step); gearbox rigid (fan inertia reflected).

## 7. Study configuration — `eigenmap/config.py` (orchestrator)

`StudyConfig(fast: bool)` gives every tunable: seeds, tolerances, grids, paths. Agents READ it,
never hard-code these numbers, never edit the file (ask via handoff). Deterministic seeds:
base seed 20261114; Monte Carlo seeds (20261114, 20261115, 20261116).
`tol` values for Layer 2: (1e-13, 1e-9, 1e-6, 1e-4). h grid: 1e-12 .. 1e-1.

## 8. Module run convention and data files

Every study module exposes `run(cfg: StudyConfig) -> dict` and is runnable as
`python -m eigenmap.<module> [--fast]`. It writes exactly one JSON:
`outputs/data/<module>.json` (plain floats/lists, `json.dump(..., indent=1, sort_keys=True)`,
complex numbers as `[re, im]`, NO timestamps, NO absolute paths, NO machine info — two runs must
give byte-identical files). Console: one line per check, starting with `PASS ` or `FAIL `
(format: `PASS <module>: <claim> (<number> vs <limit>)`), and exit code 1 if any FAIL.
The evidence agent builds all figures/tables ONLY from `outputs/data/*.json`.

Figures (later, agent E): `outputs/figures/<name>.png` (300 dpi) + `.svg`, matplotlib Agg,
black-and-white safe (distinct markers + linestyles), labelled axes with units, deterministic
(`svg.hashsalt` set, no dates in metadata). Tables: `outputs/tables/<name>.csv` + `.tex`, and one
`outputs/tables/tables.xlsx`.

## 9. PASS/FAIL conventions and tolerances

- Published turboshaft example (always labelled "turboshaft JetCat SPT5, illustration only"): 
  eigenvalues/tau within 2e-3 absolute of the brief's values; kappa within 1%; Gauss-Seidel table to 5 d.p.
- Own vs library math: 1e-8 (relative or absolute, whichever is looser) at every operating point.
- Steady state: max|r| < 1e-8 and |Ndot|/N < 1e-8 at every point (default solver tol 1e-12).
- Reference Jacobian self-consistent to 6 digits (two independent step sizes).
- Port vs replica (Gate 1): relative 1e-6 on torques, residuals, steady points.
- Port vs NASA Table 1: 1% (real-map path).
- Linear vs nonlinear, 1% fuel step: within 2% (of the peak deviation).
- Stability: diag(A) < 0 and Re(lambda) < 0 at every point, else it is a bug until proven physical.

Verdict words for hypotheses: `supported`, `not supported`, `mixed`. A verdict is written by code from
numbers with an explicit rule stored next to it in the JSON (`"rule": "..."`), never by hand.

## 10. Code style

Type hints, short functions, comments a first-year student can follow, each function docstring
ends with a line `Course topic: <topic>` (e.g. "Jacobians / partial derivatives", "Eigenvalues",
"Gauss-Seidel", "Linear ODE systems", "Taylor series", "Stability", or "Engine model (not a course topic)").

## 11. Checks list, pipeline order, CLI (amendment written before Wave 1 started)

Every `outputs/data/<module>.json` has a top-level key `"checks"`: a list of
`{"name": str, "pass": bool, "value": number|str, "limit": number|str}` — the same items the module
printed as PASS/FAIL lines. `verify.py` gathers them into one table. Hypothesis verdicts live under
a top-level key `"verdicts"`: `{"<hypothesis id>": {"verdict": "supported|not supported|mixed", "rule": str, "numbers": {...}}}`
with ids `H1_slower_at_low_power`, `H2_idle_most_fragile`, `H3_best_h_shifts_with_tol`, `H4_degradation_direction`.

Pipeline order (module names under `eigenmap/`): `anchors`, `layer1_opmap`, `steps`, `layer2_fd`,
`layer3_bauerfike`, `layer4_degradation`, `validity`, then `evidence.build`. `eigenmap/cli.py`
(`python -m eigenmap.cli <verify|full|export|clean> [--fast]`) runs them, tolerates a module that does
not exist yet (prints `SKIP`), returns nonzero if any check fails. Stages that do not depend on each
other may run in 2 worker processes. The twin lives in `twin/app.py` (`python -m twin.app --port 8050`)
and `twin/export_viewer.py` (writes `outputs/engine_viewer.html`).
