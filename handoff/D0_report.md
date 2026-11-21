# D0 report - docs and runners

Files written (all D0-owned): `run.sh`, `run.bat`, `requirements.txt`, `README.md`, `verify.py`,
`eigenmap/cli.py`, `docs/equations.tex`, `handoff/D0_note.md`, this file. Nothing else was edited.
Read `handoff/D0_note.md` for the points that differ from the brief (licence is Apache 2.0, the
forward-difference error model, symbol clashes, behaviour choices).

## Acceptance results (sandbox: Linux, Python 3.13.18, uv present, 2 CPUs)

| Check | Result |
|---|---|
| `bash -n run.sh` | OK |
| `./run.sh verify` (called from /tmp): first run built `.venv` with uv + installed, 8.5 s total; second run skipped the install | runs end-to-end; **exit 1** at that moment, correctly: `tests.test_mathcore` printed 5 FAIL lines (agent M's work in progress) and no `outputs/data/*.json` exists yet |
| `./run.sh full --fast` (real tree) | runs end-to-end: 8 stages SKIP (none written yet), table printed, exit 1 for the same two reasons |
| Failing dummy check gives exit 1, removed gives exit 0 | shown below (section 4) with a relocated `EIGENMAP_OUT`, so no dummy file ever entered the project tree |
| Stage OK / FAIL / CRASH / SKIP handling, 2 workers, `evidence.build` last | shown below in a throwaway copy `/tmp/d0_sbx` with dummy stage modules (I did not create dummy modules under other agents' file names in the real tree) |
| `./run.sh clean` leaves source intact | 355 source files before and after, 0 removed; vendor and .venv sizes unchanged |
| `python -m eigenmap.cli hashes` | works; exit 1 with a clear message when there is nothing to hash |
| Determinism helper | two dummy runs into two output folders gave identical combined hashes |
| `docs/equations.tex` | compiles with pdflatex (2 passes, exit 0, 0 warnings, 0 overfull boxes), 6 pages, 59 equations; compiled in a temp dir, no aux files left |
| Non-uv bootstrap (`python3.13 -m venv` + pip), stamp skip, reinstall after requirements change | all three verified in the throwaway copy (14.8 s first run) |
| Default `./run.sh` with no `twin/app.py` | prints the warning and exits 0 after a passing verify (throwaway copy) |
| `./run.sh app` | launches `python -m twin.app --port N` (stand-in server answered HTTP 200); browser opener does not fail when `xdg-open` is absent |
| `.venv` | left in place in the project root (Python 3.13.18, stamp written) |

NOT verified: `run.bat` (no Windows here; written carefully, CRLF line endings, never executed);
macOS (script avoids bash-4 features, uses `shasum` fallback and `open`, but was only run on Linux);
the N-CMAPSS download URL (host blocked by the sandbox proxy); real `twin.app` / `twin.export_viewer`
/ `evidence.build` (do not exist yet); the < 60 s target of `verify` with all three real test modules
(only `test_mathcore` existed: about 2 s).

## Pasted output

### Stage handling in the throwaway copy (dummy stages: steps exits 1 with a FAIL line, layer2_fd raises)
```
== Study pipeline: FAST (coarse skeleton run), 2 worker(s) ==
  OK    anchors                  0.3 s   1 PASS / 0 FAIL lines
  OK    layer1_opmap             1.0 s   1 PASS / 0 FAIL lines
  CRASH layer2_fd                1.0 s   0 PASS / 0 FAIL lines, exit code 1
        | Traceback (most recent call last):
        |   File "/tmp/d0_sbx/eigenmap/layer2_fd.py", line 3, in <module>
        |     raise RuntimeError("boom in layer2")
        | RuntimeError: boom in layer2
  FAIL  steps                    1.0 s   0 PASS / 1 FAIL lines
        | FAIL steps: dummy claim (1.2e-9 vs 1e-8)
  SKIP  layer4_degradation   (module not written yet)
  OK    layer3_bauerfike         1.0 s   1 PASS / 0 FAIL lines
  OK    validity                 1.0 s   1 PASS / 0 FAIL lines
  SKIP  evidence.build       (module not written yet)
  stages: 4 ok, 1 failed checks, 1 crashed, 2 skipped   (total 3.4 s)
...
module           | check                                 |                         value |  limit | result
layer2_fd        | pipeline stage finished without error |                exit 1 (CRASH) | exit 0 | FAIL
steps            | pipeline stage finished without error |                 exit 1 (FAIL) | exit 0 | FAIL
tests            | tests.test_mathcore                   | exit 0, 1 PASS / 0 FAIL lines | exit 0 | PASS
anchors          | dummy claim of anchors                |                   1.23457e-09 |  1e-08 | PASS
steps            | dummy claim of steps                  |                   1.23457e-09 |  1e-08 | FAIL
...
11 passed, 3 failed
HYPOTHESIS VERDICTS  (written by code from numbers; rule stored in the JSON)
  H1_slower_at_low_power: supported   [layer1_opmap]
      rule: tau_slow(idle) > tau_slow(max)
  H3_best_h_shifts_with_tol: (no verdict yet - numbers not produced)
exit(full, with failing+crashing stage)=1
(after deleting the two bad dummies)  8 passed, 0 failed   exit=0
DETERMINISM: combined hashes identical (outputs/ vs EIGENMAP_OUT=/tmp/d0_sbx_out2)
```

### Real tree
```
== 3. ./run.sh full --fast (real tree, whatever exists now) ==
exit(run.sh full --fast)=1
== Study pipeline: FAST (coarse skeleton run), 2 worker(s) ==
  SKIP  anchors              (module not written yet)
  SKIP  layer1_opmap         (module not written yet)
  SKIP  layer2_fd            (module not written yet)
  SKIP  steps                (module not written yet)
  SKIP  layer3_bauerfike     (module not written yet)
  SKIP  validity             (module not written yet)
  SKIP  layer4_degradation   (module not written yet)
  SKIP  evidence.build       (module not written yet)
  stages: 0 ok, 0 failed checks, 0 crashed, 8 skipped   (total 0.0 s)

FAIL  tests.test_mathcore: exit 1, 62 PASS / 5 FAIL lines
      | FAIL test_mathcore: diagonalize raises DefectiveMatrixError exactly for the defective edge cases (False vs True)
      | FAIL test_mathcore: step_response vs scipy solve_ivp (DOP853, rtol 1e-12), 257 systems [bound 1e-9: integrator error] (2.015e-08 vs 1.000e-09)
      | FAIL test_mathcore: step_response accepts B as (2,), (2,1), (2,m); dx(0) = 0 (False vs True)
      | FAIL test_mathcore: richardson err_est is an upper bound of the true error (truncation regime) (False vs True)
      | FAIL test_mathcore: richardson_jacobian (h=1e-2 relative, levels=2) vs analytic Jacobian (1.829e-10 vs 1.000e-10)

VERIFICATION TABLE  (source: outputs/data/*.json + test modules)

module | check               |                          value |  limit | result
-------------------------------------------------------------------------------
tests  | tests.test_mathcore | exit 1, 62 PASS / 5 FAIL lines | exit 0 | FAIL
-------------------------------------------------------------------------------
== 2b. ./run.sh verify exit code ==
exit(run.sh verify)=1
== 4. dummy checks via relocated outputs (EIGENMAP_OUT), tests excluded ==
zz_dummy_ok | dummy passing check | 1e-09 | 1e-08 | PASS
1 passed, 0 failed
exit(only passing dummy)=0
zz_dummy_bad | dummy FAILING check |   0.5 | 1e-08 | FAIL
zz_dummy_ok  | dummy passing check | 1e-09 | 1e-08 | PASS
1 passed, 1 failed
exit(with failing dummy)=1
exit(failing dummy removed)=0
== 5. hashes ==
dd3980ab6a030f027d951fa76c891827cc8f96f8138bb004402525f496e284dd  data/zz_dummy_ok.json
fa37cc341d45a7b94c28f7769947461fa1922332d10ed898e7f47899c01e6d8f  COMBINED (1 files)
exit(hashes)=0
no outputs/data/*.json or outputs/tables/*.csv found (run the study first)
exit(hashes, real outputs: no data yet)=1
== 6. ./run.sh clean leaves source intact ==
  before: outputs/verify_report.txt
  size before: vendor 6024
  size before: .venv 462552
  cleaned outputs/ (1 item(s) removed)
  removed 3 __pycache__ folder(s)
  OK    clean (source, vendor/ and .venv untouched)
exit(run.sh clean)=0
  source files removed by clean: 0   (files: 355)
  my files + config + CONTRACT: checksums unchanged
  size after:  vendor 6024
  size after:  .venv 462552
  outputs/ now: []   __pycache__ left outside .venv: 0
== 7. default ./run.sh in a throwaway copy (dummy stages, EIGENMAP_PYTHON = system python) ==
exit(full --fast, 4 dummy stages ok)=0
[run] using EIGENMAP_PYTHON=/usr/local/bin/python3 (no virtual environment is created)
== Quick verification ==
OK    tests.test_mathcore: exit 0, 1 PASS / 0 FAIL lines
OK    tests.test_plain: exit 0, 2 PASS / 0 FAIL lines
10 passed, 0 failed
[ ok] verification passed
[warn] the 3D twin (twin/app.py) is not written yet - nothing to launch.
exit(./run.sh, no twin/app.py)=0
== 8. app launch with a stand-in twin/app.py (tiny http server), browser disabled ==
http status from the app on :8061 = 200
[run] using EIGENMAP_PYTHON=/usr/local/bin/python3 (no virtual environment is created)
[run] starting the twin at http://127.0.0.1:8061/   (stop it with Ctrl-C)
stand-in twin listening on 8061
127.0.0.1 - - [14/Nov/2026 07:21:02] "GET / HTTP/1.1" 200 -
exit(app, killed by timeout after 6 s)=124
== 8b. same with the browser opener enabled (no xdg-open here: must not fail) ==
exit(app + opener, timeout)=124
```
