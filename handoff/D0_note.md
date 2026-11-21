# D0 notes for the orchestrator (things that differ from the brief, or need a decision)

1. **Licence of the vendored NASA sources is Apache 2.0, not NOSA.**
   `vendor/agtf30/LICENSE` and `vendor/tmats/LICENSE` are both "Apache License, Version 2.0"
   (201 lines each); `vendor/tmats/README.md` line 35 says "released under the Apache V2.0
   license agreement". No file under those two folders mentions the NASA Open Source Agreement.
   README.md therefore says Apache 2.0. Do not write "NASA Open Source Agreement" in the report
   unless you have another source for it.

2. **Forward-difference error model: 2 eps/h, not eps/h.**
   The brief gave `eps/h + (M2/2) h` together with `h* = 2 sqrt(eps/M2)`. Those two do not match
   (minimising eps/h + (M2/2)h gives sqrt(2 eps/M2)). `docs/equations.tex` uses
   `2 eps/h + (M2/2) h` -> `h* = 2 sqrt(eps/M2)` (two function values, each wrong by eps), which
   is also what `eigenmap/mathcore.py::optimal_h_forward` documents. Central: `eps/h + (M3/6) h^2`
   -> `h* = (3 eps/M3)^(1/3)` as in the brief.

3. **Symbol clashes inside the CONTRACT notation** (kept as in the contract, flagged in comments):
   - `D` is both the eigenvalue matrix (A = P D P^-1, E25) and the scaling matrix diag(x0)
     (A_n = D^-1 A D, E51). Suggest `D_x` for the scaling matrix in the report.
   - The brief's sensitivity formula used `x` for the right eigenvector; `x` is the state vector,
     so E56 uses `v_i` (right) and `y_i` (left).
   - Gauss-Seidel stopping tolerance is written eps_GS so that `tol` stays the inner-solver tolerance.

4. **N-CMAPSS URL still unverified.** From this sandbox the host is blocked by the proxy
   (HTTP 403, "no rule allows host phm-datasets.s3.amazonaws.com"). `fetch-data` says so, prints
   the repository page, asks for confirmation, and always exits 0.

5. **Behaviour choices I made (change if you disagree):**
   - `verify` exits 1 when there is no `outputs/data/*.json` (as specified). So on a fresh copy
     `./run.sh` stops with "run ./run.sh full first"; it does NOT auto-run the study (a silent
     `--fast` run could be mistaken for the real results).
   - Default `./run.sh`: if verification fails the twin is not started (hint: `./run.sh app`).
   - `full` also runs the test modules in the final table (`--no-tests` skips them).
   - A stage that exits non-zero adds a FAIL row "pipeline stage finished without error" to the
     table, because a crashed stage can leave an old all-PASS JSON behind.
   - A JSON without a "checks" list is a NOTE, not a failure (could be a data cache);
     `python verify.py --strict` additionally fails when a pipeline module has no JSON.
   - Test discovery is `tests/test_*.py` only, so `tests/gate1_compare.py` is not run by `verify`.
     A test file with a `__main__` block runs as `python -m tests.<name>`; otherwise a small
     built-in runner calls its `test_*` functions (pytest is not a dependency).
   - Children get `MPLBACKEND=Agg`, `PYTHONHASHSEED=0`, one BLAS thread (all `setdefault`).
   - Raw stage console text is saved to `outputs/logs/<stage>.log` (no times added by the CLI).
   - `clean` always empties `<root>/outputs`; a relocated `EIGENMAP_OUT` only loses its
     data/figures/tables/logs folders and is refused if it is `/`, home, or contains the project.
   - Extra, not in the brief: `EIGENMAP_PYTHON=/path/python ./run.sh ...` skips the venv.

6. **README**: "How the outputs map to the report" uses provisional section names; I could not
   open the report template. The reference line gives authors + NASA/TM-2017-219569 only (title
   not confirmed from any file here). "Validation status" / "Known limitations" are placeholders.
