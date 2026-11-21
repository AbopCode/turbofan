"""Command-line runner of the study:  python -m eigenmap.cli <command> [options]

Commands
    verify            quick check (< 60 s): run the test modules, print the PASS/FAIL table
    full [--fast]     run the whole pipeline (CONTRACT section 11 order), then the table
    export            write the stand-alone 3D viewer (twin/export_viewer.py), if it exists
    clean             delete generated files (outputs/ contents and __pycache__ folders)
    hashes            sha256 of outputs/data/*.json and outputs/tables/*.csv (determinism check)

Design rules
    * Every pipeline stage runs as `python -m eigenmap.<module>` in its OWN subprocess,
      so a crash in one stage is reported in the table and does not stop the others.
    * A stage whose source file does not exist yet prints SKIP (not an error).
    * `layer1_opmap` finishes before the stages that build on it; those then run
      `--jobs N` at a time (default 2 = the two CPUs we are sized for).
    * Elapsed times go to the screen only. Nothing with a timestamp is ever written
      into outputs/data, so two runs give byte-identical data files.
    * Only the Python standard library is imported here.
Course topic: (none - project plumbing)
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent           # project root

# ---- pipeline (CONTRACT section 11). Each entry: (stage name, source file under the root) ----
STAGES_FIRST = (("anchors", "eigenmap/anchors.py"),                 # published example
                ("layer1_opmap", "eigenmap/layer1_opmap.py"))       # operating-point map (needed by the rest)
STAGES_PARALLEL = (("validation", "eigenmap/validation.py"),          # port vs NASA Table 1 + Gate-1 summary
                   ("steps", "eigenmap/steps.py"),
                   ("layer2_fd", "eigenmap/layer2_fd.py"),
                   ("layer3_bauerfike", "eigenmap/layer3_bauerfike.py"),
                   ("layer4_degradation", "eigenmap/layer4_degradation.py"),
                   ("validity", "eigenmap/validity.py"))
STAGES_LAST = (("evidence.build", "eigenmap/evidence/build.py"),     # figures + tables from the JSON
               ("twin.export_viewer", "twin/export_viewer.py"))     # 3D twin: offline viewer + screenshots (extension)

STAGE_TIMEOUT_S = 3600        # safety net: no single stage may hang the whole run


def load_verify():
    """Import <root>/verify.py as a module (works from any working directory).

    Course topic: (none - project plumbing)
    """
    spec = importlib.util.spec_from_file_location("eigenmap_verify", ROOT / "verify.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


V = load_verify()             # shared helpers: paint(), child_env(), out_dir(), verify()


# --------------------------------------------------------------------------
# Running one stage
# --------------------------------------------------------------------------
def run_stage(name: str, source: str, fast: bool, timeout: int) -> dict:
    """Run `python -m eigenmap.<name> [--fast]`; never raises.

    Returns {'name', 'status': 'OK'|'FAIL'|'CRASH'|'SKIP', 'code', 'seconds', 'n_pass', 'n_fail', 'output'}.
    FAIL  = the module ran and reported failing checks (exit code 1 with FAIL lines).
    CRASH = anything else that is not exit code 0 (exception, timeout, ...).
    Course topic: (none - project plumbing)
    """
    result = {"name": name, "status": "SKIP", "code": None, "seconds": 0.0,
              "n_pass": 0, "n_fail": 0, "output": ""}
    if not (ROOT / source).is_file():
        return result
    module = name if name.startswith("twin.") else f"eigenmap.{name}"
    command = [sys.executable, "-m", module] + (["--fast"] if fast else [])
    start = time.perf_counter()
    try:
        proc = subprocess.run(command, cwd=ROOT, env=V.child_env(), text=True, timeout=timeout,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        code, output = proc.returncode, proc.stdout or ""
    except subprocess.TimeoutExpired as exc:
        partial = exc.stdout or ""
        partial = partial.decode(errors="replace") if isinstance(partial, bytes) else partial
        code, output = 124, partial + f"\n(stage timed out after {timeout} s)"
    except OSError as exc:
        code, output = 127, f"could not start the stage: {exc}"
    lines = output.splitlines()
    result.update(code=code, output=output, seconds=time.perf_counter() - start,
                  n_pass=sum(1 for ln in lines if ln.startswith("PASS ")),
                  n_fail=sum(1 for ln in lines if ln.startswith("FAIL ")))
    if code == 0:
        result["status"] = "OK"
    elif code == 1 and result["n_fail"] > 0:
        result["status"] = "FAIL"
    else:
        result["status"] = "CRASH"
    return result


def report_stage(res: dict, verbose: bool) -> None:
    """Print the one-line status of a finished stage (plus the reason if it failed).

    Course topic: (none - project plumbing)
    """
    colour = {"OK": "green", "SKIP": "yellow", "FAIL": "red", "CRASH": "red"}[res["status"]]
    tag = V.paint(f"{res['status']:<5}", colour)
    if res["status"] == "SKIP":
        print(f"  {tag} {res['name']:<20} (module not written yet)", flush=True)
        return
    detail = f"{res['n_pass']} PASS / {res['n_fail']} FAIL lines"
    if res["status"] == "CRASH":
        detail += f", exit code {res['code']}"
    print(f"  {tag} {res['name']:<20} {res['seconds']:7.1f} s   {detail}", flush=True)
    lines = res["output"].splitlines()
    if verbose:
        shown = lines
    elif res["status"] == "FAIL":
        shown = [ln for ln in lines if ln.startswith("FAIL ")]
    elif res["status"] == "CRASH":
        shown = lines[-25:]                                # the traceback is at the end
    else:
        shown = []
    for ln in shown:
        print("        | " + ln)


def save_log(res: dict) -> None:
    """Keep the raw console text of a stage in outputs/logs/<stage>.log (debugging aid).

    Only what the stage itself printed is stored: no times, no machine information.
    Course topic: (none - project plumbing)
    """
    if res["status"] == "SKIP":
        return
    try:
        logs = V.out_dir() / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        (logs / f"{res['name']}.log").write_text(res["output"], encoding="utf-8")
    except OSError:
        pass                                               # a log is a convenience, never fatal


# --------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------
def cmd_verify(args) -> int:
    """Quick check: test modules + PASS/FAIL table.

    Course topic: (none - project plumbing)
    """
    print(V.paint("== Quick verification ==", "bold"))
    return V.verify(run_test_modules=True, jobs=args.jobs)


def cmd_full(args) -> int:
    """Whole pipeline, then the verification table. Exit 1 if any stage or check failed.

    Course topic: (none - project plumbing)
    """
    mode = "FAST (coarse skeleton run)" if args.fast else "FULL"
    print(V.paint(f"== Study pipeline: {mode}, {args.jobs} worker(s) ==", "bold"))
    start = time.perf_counter()
    results: list[dict] = []

    def finish(res: dict) -> None:
        save_log(res)
        report_stage(res, args.verbose)
        results.append(res)

    for name, source in STAGES_FIRST:                      # one after the other
        finish(run_stage(name, source, args.fast, args.timeout))
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:   # independent stages
        futures = [pool.submit(run_stage, name, source, args.fast, args.timeout)
                   for name, source in STAGES_PARALLEL]
        for future in as_completed(futures):
            finish(future.result())
    for name, source in STAGES_LAST:                       # needs every JSON: runs last
        finish(run_stage(name, source, args.fast, args.timeout))

    counts = {s: sum(1 for r in results if r["status"] == s) for s in ("OK", "FAIL", "CRASH", "SKIP")}
    print(f"  stages: {counts['OK']} ok, {counts['FAIL']} failed checks, {counts['CRASH']} crashed, "
          f"{counts['SKIP']} skipped   (total {time.perf_counter() - start:.1f} s)\n")

    # A stage that did not exit cleanly becomes a FAIL row of the table. This matters because
    # a crashed stage can leave an OLD json behind whose checks all say PASS.
    bad_rows = [{"module": r["name"], "name": "pipeline stage finished without error", "pass": False,
                 "value": f"exit {r['code']} ({r['status']})", "limit": "exit 0"}
                for r in results if r["status"] in ("FAIL", "CRASH")]
    code = V.verify(run_test_modules=not args.no_tests, jobs=args.jobs, extra_rows=bad_rows)
    return 1 if (bad_rows or code) else 0


def cmd_export(args) -> int:
    """Write outputs/engine_viewer.html through twin/export_viewer.py (SKIP if absent).

    Course topic: (none - project plumbing)
    """
    if not (ROOT / "twin" / "export_viewer.py").is_file():
        print("  " + V.paint("SKIP ", "yellow") + " export: twin/export_viewer.py is not written yet")
        return 0
    start = time.perf_counter()
    code = subprocess.run([sys.executable, "-m", "twin.export_viewer"], cwd=ROOT, env=V.child_env()).returncode
    tag = V.paint("OK   ", "green") if code == 0 else V.paint("FAIL ", "red")
    print(f"  {tag} export {time.perf_counter() - start:7.1f} s   (exit code {code})")
    return 0 if code == 0 else 1


def _empty_folder(folder: Path, keep: tuple = (".gitkeep",)) -> int:
    """Delete everything inside `folder` (the folder itself stays). Returns the number removed."""
    removed = 0
    if not folder.is_dir():
        return 0
    for item in folder.iterdir():
        if item.name in keep:
            continue
        if item.is_dir() and not item.is_symlink():
            shutil.rmtree(item, ignore_errors=True)
        else:
            item.unlink(missing_ok=True)
        removed += 1
    return removed


def cmd_clean(args) -> int:
    """Remove generated files only: outputs/ contents and __pycache__ folders.

    Never touched: .venv, vendor, data (downloads), .git and every source file.
    Course topic: (none - project plumbing)
    """
    removed = _empty_folder(ROOT / "outputs")
    print(f"  cleaned outputs/ ({removed} item(s) removed)")

    # EIGENMAP_OUT may point somewhere else. Be careful there: only remove the
    # sub-folders/files this project creates, and refuse obviously dangerous targets.
    other = V.out_dir().resolve()
    if other != (ROOT / "outputs").resolve():
        unsafe = (other == ROOT.resolve() or other == Path.home().resolve()
                  or other.parent == other or ROOT.resolve().is_relative_to(other))
        if unsafe:
            print(V.paint(f"  refusing to clean EIGENMAP_OUT={other} (unsafe location)", "yellow"))
        else:
            for sub in ("data", "figures", "tables", "logs"):
                shutil.rmtree(other / sub, ignore_errors=True)
            for name in ("verify_report.txt", "engine_viewer.html"):
                (other / name).unlink(missing_ok=True)
            print(f"  cleaned generated folders in {other}")

    skip = {".venv", "vendor", ".git", "data", "node_modules"}
    n_cache = 0
    for current, folders, _files in os.walk(ROOT):
        folders[:] = [f for f in folders if not (Path(current) == ROOT and f in skip)]
        for f in list(folders):
            if f in ("__pycache__", ".pytest_cache"):
                shutil.rmtree(Path(current) / f, ignore_errors=True)
                folders.remove(f)
                n_cache += 1
    print(f"  removed {n_cache} __pycache__ folder(s)")
    print("  " + V.paint("OK   ", "green") + " clean (source, vendor/ and .venv untouched)")
    return 0


def cmd_hashes(args) -> int:
    """Print sha256 of every data JSON and table CSV, plus one combined hash.

    Determinism check: run the study twice (e.g. the second time with EIGENMAP_OUT
    pointing to another folder) and compare the combined hashes; they must be equal.
    Course topic: (none - project plumbing)
    """
    out = V.out_dir()
    files = sorted(list((out / "data").glob("*.json")) + list((out / "tables").glob("*.csv")),
                   key=lambda p: p.relative_to(out).as_posix())
    if not files:
        print("no outputs/data/*.json or outputs/tables/*.csv found (run the study first)")
        return 1
    combined = hashlib.sha256()
    for path in files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rel = path.relative_to(out).as_posix()
        combined.update(f"{digest}  {rel}\n".encode())     # combined hash covers names AND contents
        print(f"{digest}  {rel}")
    print(f"{combined.hexdigest()}  COMBINED ({len(files)} files)")
    return 0


def main(argv: list[str] | None = None) -> int:
    """Parse the command line and dispatch.

    Course topic: (none - project plumbing)
    """
    parser = argparse.ArgumentParser(prog="python -m eigenmap.cli",
                                     description="Runner of the turbofan eigenvalue-map study.")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("verify", help="run the tests and print the PASS/FAIL table (< 60 s)")
    p.add_argument("--jobs", type=int, default=2, help="test modules run at once (default 2)")
    p.add_argument("--fast", action="store_true", help="accepted for symmetry; verify is always quick")
    p.set_defaults(func=cmd_verify)

    p = sub.add_parser("full", help="run the whole pipeline, then the table")
    p.add_argument("--fast", action="store_true", help="coarse skeleton run (fewer points/samples)")
    p.add_argument("--jobs", type=int, default=2, help="independent stages run at once (default 2)")
    p.add_argument("--timeout", type=int, default=STAGE_TIMEOUT_S, help="seconds allowed per stage")
    p.add_argument("--no-tests", action="store_true", help="skip tests/test_*.py in the final table")
    p.add_argument("--verbose", action="store_true", help="print every line each stage printed")
    p.set_defaults(func=cmd_full)

    for name, func, text in (("export", cmd_export, "write the stand-alone 3D viewer"),
                             ("clean", cmd_clean, "delete generated files"),
                             ("hashes", cmd_hashes, "sha256 of the data files (determinism check)")):
        p = sub.add_parser(name, help=text)
        p.add_argument("--fast", action="store_true", help=argparse.SUPPRESS)
        p.set_defaults(func=func)

    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
