#!/usr/bin/env python3
"""verify.py - gather every PASS/FAIL check of the study into one table.

Usage (from anywhere):
    python verify.py              # run the test modules, then print the table
    python verify.py --no-tests   # only read outputs/data/*.json
    python verify.py --strict     # also FAIL when a pipeline module has no JSON yet

What it does
    1. (unless --no-tests) runs every tests/test_*.py in its own subprocess.
    2. Reads every outputs/data/*.json and collects
         "checks":   [{"name", "pass", "value", "limit"}, ...]      (CONTRACT sec. 11)
         "verdicts": {"H1_...": {"verdict", "rule", "numbers"}, ...}
    3. Prints an aligned table  module | check | value | limit | PASS/FAIL,
       a summary line "N passed, M failed" and the hypothesis verdicts.
    4. Writes the same text (no colours, no timestamps) to outputs/verify_report.txt.
    5. Exit code 1 if any check failed or if no JSON file exists, else 0.

Only the Python standard library is used, so this file works before the
study modules exist and even outside the virtual environment.
Course topic: (none - project plumbing)
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

# Pipeline modules in CONTRACT section 11 order (their JSON is outputs/data/<name>.json).
PIPELINE = ("anchors", "layer1_opmap", "steps", "layer2_fd",
            "layer3_bauerfike", "layer4_degradation", "validity")

# Hypothesis ids fixed by CONTRACT section 11.
HYPOTHESES = ("H1_slower_at_low_power", "H2_idle_most_fragile",
              "H3_best_h_shifts_with_tol", "H4_degradation_direction")

# Tests named in the brief come first; any other tests/test_*.py follows alphabetically.
KNOWN_TESTS = ("test_mathcore", "test_anchors", "test_engine")

TEST_TIMEOUT_S = 600          # one test module may not run longer than this


def out_dir() -> Path:
    """Folder that holds the generated outputs.

    Same rule as eigenmap/config.py: the EIGENMAP_OUT environment variable
    relocates it (used for clean-room / determinism runs), default <root>/outputs.
    Course topic: (none - project plumbing)
    """
    return Path(os.environ.get("EIGENMAP_OUT", ROOT / "outputs"))


# --------------------------------------------------------------------------
# Colours (only when writing to a real terminal)
# --------------------------------------------------------------------------
_CODES = {"green": "32", "red": "31", "yellow": "33", "cyan": "36", "bold": "1", "dim": "2"}


def use_colour(stream=None) -> bool:
    """True when `stream` is an interactive terminal and NO_COLOR is not set.

    Course topic: (none - project plumbing)
    """
    stream = sys.stdout if stream is None else stream
    if os.environ.get("NO_COLOR"):
        return False
    return bool(getattr(stream, "isatty", lambda: False)())


def paint(text: str, colour: str, enabled: bool | None = None) -> str:
    """Wrap `text` in an ANSI colour code (or return it unchanged).

    Course topic: (none - project plumbing)
    """
    enabled = use_colour() if enabled is None else enabled
    if not enabled or colour not in _CODES:
        return text
    return f"\033[{_CODES[colour]}m{text}\033[0m"


# --------------------------------------------------------------------------
# Running the test modules
# --------------------------------------------------------------------------
def find_test_modules() -> list[str]:
    """Names (without .py) of the test modules that exist right now.

    Course topic: (none - project plumbing)
    """
    tests = ROOT / "tests"
    found = sorted(p.stem for p in tests.glob("test_*.py")) if tests.is_dir() else []
    first = [name for name in KNOWN_TESTS if name in found]
    return first + [name for name in found if name not in first]


def child_env() -> dict:
    """Environment for every subprocess: headless plots, repeatable hashing,
    one BLAS thread per process (we already run 2 processes on 2 CPUs).

    Course topic: (none - project plumbing)
    """
    env = dict(os.environ)
    env.setdefault("MPLBACKEND", "Agg")
    env.setdefault("PYTHONHASHSEED", "0")
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env.setdefault(var, "1")
    # Make `import eigenmap`, `import tests`, `import twin` work from any cwd.
    env["PYTHONPATH"] = os.pathsep.join(
        [str(ROOT)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
    return env


def _test_command(name: str) -> list[str]:
    """Choose how to run tests/<name>.py.

    * a file with an `if __name__ == "__main__"` block is a script: `python -m tests.<name>`;
    * otherwise the tiny built-in runner below calls its test_* functions
      (pytest is NOT a project dependency, so we never require it).
    Course topic: (none - project plumbing)
    """
    text = (ROOT / "tests" / f"{name}.py").read_text(encoding="utf-8", errors="replace")
    if "__main__" in text:
        return [sys.executable, "-m", f"tests.{name}"]
    return [sys.executable, str(ROOT / "verify.py"), "--run-test-module", f"tests.{name}"]


def run_one_test(name: str) -> dict:
    """Run one test module in a subprocess; return a check row for the table.

    Course topic: (none - project plumbing)
    """
    try:
        proc = subprocess.run(_test_command(name), cwd=ROOT, env=child_env(),
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, timeout=TEST_TIMEOUT_S)
        code, output = proc.returncode, proc.stdout or ""
    except subprocess.TimeoutExpired as exc:
        code, output = 124, (exc.stdout or "") + f"\n(timed out after {TEST_TIMEOUT_S} s)"
    lines = output.splitlines()
    n_pass = sum(1 for ln in lines if ln.startswith("PASS "))
    n_fail = sum(1 for ln in lines if ln.startswith("FAIL "))
    ok = (code == 0 and n_fail == 0)
    return {"module": "tests", "name": f"tests.{name}", "pass": ok,
            "value": f"exit {code}, {n_pass} PASS / {n_fail} FAIL lines", "limit": "exit 0",
            "output": output}


def run_tests(jobs: int = 2, echo: bool = True) -> list[dict]:
    """Run every existing test module (up to `jobs` at a time).

    Course topic: (none - project plumbing)
    """
    from concurrent.futures import ThreadPoolExecutor
    names = find_test_modules()
    if not names:
        if echo:
            print(paint("SKIP", "yellow") + "  tests: no tests/test_*.py found yet")
        return []
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        rows = list(pool.map(run_one_test, names))
    if echo:
        for row in rows:
            tag = paint("OK  ", "green") if row["pass"] else paint("FAIL", "red")
            print(f"{tag}  {row['name']}: {row['value']}")
            if not row["pass"]:                      # show why: the FAIL lines, else the end of the output
                lines = row["output"].splitlines()
                shown = [ln for ln in lines if ln.startswith("FAIL ")] or lines[-25:]
                for ln in shown:
                    print("      | " + ln)
    return rows


def _run_test_module_inline(dotted: str) -> int:
    """Mini test runner (fallback when a test file has no __main__ block).

    Imports the module, runs unittest.TestCase classes and plain `test_*()`
    functions that take no arguments. Prints one PASS/FAIL line per test.
    Course topic: (none - project plumbing)
    """
    import importlib
    import traceback
    import unittest
    sys.path.insert(0, str(ROOT))
    try:
        mod = importlib.import_module(dotted)
    except Exception:                                # import error = failure, with the reason
        traceback.print_exc()
        print(f"FAIL {dotted}: import failed")
        return 1
    failures = 0
    suite = unittest.defaultTestLoader.loadTestsFromModule(mod)
    if suite.countTestCases():
        result = unittest.TextTestRunner(stream=sys.stdout, verbosity=1).run(suite)
        failures += len(result.failures) + len(result.errors)
    for attr in sorted(vars(mod)):
        func = getattr(mod, attr)
        if not (attr.startswith("test_") and callable(func)) or isinstance(func, type):
            continue
        try:
            func()
            print(f"PASS {dotted}: {attr}")
        except TypeError as exc:                     # needs fixtures/arguments: cannot run here
            if "argument" in str(exc):
                print(f"SKIP {dotted}: {attr} needs arguments ({exc})")
            else:
                failures += 1
                traceback.print_exc()
                print(f"FAIL {dotted}: {attr}")
        except Exception:
            failures += 1
            traceback.print_exc()
            print(f"FAIL {dotted}: {attr}")
    return 1 if failures else 0


# --------------------------------------------------------------------------
# Reading outputs/data/*.json
# --------------------------------------------------------------------------
def fmt(value) -> str:
    """Short, stable text for a check value or limit.

    Course topic: (none - project plumbing)
    """
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        return f"{value:.6g}"
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(fmt(v) for v in value) + "]"
    return str(value)


def collect(data_dir: Path) -> tuple[list[dict], dict, list[str]]:
    """Read every JSON; return (check rows, verdicts, notes).

    A JSON that cannot be parsed, or a check entry without a boolean "pass",
    becomes a FAIL row (a broken file must never look like a success).
    A JSON without a "checks" list is only a note (it may be a data cache).
    Course topic: (none - project plumbing)
    """
    rows: list[dict] = []
    verdicts: dict = {}
    notes: list[str] = []
    for path in sorted(data_dir.glob("*.json")):
        module = path.stem
        try:
            with open(path, encoding="utf-8") as fh:
                doc = json.load(fh)
        except (OSError, ValueError) as exc:
            rows.append({"module": module, "name": "JSON file is readable", "pass": False,
                         "value": type(exc).__name__, "limit": "valid JSON"})
            continue
        checks = doc.get("checks") if isinstance(doc, dict) else None
        if not isinstance(checks, list):
            notes.append(f"{path.name}: no \"checks\" list (not counted)")
        else:
            if not checks:
                notes.append(f"{path.name}: \"checks\" list is empty")
            for item in checks:
                item = item if isinstance(item, dict) else {"name": fmt(item)}
                rows.append({"module": module, "name": str(item.get("name", "(unnamed check)")),
                             "pass": item.get("pass") is True,       # anything else counts as FAIL
                             "value": item.get("value", ""), "limit": item.get("limit", "")})
        found = doc.get("verdicts") if isinstance(doc, dict) else None
        if isinstance(found, dict):
            for hyp, body in found.items():
                body = body if isinstance(body, dict) else {"verdict": fmt(body)}
                verdicts[hyp] = {"module": module, "verdict": str(body.get("verdict", "?")),
                                 "rule": str(body.get("rule", "(no rule stored)"))}
    return rows, verdicts, notes


# --------------------------------------------------------------------------
# Building the report text
# --------------------------------------------------------------------------
def _clip(text: str, width: int) -> str:
    """Cut long text so the table stays aligned."""
    return text if len(text) <= width else text[: width - 3] + "..."


def build_report(rows: list[dict], verdicts: dict, notes: list[str],
                 colour: bool, strict_missing: list[str]) -> tuple[str, int, int]:
    """Return (report text, number passed, number failed).

    Course topic: (none - project plumbing)
    """
    caps = {"module": 20, "name": 62, "value": 30, "limit": 30}
    cells = [{"module": _clip(r["module"], caps["module"]), "name": _clip(r["name"], caps["name"]),
              "value": _clip(fmt(r["value"]), caps["value"]),
              "limit": _clip(fmt(r["limit"]), caps["limit"]), "pass": r["pass"]} for r in rows]
    heads = {"module": "module", "name": "check", "value": "value", "limit": "limit"}
    width = {k: max([len(heads[k])] + [len(c[k]) for c in cells]) for k in heads}

    def line(c: dict, status: str) -> str:
        return " | ".join([c["module"].ljust(width["module"]), c["name"].ljust(width["name"]),
                           c["value"].rjust(width["value"]), c["limit"].rjust(width["limit"]),
                           status])

    out = ["VERIFICATION TABLE  (source: outputs/data/*.json + test modules)", ""]
    header = line(heads, "result")
    out += [paint(header, "bold", colour), "-" * len(header)]
    for c in cells:
        status = paint("PASS", "green", colour) if c["pass"] else paint("FAIL", "red", colour)
        out.append(line(c, status))
    if not cells:
        out.append("(no checks found)")
    n_pass = sum(1 for c in cells if c["pass"])
    n_fail = len(cells) - n_pass + len(strict_missing)
    out += ["-" * len(header)]
    summary = f"{n_pass} passed, {n_fail} failed"
    out.append(paint(summary, "red" if n_fail else "green", colour))

    for name in strict_missing:
        out.append(paint("FAIL", "red", colour) + f" (--strict) no outputs/data/{name}.json")

    out += ["", "HYPOTHESIS VERDICTS  (written by code from numbers; rule stored in the JSON)"]
    shade = {"supported": "green", "not supported": "red", "mixed": "yellow"}
    for hyp in list(HYPOTHESES) + sorted(h for h in verdicts if h not in HYPOTHESES):
        if hyp in verdicts:
            v = verdicts[hyp]
            out.append(f"  {hyp}: " + paint(v["verdict"], shade.get(v["verdict"], "cyan"), colour)
                       + f"   [{v['module']}]")
            out.append(f"      rule: {v['rule']}")
        else:
            out.append(f"  {hyp}: " + paint("(no verdict yet - numbers not produced)", "dim", colour))
    if notes:
        out += ["", "NOTES"] + [f"  {n}" for n in notes]
    return "\n".join(out) + "\n", n_pass, n_fail


def verify(run_test_modules: bool = True, jobs: int = 2, strict: bool = False,
           extra_rows: list[dict] | None = None) -> int:
    """Whole verification: tests (optional) + table + report file. Returns the exit code.

    `extra_rows` lets the CLI add its own rows (for example pipeline stages that crashed).
    Course topic: (none - project plumbing)
    """
    out = out_dir()
    data = out / "data"
    rows: list[dict] = list(extra_rows or [])
    if run_test_modules:
        rows += run_tests(jobs=jobs)
        print()
    json_files = sorted(data.glob("*.json")) if data.is_dir() else []
    json_rows, verdicts, notes = collect(data) if json_files else ([], {}, [])
    rows += json_rows

    present = {p.stem for p in json_files}
    missing = [m for m in PIPELINE if m not in present]
    if missing:
        notes.append("pipeline modules with no JSON yet: " + ", ".join(missing))
    strict_missing = missing if strict else []

    text, _, n_fail = build_report(rows, verdicts, notes, use_colour(), strict_missing)
    plain, _, _ = build_report(rows, verdicts, notes, False, strict_missing)
    print(text, end="")

    code = 1 if n_fail else 0
    if not json_files:
        code = 1
        msg = ("No outputs/data/*.json found: the study has not been run yet.\n"
               "Run it first:  ./run.sh full --fast   (quick)   or   ./run.sh full   (complete)")
        print(paint(msg, "red"))
        plain += msg + "\n"
    try:                                             # report file: plain text, no timestamps
        out.mkdir(parents=True, exist_ok=True)
        (out / "verify_report.txt").write_text(plain, encoding="utf-8")
    except OSError as exc:
        print(paint(f"(could not write verify_report.txt: {exc})", "yellow"))
    return code


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point.

    Course topic: (none - project plumbing)
    """
    import argparse
    parser = argparse.ArgumentParser(description="Print the PASS/FAIL table of the study.")
    parser.add_argument("--no-tests", action="store_true", help="do not run tests/test_*.py")
    parser.add_argument("--strict", action="store_true",
                        help="fail when a pipeline module has produced no JSON")
    parser.add_argument("--jobs", type=int, default=2, help="test modules run at once (default 2)")
    parser.add_argument("--run-test-module", metavar="MODULE", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.run_test_module:                         # internal: fallback mini test runner
        return _run_test_module_inline(args.run_test_module)
    return verify(run_test_modules=not args.no_tests, jobs=args.jobs, strict=args.strict)


if __name__ == "__main__":
    sys.exit(main())
