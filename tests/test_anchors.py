"""Re-assert the published-example anchor: turboshaft JetCat SPT5 (illustration).

Run:  cd $P && python -m tests.test_anchors
Runs `eigenmap.anchors` with the FAST configuration into a temporary folder (the
real `outputs/data/anchors.json` is not touched), then re-checks the stored
numbers independently with numpy / scipy, the JSON conventions of CONTRACT.md
sections 8 and 11, and byte-for-byte determinism (twice in this process and once
in a separate Python process). Prints `PASS ...` / `FAIL ...`; exit code 1 on FAIL.

Course topic: (none - test harness)
"""
from __future__ import annotations

import contextlib
import hashlib
import io
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import scipy.linalg

from eigenmap import anchors
from eigenmap.config import ROOT, StudyConfig

RESULTS: list[bool] = []


def check(ok: bool, claim: str, value, limit) -> None:
    """Record and print one PASS/FAIL line.

    Course topic: (none - test harness)
    """
    ok = bool(ok)
    v = f"{value:.6g}" if isinstance(value, (float, np.floating)) else str(value)
    lim = f"{limit:.6g}" if isinstance(limit, (float, np.floating)) else str(limit)
    print(f"{'PASS' if ok else 'FAIL'} test_anchors: {claim} ({v} vs {lim})")
    RESULTS.append(ok)


def sha256(path: Path) -> str:
    """SHA-256 of a file.

    Course topic: (none - test harness)
    """
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cplx(pairs) -> np.ndarray:
    """JSON [[re, im], ...] -> complex array.

    Course topic: (none - test harness)
    """
    return np.array([complex(p[0], p[1]) for p in pairs])


def test_run_and_determinism(tmp: Path) -> dict:
    """Run the module (fast config), check exit status, checks list and byte-identical output.

    Course topic: (none - test harness)
    """
    cfg = StudyConfig(fast=True)
    f1, f2 = tmp / "run1" / "anchors.json", tmp / "run2" / "anchors.json"
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        res = anchors.run(cfg, out_path=f1)
        anchors.run(cfg, out_path=f2)
    lines = [ln for ln in log.getvalue().splitlines() if ln.startswith(("PASS ", "FAIL "))]
    n_checks = len(res["checks"])
    n_fail = sum(1 for c in res["checks"] if not c["pass"])
    check(n_fail == 0 and n_checks > 0, f"anchors.run(fast): all {n_checks} checks PASS", n_fail, 0)
    check(len(lines) == 2 * n_checks and all(ln.split(" ", 2)[1] == "anchors:" for ln in lines),
          "one 'PASS anchors: ...' console line per check", len(lines) // 2, n_checks)
    check(all(anchors.LABEL in ln for ln in lines), f"every console line says '{anchors.LABEL}'", True, True)
    check(f1.read_bytes() == f2.read_bytes(), "two runs in one process give byte-identical anchors.json",
          sha256(f1)[:16], sha256(f2)[:16])
    env = dict(os.environ, EIGENMAP_OUT=str(tmp / "run3"), PYTHONHASHSEED="12345")
    proc = subprocess.run([sys.executable, "-m", "eigenmap.anchors", "--fast"], cwd=str(ROOT), env=env,
                          capture_output=True, text=True)
    f3 = tmp / "run3" / "data" / "anchors.json"
    check(proc.returncode == 0 and f3.exists(), "`python -m eigenmap.anchors --fast` exits 0 and writes "
          "data/anchors.json", proc.returncode, 0)
    check(f3.exists() and f3.read_bytes() == f1.read_bytes(),
          "a separate Python process gives the byte-identical file", sha256(f3)[:16] if f3.exists() else "missing",
          sha256(f1)[:16])
    text = f1.read_text(encoding="utf-8")
    data = json.loads(text)
    check(data == res, "run() returns exactly what it wrote", True, True)
    check(text == json.dumps(data, indent=1, sort_keys=True) + "\n", "JSON is written with indent=1, sort_keys=True",
          True, True)
    bad = [w for w in ("/home", "\\\\", "timestamp", "hostname", "NaN", "Infinity") if w in text]
    check(not bad, "JSON holds no paths, time stamps, machine info or non-finite numbers", str(bad), "[]")
    ok = all(set(c) == {"name", "pass", "value", "limit"} and isinstance(c["name"], str)
             and isinstance(c["pass"], bool) and isinstance(c["value"], (int, float, str, bool))
             and isinstance(c["limit"], (int, float, str, bool)) for c in data["checks"])
    check(ok, "checks list items are {name, pass, value, limit} (CONTRACT section 11)", ok, True)
    return data


def test_labels(data: dict) -> None:
    """The example is always labelled as the turboshaft illustration, never as the turbofan.

    Course topic: (none - test harness)
    """
    label = anchors.LABEL
    labelled = [data["label"], data["engine"], data["disclaimer"], data["gauss_seidel"]["label"],
                data["flip"]["label"], data["monte_carlo"]["label"]]
    labelled += [data["points"][k]["label"] for k in data["order"]]
    labelled += [c["name"] for c in data["checks"]]
    check(label == "turboshaft JetCat SPT5 (illustration)" and all(label in s for s in labelled),
          f"all {len(labelled)} labels / check names contain '{label}'", True, True)
    text = json.dumps(data)
    n_agtf = text.count("AGTF30")
    n_neg = text.count("NOT the AGTF30 turbofan")
    check(n_agtf == n_neg and "turbofan" not in text.replace("NOT the AGTF30 turbofan", ""),
          "the word 'turbofan' appears only in the disclaimer 'NOT the AGTF30 turbofan'", n_agtf - n_neg, 0)
    check("arXiv:1403.3738" in data["source"] and "(57)-(59)" in data["source"], "source citation stored",
          True, True)


def test_numbers(data: dict) -> None:
    """Recompute the stored anchor numbers with numpy / scipy only and compare.

    Course topic: Eigenvalues
    """
    brief = {"idle": ((-0.346, -0.374), (2.89, 2.67), 18.9), "mid": ((-0.599, -0.891), (1.67, 1.12), 2.39),
             "cruise": ((-1.013, -1.787), (0.99, 0.56), 1.84)}
    printed = {"idle": [[-0.38, -0.0008], [0.26, -0.34]], "mid": [[-0.85, 0.032], [0.32, -0.64]],
               "cruise": [[-1.7, 0.1], [0.6, -1.1]]}
    printed_b = {"idle": [[0.7, 0.0], [0.1, -0.0024]], "mid": [[1.0, 0.0], [0.17, -0.011]],
                 "cruise": [[1.2, 0.0], [0.3, -0.023]]}
    check(data["order"] == ["idle", "mid", "cruise"]
          and all(data["points"][k]["A"] == printed[k] and data["points"][k]["B"] == printed_b[k] for k in printed),
          "stored A and B equal the printed matrices of the brief", True, True)
    for key, (lam_b, tau_b, kap_b) in brief.items():
        p = data["points"][key]
        A, B = np.array(p["A"]), np.array(p["B"])
        w, V = np.linalg.eig(A)
        w = np.sort(w.real)[::-1]                                   # slow (least negative) first
        own = cplx(p["eig"]["own"])
        e1 = float(np.max(np.abs(own - w)))
        e2 = float(np.max(np.abs(w - np.array(lam_b))))
        check(e1 <= 1e-10 and e2 <= 2e-3, f"{key}: stored eigenvalues = numpy ({w[0]:.5f}, {w[1]:.5f}); "
              f"within 2e-3 of the brief", max(e1, e2), 2e-3)
        tau = -1.0 / w
        e1 = max(abs(p["tau"]["slow"] - tau[0]), abs(p["tau"]["fast"] - tau[1]))
        rounded_ok = round(tau[0], 2) == tau_b[0] and round(tau[1], 2) == tau_b[1]
        check(e1 <= 1e-10 and rounded_ok, f"{key}: stored tau = -1/lambda ({tau[0]:.4f}, {tau[1]:.4f}); rounds "
              f"to the brief's {tau_b[0]} / {tau_b[1]}", e1, 1e-10)
        kap = float(np.linalg.cond(V, 2))
        e1 = abs(p["kappaV"]["own"] - kap) / kap
        e2 = abs(kap - kap_b) / kap_b
        check(e1 <= 1e-10 and e2 <= 0.01, f"{key}: stored kappa(V) = numpy cond ({kap:.4f}); within 1% of the "
              f"brief's {kap_b}", max(e1, e2), 0.01)
        t = np.array(p["step_fuel"]["t"])
        dx = np.column_stack([p["step_fuel"]["dx1"], p["step_fuel"]["dx2"]])
        aug = np.zeros((3, 3))
        aug[:2, :2], aug[:2, 2] = A, B[:, 0]
        ref = np.array([scipy.linalg.expm(aug * tk)[:2, 2] for tk in t])
        ss = -np.linalg.solve(A, B[:, 0])
        e1 = float(np.max(np.abs(dx - ref)))
        e2 = float(np.max(np.abs(np.array(p["steady_gain"]["fuel flow"]) - ss)))
        ok = t[0] == 0.0 and len(t) == 201 and bool(np.all(np.diff(t) > 0)) and bool(np.all(dx[0] == 0.0))
        check(e1 <= 1e-10 and e2 <= 1e-10 and ok, f"{key}: stored fuel-step response (201 points) and steady "
              "gain = scipy/numpy", max(e1, e2), 1e-10)
    tau_only = sorted((r["point"], r["quantity"]) for r in data["matched_to_printed_rounding_only"])
    check(tau_only == [("cruise", "tau_slow"), ("idle", "tau_fast"), ("mid", "tau_fast")],
          "values matched to printed rounding only are listed (3 tau values, |diff| between 2e-3 and 5e-3)",
          str(tau_only), "idle tau_fast, mid tau_fast, cruise tau_slow")

    # Gauss-Seidel table recomputed with a plain loop
    gs = data["gauss_seidel"]
    A, b = np.array(gs["A"]), np.array(gs["b"])
    x = np.zeros(2)
    brief_gs = [(0.70588, 0.65775), (0.74457, 0.67886), (0.74582, 0.67954), (0.74586, 0.67956)]
    worst_tab = 0.0
    ok5 = True
    for k in range(1, 5):
        x[0] = (b[0] - A[0, 1] * x[1]) / A[0, 0]
        x[1] = (b[1] - A[1, 0] * x[0]) / A[1, 1]
        row = gs["table"][k]
        worst_tab = max(worst_tab, abs(row["x1"] - x[0]), abs(row["x2"] - x[1]))
        ok5 &= (round(x[0], 5), round(x[1], 5)) == brief_gs[k - 1] and row["k"] == k
    exact = np.linalg.solve(A, b)
    check(b.tolist() == [-1.2, -0.3] and worst_tab <= 1e-14 and ok5,
          "Gauss-Seidel iterates 1-4 recomputed; rounded to 5 d.p. they equal the brief's table", worst_tab, 1e-14)
    check((round(exact[0], 5), round(exact[1], 5)) == (0.74586, 0.67956)
          and np.max(np.abs(np.array(gs["exact_own"]) - exact)) <= 1e-12,
          "Gauss-Seidel exact solution 0.74586 / 0.67956", f"{exact[0]:.5f} / {exact[1]:.5f}", "0.74586 / 0.67956")
    rho = abs(A[0, 1] * A[1, 0] / (A[0, 0] * A[1, 1]))
    check(abs(gs["spectral_radius"] - rho) <= 1e-15 and round(rho, 3) == 0.032,
          "Gauss-Seidel spectral radius about 0.032", rho, 0.032)
    errs = [r["err_inf"] for r in gs["table"]]
    check(all(errs[k + 1] < errs[k] for k in range(len(errs) - 2)) and gs["converged"],
          "Gauss-Seidel table errors decrease monotonically; converged", errs[-1], 1e-10)

    # flip
    fl = data["flip"]
    w = np.linalg.eigvals(np.array(fl["A_changed"]))
    check(np.all(np.linalg.eigvals(np.array(fl["A_printed"])).imag == 0.0) and abs(w[0].imag) > 0.0
          and abs(w[0].real + 0.36) <= 2e-3 and round(abs(w[0].imag), 3) == 0.008
          and fl["A_changed"][0][1] == -0.0018,
          "idle flip: real pair -> complex pair about -0.36 +/- 0.008i when a12 = -0.0018",
          f"{w[0].real:.4f} +/- {abs(w[0].imag):.5f}i", "-0.36 +/- 0.008i")

    # Monte Carlo
    cfg = StudyConfig(fast=True)
    mcd = data["monte_carlo"]
    seeds_ok = [r["seed"] for r in mcd["per_seed"]] == list(cfg.mc_seeds) and mcd["n"] == cfg.n_mc \
        and mcd["E_norm"] == cfg.E_norm
    check(seeds_ok, "Monte Carlo uses cfg.mc_seeds, cfg.n_mc, cfg.E_norm", str([r["seed"] for r in mcd["per_seed"]]),
          str(list(cfg.mc_seeds)))
    for r in mcd["per_seed"]:
        ok = (abs(r["worst_shift"]["idle"] - 0.016) <= 0.15 * 0.016
              and abs(r["worst_shift"]["cruise"] - 0.0012) <= 0.15 * 0.0012
              and 10.0 <= r["ratio_idle_over_cruise"] <= 16.0
              and all(r["worst_shift"][k] <= mcd["bauer_fike_bound"][k] for k in data["order"]))
        check(ok, f"seed {r['seed']}: worst shift idle {r['worst_shift']['idle']:.5f} (~0.016), cruise "
              f"{r['worst_shift']['cruise']:.5f} (~0.0012), ratio in [10, 16], below the Bauer-Fike bounds",
              r["ratio_idle_over_cruise"], "10 .. 16")
    # independent numpy re-computation of the first seed (same generator recipe, numpy eigenvalues)
    r = mcd["per_seed"][0]
    rng = np.random.default_rng(r["seed"])
    E = rng.standard_normal((mcd["n"], 2, 2))
    for _ in range(2):
        E = E * (mcd["E_norm"] / np.linalg.norm(E, 2, axis=(1, 2)))[:, None, None]
    worst = 0.0
    for k in data["order"]:
        A = np.array(data["points"][k]["A"])
        w, mu = np.linalg.eigvals(A), np.linalg.eigvals(A[None] + E)
        d = np.abs(mu[:, :, None] - w[None, None, :]).min(axis=2).max(axis=1)
        worst = max(worst, abs(float(d.max()) - r["worst_shift"][k]))
    check(worst <= 1e-10, "seed 1 worst shifts recomputed with numpy only (numpy norms and eigenvalues)", worst, 1e-10)
    kb = [mcd["bauer_fike_bound"][k] for k in data["order"]]
    ref = [np.linalg.cond(np.linalg.eig(np.array(data["points"][k]["A"]))[1], 2) * mcd["E_norm"]
           for k in data["order"]]
    check(np.allclose(kb, ref, rtol=1e-10, atol=0.0), "stored Bauer-Fike bounds = numpy cond(V) * ||E||_2",
          str([round(v, 6) for v in kb]), "about [0.018875, 0.002393, 0.001836]")


def main() -> int:
    """Run all anchor tests; return the process exit code.

    Course topic: (none - test harness)
    """
    with tempfile.TemporaryDirectory(prefix="eigenmap_anchors_") as tmp:
        data = test_run_and_determinism(Path(tmp))
    test_labels(data)
    test_numbers(data)
    n_fail = sum(1 for ok in RESULTS if not ok)
    print(f"test_anchors: {len(RESULTS) - n_fail} PASS, {n_fail} FAIL")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
