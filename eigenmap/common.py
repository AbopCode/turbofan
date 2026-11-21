"""Shared helpers for the study modules: JSON output, PASS/FAIL lines, eigen-summaries.

Course topic: (none - project plumbing)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from . import config
from . import mathcore as mc


class Checks:
    """Collects PASS/FAIL checks, prints them as CONTRACT section 8 lines, stores them in the JSON."""

    def __init__(self, module: str):
        self.module = module
        self.items: list[dict] = []

    def add(self, name: str, ok: bool, value, limit) -> bool:
        ok = bool(ok)
        self.items.append({"name": name, "pass": ok, "value": _plain(value), "limit": _plain(limit)})
        vs = "%.3g" % value if isinstance(value, (float, np.floating)) else str(value)
        ls = "%.3g" % limit if isinstance(limit, (float, np.floating)) else str(limit)
        print("%s %s: %s (%s vs %s)" % ("PASS" if ok else "FAIL", self.module, name, vs, ls), flush=True)
        return ok

    @property
    def failed(self) -> int:
        return sum(1 for c in self.items if not c["pass"])


def _plain(v):
    """Convert numpy / complex values to JSON-friendly plain Python (complex -> [re, im])."""
    if isinstance(v, dict):
        return {str(k): _plain(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_plain(x) for x in v]
    if isinstance(v, np.ndarray):
        return _plain(v.tolist())
    if isinstance(v, (complex, np.complexfloating)):
        return [float(np.real(v)), float(np.imag(v))]
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    if isinstance(v, float):
        return float(repr(v)) if v == v else None
    return v


def write_json(module: str, payload: dict, checks: Checks) -> Path:
    """Write outputs/data/<module>.json deterministically (sorted keys, no timestamps)."""
    config.DATA.mkdir(parents=True, exist_ok=True)
    payload = dict(payload)
    payload["checks"] = checks.items
    path = config.DATA / ("%s.json" % module)
    path.write_text(json.dumps(_plain(payload), indent=1, sort_keys=True) + "\n")
    return path


def finish(checks: Checks) -> None:
    """Exit with code 1 if any check failed."""
    print("%s: %d checks, %d failed" % (checks.module, len(checks.items), checks.failed), flush=True)
    sys.exit(1 if checks.failed else 0)


def eig_summary(A: np.ndarray) -> dict:
    """Eigen-analysis of a 2x2 state matrix with the project's OWN formulas (mathcore).

    Course topic: Eigenvalues, eigenvectors, stability
    """
    lam, info = mc.eig2x2(A)
    V = mc.eigvec2x2(A, lam)
    tau_s, tau_f = mc.time_constants(lam)
    return {"lam": [complex(l) for l in lam], "tr": float(info["tr"]), "det": float(info["det"]),
            "disc": float(info["disc"]), "tau_slow": float(tau_s), "tau_fast": float(tau_f),
            "V": V, "kappaV": float(mc.kappa2(V)), "stable": bool(mc.is_stable2x2(A)),
            "complex_pair": bool(info["disc"] < 0)}
