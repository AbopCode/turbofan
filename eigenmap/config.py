"""Study configuration: every tunable number of the study lives here.

Owned by the orchestrator. Modules READ this; they never hard-code these values.
Course topic: (none - project plumbing)
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent          # project root
VENDOR = ROOT / "vendor"
OUT = Path(os.environ.get("EIGENMAP_OUT", ROOT / "outputs"))   # override for clean-room / determinism runs
DATA = OUT / "data"
FIGS = OUT / "figures"
TABLES = OUT / "tables"

BASE_SEED = 20261114
MC_SEEDS = (20261114, 20261115, 20261116)


@dataclass(frozen=True)
class StudyConfig:
    """All settings of one study run. `fast=True` is the coarse 'skeleton' run.

    Course topic: (none - project plumbing)
    """
    fast: bool = False

    # ---- solver tolerances (tol = inner flow-matching Newton tolerance) ----
    tol_tight: float = 1e-13          # used for steady states and reference Jacobians
    tol_default: float = 1e-12
    fd_tols: tuple = (1e-13, 1e-9, 1e-6, 1e-4)   # Layer 2 inner-solver tolerances

    # ---- seeds ----
    seed: int = BASE_SEED
    mc_seeds: tuple = MC_SEEDS

    # ---- Bauer-Fike (Layer 3) ----
    E_norm: float = 1e-3              # ||E||_2 of the random perturbations
    E_scales: tuple = (1e-5, 1e-4, 1e-3, 1e-2, 1e-1)

    # ---- degradation (Layer 4) ----
    deg_max: float = 0.05             # 0 .. 5 %

    # ---- validity study ----
    step_fracs: tuple = (0.01, 0.03, 0.10, 0.25)

    @property
    def n_sls(self) -> int:
        """Number of sea-level-static power settings (idle .. max)."""
        return 4 if self.fast else 10

    @property
    def n_cruise(self) -> int:
        """Number of power settings in the 35,000 ft / Mach 0.8 sweep."""
        return 3 if self.fast else 8

    @property
    def n_mc(self) -> int:
        """Monte Carlo samples per operating point."""
        return 500 if self.fast else 20000

    @property
    def h_grid(self) -> tuple:
        """Relative finite-difference steps h, 1e-12 .. 1e-1 (log spaced)."""
        per_decade = 1 if self.fast else 4
        n = 11 * per_decade + 1
        return tuple(10.0 ** (-12.0 + 11.0 * i / (n - 1)) for i in range(n))

    @property
    def deg_grid(self) -> tuple:
        """Degradation levels as fractions (0 .. deg_max)."""
        n = 3 if self.fast else 6
        return tuple(self.deg_max * i / (n - 1) for i in range(n))

    def ensure_dirs(self) -> None:
        """Create the output folders."""
        for p in (OUT, DATA, FIGS, TABLES):
            p.mkdir(parents=True, exist_ok=True)


def from_argv(argv: list[str] | None = None) -> StudyConfig:
    """Build a StudyConfig from command-line flags (`--fast`)."""
    import sys
    argv = sys.argv[1:] if argv is None else argv
    return StudyConfig(fast=("--fast" in argv))
