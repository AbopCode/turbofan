"""LAYER 3 - idle fragility: how far can a small error in A move the eigenvalues?  (Bauer-Fike)

If A = V D V^-1 and E is any perturbation, every eigenvalue mu of A + E satisfies
    min_j |mu - lam_j| <= kappa(V) ||E||_2,       kappa(V) = ||V||_2 ||V^-1||_2.
So the eigenvector condition number kappa(V) amplifies errors in the Jacobian (e.g. from finite
differences, Layer 2). Monte Carlo: 20,000 random E with ||E||_2 = 1e-3 at every operating point,
three seeds, plus a sweep of ||E||_2.

The matrix studied is the NORMALIZED state matrix A_n = D^-1 A D, D = diag(x0) (states as fractions
of the operating speed), so that ||E||_2 = 1e-3 means the same relative error at every point.

Tests hypothesis H2: "the slowest operating point is also the most fragile".

Run: python -m eigenmap.layer3_bauerfike [--fast]
Course topic: Eigenvalues; Matrix norms and conditioning (Bauer-Fike theorem)
"""
from __future__ import annotations

import numpy as np

from . import config
from . import mathcore as mc
from .common import Checks, finish, write_json
from .layer1_opmap import point_record
from .opoints import standard_points


def monte_carlo(An: np.ndarray, n: int, norm: float, seed: int) -> dict:
    """Eigenvalue shifts of A_n + E for n random E with ||E||_2 = norm."""
    rng = np.random.default_rng(seed)
    Es = mc.random_perturbations(rng, n, norm)
    sh = mc.eig_shift_batch(An, Es)
    lam, _ = mc.eig2x2(An)
    V = mc.eigvec2x2(An, lam)
    bound = mc.kappa2(V) * norm
    pert = np.linalg.eigvals(An[None, :, :] + Es)          # numpy cross-check on a subset
    k = min(n, 200)
    sh_np = np.array([max(min(abs(m - l) for l in np.linalg.eigvals(An)) for m in pert[i]) for i in range(k)])
    # real-to-complex flips: perturbations that change the type of the eigenvalue pair
    disc0 = (An[0, 0] - An[1, 1]) ** 2 + 4 * An[0, 1] * An[1, 0]
    M = An[None] + Es
    disc = (M[:, 0, 0] - M[:, 1, 1]) ** 2 + 4 * M[:, 0, 1] * M[:, 1, 0]
    return {"worst": float(np.max(sh)), "median": float(np.median(sh)), "p99": float(np.quantile(sh, 0.99)),
            "bound": float(bound), "worst_over_bound": float(np.max(sh) / bound),
            "own_vs_numpy": float(np.max(np.abs(sh[:k] - sh_np))), "flip_fraction": float(np.mean(np.sign(disc) != np.sign(disc0)))}


def h2_verdict(recs: list) -> dict:
    """H2 'the slowest point is also the most fragile', decided by an explicit rule for every seed."""
    per_seed = {}
    for seed in recs[0]["mc"]:
        slow = max(recs, key=lambda r: r["tau_slow"])
        frag = max(recs, key=lambda r: r["mc"][seed]["worst"])
        rank = sorted(recs, key=lambda r: -r["mc"][seed]["worst"]).index(slow) + 1
        per_seed[seed] = {"slowest_point": slow["id"], "most_fragile_point": frag["id"], "same": slow["id"] == frag["id"],
                          "fragility_rank_of_slowest": rank, "of_points": len(recs),
                          "worst_shift_slowest": slow["mc"][seed]["worst"], "worst_shift_max": frag["mc"][seed]["worst"]}
    same = [v["same"] for v in per_seed.values()]
    verdict = "supported" if all(same) else ("not supported" if not any(same) else "mixed")
    return {"verdict": verdict, "per_seed": per_seed,
            "rule": "over all SLS and cruise operating points (normalized A_n, ||E||_2 = 1e-3): supported if the point with the "
                    "largest tau_slow is also the point with the largest worst-case eigenvalue shift for every Monte Carlo seed; "
                    "not supported if it is not for any seed; otherwise mixed"}


def run(cfg: config.StudyConfig) -> Checks:
    ck = Checks("layer3_bauerfike")
    engine, pts = standard_points(cfg)
    recs = []
    for p in pts:
        if p.group == "ANCHOR":
            continue
        r = point_record(engine, p)
        An = np.array(r["A_n"])
        rec = {"id": r["id"], "label": r["label"], "group": r["group"], "power": r["power"], "tau_slow": r["tau_slow"],
               "tau_fast": r["tau_fast"], "kappaV": r["kappaV"], "kappaV_raw": r["kappaV_raw"], "lam": r["lam"],
               "complex_pair": r["complex_pair"], "A_n": An, "mc": {}}
        for seed in cfg.mc_seeds:
            m = monte_carlo(An, cfg.n_mc, cfg.E_norm, seed)
            rec["mc"][str(seed)] = m
            ck.add("%s seed %d Bauer-Fike bound holds (worst shift / bound)" % (r["id"], seed), m["worst_over_bound"] <= 1.0 + 1e-12,
                   m["worst_over_bound"], 1.0)
            ck.add("%s seed %d own batch eigenvalues vs numpy (abs)" % (r["id"], seed), m["own_vs_numpy"] < 1e-10, m["own_vs_numpy"], 1e-10)
        rec["scale_sweep"] = {"%g" % s: monte_carlo(An, max(500, cfg.n_mc // 10), s, cfg.seed) for s in cfg.E_scales}
        recs.append(rec)
    seeds = [str(s) for s in cfg.mc_seeds]
    spread = {r["id"]: {"worst_min": min(r["mc"][s]["worst"] for s in seeds), "worst_max": max(r["mc"][s]["worst"] for s in seeds),
                        "rel_spread": (max(r["mc"][s]["worst"] for s in seeds) - min(r["mc"][s]["worst"] for s in seeds))
                        / np.mean([r["mc"][s]["worst"] for s in seeds])} for r in recs}
    idle = [r for r in recs if r["group"] == "SLS"][0]
    cruise = [r for r in recs if r["group"] == "CRUISE"][-1]
    ratio = {s: idle["mc"][s]["worst"] / cruise["mc"][s]["worst"] for s in seeds}
    h2 = h2_verdict(recs)
    seed0 = str(cfg.mc_seeds[0])
    # sensitivity: the OTHER idle steady state (layer1_opmap found it; same Nc_fan, VAFN, VBV, different branch)
    import json
    l1p = config.DATA / "layer1_opmap.json"
    other = {}
    if l1p.exists():
        for pid, b in json.loads(l1p.read_text()).get("idle_branches", {}).items():
            if b.get("other_found"):
                o = b["other"]
                x0 = np.array(o["x0"])
                An = np.diag(1.0 / x0) @ np.array(o["A"]) @ np.diag(x0)
                m = monte_carlo(An, cfg.n_mc, cfg.E_norm, cfg.seed)
                slowest_all = max([r["tau_slow"] for r in recs if r["id"] != pid] + [o["tau_slow"]])
                worst_all = max([r["mc"][seed0]["worst"] for r in recs if r["id"] != pid] + [m["worst"]])
                other[pid] = {"tau_slow": o["tau_slow"], "kappaV": m["bound"] / cfg.E_norm, "worst": m["worst"], "bound": m["bound"],
                              "would_be_slowest_overall": o["tau_slow"] >= slowest_all, "would_be_most_fragile_overall": m["worst"] >= worst_all,
                              "H2_same_point_if_used": (o["tau_slow"] >= slowest_all) == (m["worst"] >= worst_all) and o["tau_slow"] >= slowest_all}
    h2["idle_branch_sensitivity"] = other
    h2["within_group"] = {g: {"slowest": max([r for r in recs if r["group"] == g], key=lambda r: r["tau_slow"])["id"],
                              "most_fragile": max([r for r in recs if r["group"] == g], key=lambda r: r["mc"][seed0]["worst"])["id"]}
                          for g in ("SLS", "CRUISE")}
    kap = [r["kappaV"] for r in recs]
    tau = [r["tau_slow"] for r in recs]
    from scipy.stats import spearmanr
    res = {"points": recs, "seed_spread": spread, "idle_over_cruise_worst_shift": ratio,
           "idle_id": idle["id"], "cruise_id": cruise["id"], "E_norm": cfg.E_norm, "n_mc": cfg.n_mc,
           "spearman_kappa_vs_tau_slow": float(spearmanr(kap, tau).statistic),
           "verdicts": {"H2_idle_most_fragile": h2}, "fast": cfg.fast,
           "form": "normalized states: A_n = D^-1 A D, D = diag(x0)"}
    print("H2 idle most fragile: %s  (slowest %s, most fragile %s; SLS idle / max-cruise worst-shift ratio %s)" % (
        h2["verdict"], h2["per_seed"][seeds[0]]["slowest_point"], h2["per_seed"][seeds[0]]["most_fragile_point"],
        ", ".join("%.2f" % v for v in ratio.values())))
    write_json("layer3_bauerfike", res, ck)
    return ck


if __name__ == "__main__":
    c = config.from_argv()
    c.ensure_dirs()
    finish(run(c))
