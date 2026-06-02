#!/usr/bin/env python3
"""Run IV + PC discovery at multiple intermediate timesteps for each
(scale, scenario, seed) combination.

Produces output/multiseed_timeseries_results.json with per-timestep
hypoxic fraction, IV/OLS estimates, and discovery F1, so the manuscript
can characterize how confounder visibility evolves with tumor growth.

Requires the full per-hour snapshot time series (i.e., the simulations
must have been run with the cleanup-off default introduced in commit
0970608).
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from causanta.analyze.discovery import (
    compare_to_ground_truth,
    pc_algorithm,
)

GT_WITH_HIF = [
    ("ecDNA_count", "egfr_expression"),
    ("is_hypoxic", "egfr_expression"),
    ("egfr_expression", "VEGF_secretion"),
    ("is_hypoxic", "VEGF_secretion"),
    ("egfr_expression", "migration_rate"),
]
GT_REMOVED = [
    ("ecDNA_count", "egfr_expression"),
    ("egfr_expression", "VEGF_secretion"),
    ("is_hypoxic", "VEGF_secretion"),
    ("egfr_expression", "migration_rate"),
]

VARS = ["ecDNA_count", "egfr_expression", "is_hypoxic", "VEGF_secretion", "migration_rate"]
KEEP_COLS = set(VARS + ["cell_type"])


def load_tumor_subset(path: Path) -> np.ndarray | None:
    """Stream-load only the columns we need and only tumor cells."""
    if not path.exists():
        return None
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if r.get("cell_type") != "6":
                continue
            row = []
            for v in VARS:
                val = r[v]
                if val in ("True", "False"):
                    row.append(1.0 if val == "True" else 0.0)
                else:
                    row.append(float(val))
            rows.append(row)
    if not rows:
        return None
    return np.array(rows)


def fit_ols_iv(Z: np.ndarray, D: np.ndarray, Y: np.ndarray) -> tuple[float, float, float]:
    """Simple OLS regression of Y on D, and 2SLS using Z as instrument."""
    # OLS
    X_ols = np.column_stack([np.ones(len(D)), D])
    b_ols, *_ = np.linalg.lstsq(X_ols, Y, rcond=None)
    ols = float(b_ols[1])
    # 2SLS first stage
    X1 = np.column_stack([np.ones(len(Z)), Z])
    pi, *_ = np.linalg.lstsq(X1, D, rcond=None)
    D_hat = X1 @ pi
    # Second stage
    X2 = np.column_stack([np.ones(len(D_hat)), D_hat])
    b_iv, *_ = np.linalg.lstsq(X2, Y, rcond=None)
    iv = float(b_iv[1])
    # First-stage F
    resid = D - D_hat
    ss_res = float(np.sum(resid**2))
    ss_tot = float(np.sum((D - D.mean())**2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    n = len(Z); k = 1
    fstat = (r2 / k) / ((1 - r2) / (n - k - 1)) if r2 < 1 else float("inf")
    return ols, iv, float(fstat)


def analyze_timestep(run_dir: Path, t: int, scenario: str) -> dict | None:
    path = run_dir / "data" / f"cells_t{t:06d}.tsv"
    M = load_tumor_subset(path)
    if M is None or len(M) < 100:
        return None
    Z = M[:, 0]  # ecDNA_count
    D = M[:, 1]  # egfr_expression
    H = M[:, 2]  # is_hypoxic (0/1)
    V = M[:, 3]  # VEGF_secretion
    ols, iv, fstat = fit_ols_iv(Z, D, V)
    bias_pct = (ols / iv - 1) * 100 if abs(iv) > 1e-10 else 0
    # Discovery
    std = M.std(axis=0); std[std == 0] = 1
    Ms = (M - M.mean(axis=0)) / std
    dag = pc_algorithm(Ms, VARS, alpha=0.05, max_cond_size=3)
    gt = GT_REMOVED if scenario == "removed" else GT_WITH_HIF
    metrics = compare_to_ground_truth(dag, gt)
    return {
        "t": t,
        "n_tumor": int(len(M)),
        "hypoxic_fraction": float(H.mean()),
        "ols_beta": round(ols, 3),
        "iv_beta": round(iv, 3),
        "ols_bias_pct": round(bias_pct, 1),
        "first_stage_F": round(fstat, 0),
        "precision": round(metrics["precision"], 3),
        "recall": round(metrics["recall"], 3),
        "f1": round(metrics["f1_score"], 3),
    }


def main() -> int:
    # (scale, scenario) -> [(run_dir, seed)]
    cells = {
        ("2mm", "baseline"): [("output/large_baseline", 42)] + [(f"output/multiseed_2mm/seed_{s}", s) for s in [43, 44, 45, 46]],
        ("2mm", "reduced"):  [("output/large_reduced", 42)] + [(f"output/multiseed_2mm/reduced_seed_{s}", s) for s in [43, 44, 45, 46]],
        ("2mm", "removed"):  [("output/large_removed", 42)] + [(f"output/multiseed_2mm/removed_seed_{s}", s) for s in [43, 44, 45, 46]],
        ("6mm", "baseline"): [("output/xlarge_baseline", 42)] + [(f"output/multiseed_6mm/baseline_seed_{s}", s) for s in [43, 44, 45, 46]],
        ("6mm", "reduced"):  [("output/xlarge_reduced", 42)] + [(f"output/multiseed_6mm/reduced_seed_{s}", s) for s in [43, 44, 45, 46]],
        ("6mm", "removed"):  [("output/xlarge_removed", 42)] + [(f"output/multiseed_6mm/removed_seed_{s}", s) for s in [43, 44, 45, 46]],
    }
    timesteps = {
        "2mm": [80, 120, 160, 200, 240],
        "6mm": [100, 150, 200, 250, 300],
    }
    out: dict[str, dict] = {}
    for (scale, scenario), runs in cells.items():
        cell_key = f"{scale}_{scenario}"
        out[cell_key] = {"scale": scale, "scenario": scenario, "by_t": {}}
        for t in timesteps[scale]:
            per_seed = []
            for run_dir, seed in runs:
                r = analyze_timestep(Path(run_dir), t, scenario)
                if r is not None:
                    per_seed.append({"seed": seed, **r})
            if not per_seed:
                continue
            # aggregate
            def agg(k):
                vals = np.array([s[k] for s in per_seed])
                return {"mean": float(vals.mean()), "sd": float(vals.std(ddof=1)) if len(vals) > 1 else 0.0}
            out[cell_key]["by_t"][str(t)] = {
                "n_seeds": len(per_seed),
                "per_seed": per_seed,
                "agg": {k: agg(k) for k in ["n_tumor", "hypoxic_fraction", "ols_beta", "iv_beta", "ols_bias_pct", "first_stage_F", "precision", "recall", "f1"]},
            }
        print(f"{cell_key}:")
        for t_str, d in out[cell_key]["by_t"].items():
            a = d["agg"]
            print(f"  t={t_str}  hypoxic={a['hypoxic_fraction']['mean']:.3f}±{a['hypoxic_fraction']['sd']:.3f}  "
                  f"IVβ={a['iv_beta']['mean']:.3f}±{a['iv_beta']['sd']:.3f}  "
                  f"bias={a['ols_bias_pct']['mean']:+.1f}±{a['ols_bias_pct']['sd']:.1f}%  "
                  f"F1={a['f1']['mean']:.3f}±{a['f1']['sd']:.3f}")

    out_path = Path("output/multiseed_timeseries_results.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
