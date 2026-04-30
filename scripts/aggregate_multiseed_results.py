#!/usr/bin/env python3
"""Aggregate multi-seed reliability results across (scale, scenario) cells.

Produces output/multiseed_full_results.json with per-seed and aggregate
(mean, SD) statistics for IV β, OLS β, OLS bias %, PC F1/precision/recall,
first-stage F, and tumor cell count, across:

  - 2mm × {baseline, reduced, removed}
  - 6mm × {baseline, reduced, removed}

Each cell uses 5 seeds: the existing seed-42 run from the
{large,xlarge}_{baseline,reduced,removed} directories, plus seeds 43-46
from the multiseed_{2mm,6mm}_{scenario}_seed{N} runs.

Run after both batch controllers have written their .done markers:
    output/multiseed_logs/2mm_batch.done
    output/multiseed_logs/6mm_batch.done
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from causanta.analyze.discovery import compare_to_ground_truth, discover_causal_structure
from causanta.analyze.iv import estimate_iv_effects
from causanta.analyze.loader import load_simulation_output

GROUND_TRUTH = [
    ("ecDNA_count", "egfr_expression"),
    ("egfr_expression", "VEGF_secretion"),
    ("egfr_expression", "migration_rate"),
    ("O2_local", "VEGF_secretion"),
]


def analyze_run(run_dir: str, seed: int) -> dict | None:
    if not Path(run_dir).exists():
        return None
    data = load_simulation_output(run_dir)
    n_cells = len(data.tumor_cells)
    if n_cells < 100:
        return None
    iv_results = estimate_iv_effects(data, outcomes=["VEGF_secretion"])
    if not iv_results:
        return None
    iv = iv_results[0]
    dag = discover_causal_structure(data, algorithm="pc", alpha=0.05)
    metrics = compare_to_ground_truth(dag, GROUND_TRUTH)
    bias_pct = (iv.bias_ratio - 1) * 100 if iv.bias_ratio != float("inf") else 0
    return {
        "seed": seed,
        "n_cells": n_cells,
        "ols_beta": round(iv.ols_coef, 3),
        "iv_beta": round(iv.second_stage_coef, 3),
        "ols_bias_pct": round(bias_pct, 1),
        "first_stage_F": round(iv.first_stage_f_stat, 0),
        "precision": round(metrics["precision"], 3),
        "recall": round(metrics["recall"], 3),
        "f1_score": round(metrics["f1_score"], 3),
    }


def aggregate(per_seed: list[dict]) -> dict:
    if not per_seed:
        return {}
    keys = ["n_cells", "ols_beta", "iv_beta", "ols_bias_pct", "precision", "recall", "f1_score"]
    out = {}
    for k in keys:
        vals = np.array([r[k] for r in per_seed])
        out[f"{k}_mean"] = float(vals.mean())
        out[f"{k}_sd"] = float(vals.std(ddof=1)) if len(vals) > 1 else 0.0
    return out


def main() -> int:
    # (scale, scenario) -> [(run_dir, seed)]
    cells: dict[tuple[str, str], list[tuple[str, int]]] = {
        ("2mm", "baseline"): [
            ("output/large_baseline", 42),
            *(("output/multiseed_2mm/seed_%d" % s, s) for s in [43, 44, 45, 46]),
        ],
        ("2mm", "reduced"): [
            ("output/large_reduced", 42),
            *(("output/multiseed_2mm/reduced_seed_%d" % s, s) for s in [43, 44, 45, 46]),
        ],
        ("2mm", "removed"): [
            ("output/large_removed", 42),
            *(("output/multiseed_2mm/removed_seed_%d" % s, s) for s in [43, 44, 45, 46]),
        ],
        ("6mm", "baseline"): [
            ("output/xlarge_baseline", 42),
            *(("output/multiseed_6mm/baseline_seed_%d" % s, s) for s in [43, 44, 45, 46]),
        ],
        ("6mm", "reduced"): [
            ("output/xlarge_reduced", 42),
            *(("output/multiseed_6mm/reduced_seed_%d" % s, s) for s in [43, 44, 45, 46]),
        ],
        ("6mm", "removed"): [
            ("output/xlarge_removed", 42),
            *(("output/multiseed_6mm/removed_seed_%d" % s, s) for s in [43, 44, 45, 46]),
        ],
    }

    results: dict[str, dict] = {}
    for (scale, scenario), runs in cells.items():
        per_seed = []
        for run_dir, seed in runs:
            r = analyze_run(run_dir, seed)
            if r is not None:
                per_seed.append(r)
        cell_key = f"{scale}_{scenario}"
        results[cell_key] = {
            "scale": scale,
            "scenario": scenario,
            "n_seeds": len(per_seed),
            "per_seed": per_seed,
            "aggregate": aggregate(per_seed),
        }

    out_path = Path("output/multiseed_full_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print("=" * 80)
    print(f"{'Cell':<14} {'n':>3} {'IV β':>14} {'OLS β':>14} {'F1':>14}")
    print("-" * 80)
    for key, cell in results.items():
        a = cell["aggregate"]
        if not a:
            print(f"{key:<14} {cell['n_seeds']:>3} -- no usable runs --")
            continue
        iv_str = f"{a['iv_beta_mean']:.3f}±{a['iv_beta_sd']:.3f}"
        ols_str = f"{a['ols_beta_mean']:.3f}±{a['ols_beta_sd']:.3f}"
        f1_str = f"{a['f1_score_mean']:.3f}±{a['f1_score_sd']:.3f}"
        print(f"{key:<14} {cell['n_seeds']:>3} {iv_str:>14} {ols_str:>14} {f1_str:>14}")
    print("=" * 80)
    print(f"Saved: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
