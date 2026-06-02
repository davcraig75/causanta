#!/usr/bin/env python3
"""Aggregate discovery_sweep results across the Tier C multi-seed runs.

Loads each output/run_<TS> created during the multi-seed sweep, runs PC and
GES on each, and reports per-seed metrics plus mean/std summary.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from causanta.analyze.discovery import (
    compare_to_ground_truth,
    ges_algorithm,
    pc_algorithm,
)
from causanta.analyze.loader import cells_to_array, load_simulation_output

OBSERVED = [
    "ecDNA_count",
    "egfr_expression",
    "O2_local",
    "glucose_local",
    "VEGF_secretion",
    "migration_rate",
]

GT_PROJECTED = [
    ("ecDNA_count", "egfr_expression"),
    ("O2_local", "egfr_expression"),
    ("egfr_expression", "VEGF_secretion"),
    ("O2_local", "VEGF_secretion"),
    ("egfr_expression", "migration_rate"),
    ("O2_local", "migration_rate"),
]


def discover_one(run_dir: Path) -> dict | None:
    data_dir = run_dir / "data"
    if not (data_dir / "lineage.tsv").exists():
        return None
    data = load_simulation_output(data_dir)
    if not data.cells:
        return None
    final_t = max(data.cells.keys())
    tumor = [c for c in data.cells[final_t] if int(c.get("cell_type", -1)) == 6]
    if len(tumor) < 30:
        return {"run": str(run_dir.name), "error": "few_tumor", "n_tumor": len(tumor)}
    X = cells_to_array(tumor, OBSERVED).astype(float)
    std = X.std(axis=0)
    std[std == 0] = 1.0
    Xs = (X - X.mean(axis=0)) / std

    # O2 hypoxia fraction
    o2_idx = OBSERVED.index("O2_local")
    hypoxic_frac = float((X[:, o2_idx] < 18.0).mean())
    o2_egfr_corr = float(
        np.corrcoef(X[:, o2_idx], X[:, OBSERVED.index("egfr_expression")])[0, 1]
    )

    pc = pc_algorithm(Xs, OBSERVED, alpha=0.05, max_cond_size=3)
    ges = ges_algorithm(Xs, OBSERVED, penalty=1.0)
    pc_m = compare_to_ground_truth(pc, GT_PROJECTED)
    ges_m = compare_to_ground_truth(ges, GT_PROJECTED)

    return {
        "run": run_dir.name,
        "n_tumor": len(tumor),
        "hypoxic_frac": hypoxic_frac,
        "corr_O2_EGFR": o2_egfr_corr,
        "pc": pc_m,
        "ges": ges_m,
    }


def fmt_summary(results: list[dict], algo: str) -> str:
    f1s = [r[algo]["f1_score"] for r in results if algo in r]
    precs = [r[algo]["precision"] for r in results if algo in r]
    recs = [r[algo]["recall"] for r in results if algo in r]
    shds = [r[algo]["structural_hamming_distance"] for r in results if algo in r]
    if not f1s:
        return f"{algo}: no valid runs"
    arr = np.array
    return (
        f"{algo}: n={len(f1s)}  "
        f"F1 = {arr(f1s).mean():.3f} +/- {arr(f1s).std():.3f}  "
        f"P = {arr(precs).mean():.3f} +/- {arr(precs).std():.3f}  "
        f"R = {arr(recs).mean():.3f} +/- {arr(recs).std():.3f}  "
        f"SHD = {arr(shds).mean():.2f} +/- {arr(shds).std():.2f}"
    )


def main() -> None:
    base = Path("output")
    runs = sorted(p for p in base.glob("run_*") if p.is_dir())

    # Restrict to runs that started during the multi-seed sweep window if a
    # since-marker file is present; otherwise take all run_* dirs and let the
    # caller filter.
    print(f"Found {len(runs)} candidate run dirs")

    results: list[dict] = []
    print(f"\n{'run':<28} {'n_tumor':>7} {'hypox':>6} {'r(O2,E)':>8}  {'PC F1':>6} {'PC P':>6} {'PC R':>6}  {'GES F1':>6}")
    for run in runs:
        r = discover_one(run)
        if r is None:
            continue
        if "error" in r:
            print(f"{r['run']:<28} {r.get('n_tumor', 0):>7} (skip: {r['error']})")
            continue
        results.append(r)
        print(
            f"{r['run']:<28} {r['n_tumor']:>7} {r['hypoxic_frac']*100:>5.1f}% "
            f"{r['corr_O2_EGFR']:>+8.3f}  "
            f"{r['pc']['f1_score']:>6.3f} "
            f"{r['pc']['precision']:>6.3f} "
            f"{r['pc']['recall']:>6.3f}  "
            f"{r['ges']['f1_score']:>6.3f}"
        )

    print("\n## Summary across runs ##")
    print(fmt_summary(results, "pc"))
    print(fmt_summary(results, "ges"))

    out = Path("output/multiseed_discovery.json")
    out.write_text(json.dumps(results, indent=2, default=str))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
