#!/usr/bin/env python3
"""Tier A: discovery hyperparameter sweep + bootstrap stability against the
projected ground-truth DAG.

The default GT in run_comprehensive_analysis.py omits the O2 -> EGFR edge that
arises from HIF-1alpha upregulation (EGFR_HYPOXIA_UPREGULATION=1.5 in
ecdna.py). It also penalises the algorithm for any glucose edge despite the
fact that glucose only gates proliferation (no direct edge to VEGF/migration
in cell-record space). This sweep uses a defensible 6-edge projected GT and
reports F1 across PC alpha and GES penalty values, with optional bootstrap
edge stability.

Usage:
    python scripts/discovery_sweep.py output/run_<TS>/data/
"""

from __future__ import annotations

import argparse
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

# Projected ground truth in observed-variable subspace.
# Direct mechanistic edges per causanta/simulate/ecdna.py and behaviors.py:
#   ecDNA -> EGFR via gene dosage (kappa)
#   O2    -> EGFR via HIF-1a upregulation (hypoxia indicator)
#   EGFR  -> VEGF (sqrt, beta)
#   O2    -> VEGF via hypoxia (multiplicative phi)
#   EGFR  -> migration (multiplicative, delta)
#   O2    -> migration via hypoxia (Go-or-Grow boost)
GT_PROJECTED = [
    ("ecDNA_count", "egfr_expression"),
    ("O2_local", "egfr_expression"),
    ("egfr_expression", "VEGF_secretion"),
    ("O2_local", "VEGF_secretion"),
    ("egfr_expression", "migration_rate"),
    ("O2_local", "migration_rate"),
]


def get_tumor_data(run_data_dir: Path) -> tuple[list[dict], np.ndarray]:
    data = load_simulation_output(run_data_dir)
    final_t = max(data.cells.keys())
    cells = data.cells[final_t]
    tumor = [c for c in cells if int(c.get("cell_type", -1)) == 6]
    if not tumor:
        return [], np.zeros((0, len(OBSERVED)))
    X = cells_to_array(tumor, OBSERVED).astype(float)
    # Standardise (matches discover_causal_structure)
    std = X.std(axis=0)
    std[std == 0] = 1.0
    X = (X - X.mean(axis=0)) / std
    return tumor, X


def score(disc, gt) -> dict:
    return compare_to_ground_truth(disc, gt)


def fmt_row(label, m, width=14):
    return (
        f"{label:<{width}} "
        f"{m['true_positives']:>3} "
        f"{m['false_positives']:>3} "
        f"{m['false_negatives']:>3} "
        f"{m['precision']:>6.3f} "
        f"{m['recall']:>6.3f} "
        f"{m['f1_score']:>6.3f} "
        f"{m['structural_hamming_distance']:>4}"
    )


HEADER = (
    f"{'config':<14} "
    f"{'TP':>3} {'FP':>3} {'FN':>3} "
    f"{'P':>6} {'R':>6} {'F1':>6} {'SHD':>4}"
)


def bootstrap_pc_stability(
    X: np.ndarray,
    alpha: float,
    n_boot: int = 50,
    threshold: float = 0.7,
    seed: int = 0,
) -> list[tuple[str, str]]:
    rng = np.random.default_rng(seed)
    n = X.shape[0]
    edge_counts: dict[tuple[str, str], int] = {}
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        Xb = X[idx]
        result = pc_algorithm(Xb, OBSERVED, alpha=alpha, max_cond_size=3)
        seen: set[tuple[str, str]] = set()
        for e in result.edges:
            seen.add(e)
        for a, b in result.undirected:
            seen.add((a, b))
            seen.add((b, a))
        for e in seen:
            edge_counts[e] = edge_counts.get(e, 0) + 1
    return [e for e, c in edge_counts.items() if c / n_boot >= threshold]


def score_edges(edges, gt):
    gt_set = set(gt)
    disc_set = set(edges)
    tp = len(gt_set & disc_set)
    fp = len(disc_set - gt_set)
    fn = len(gt_set - disc_set)
    prec = tp / (tp + fp) if (tp + fp) else 0
    rec = tp / (tp + fn) if (tp + fn) else 0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0
    return {
        "true_positives": tp,
        "false_positives": fp,
        "false_negatives": fn,
        "precision": prec,
        "recall": rec,
        "f1_score": f1,
        "structural_hamming_distance": fp + fn,
    }


def run_sweep(run_data_dir: Path, n_boot: int = 50) -> dict:
    tumor, X = get_tumor_data(run_data_dir)
    print(f"# Run: {run_data_dir}")
    print(f"# Tumor cells: {len(tumor)} | data matrix {X.shape}")
    print(f"# Projected GT ({len(GT_PROJECTED)} edges):")
    for s, t in GT_PROJECTED:
        print(f"#   {s} -> {t}")
    print()

    print(HEADER)
    out: dict = {"pc": [], "ges": [], "pc_bootstrap": None}

    for alpha in [0.001, 0.01, 0.05]:
        result = pc_algorithm(X, OBSERVED, alpha=alpha, max_cond_size=3)
        m = score(result, GT_PROJECTED)
        print(fmt_row(f"PC a={alpha}", m))
        out["pc"].append({"alpha": alpha, **m})

    for penalty in [1.0, 2.0, 3.0, 5.0, 8.0]:
        result = ges_algorithm(X, OBSERVED, penalty=penalty)
        m = score(result, GT_PROJECTED)
        print(fmt_row(f"GES p={penalty}", m))
        out["ges"].append({"penalty": penalty, **m})

    if n_boot > 0 and X.shape[0] >= 30:
        # Pick best PC alpha by F1 to bootstrap
        best = max(out["pc"], key=lambda r: r["f1_score"])
        alpha_best = best["alpha"]
        stable = bootstrap_pc_stability(X, alpha_best, n_boot=n_boot, threshold=0.7)
        m = score_edges(stable, GT_PROJECTED)
        print(fmt_row(f"PC boot a={alpha_best}", m))
        print(f"#   stable edges (>=70% in {n_boot} resamples): {sorted(stable)}")
        out["pc_bootstrap"] = {
            "alpha": alpha_best,
            "n_boot": n_boot,
            "threshold": 0.7,
            "stable_edges": stable,
            **m,
        }

    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("run_dir", help="Path to a run's data directory")
    p.add_argument("--n_boot", type=int, default=50)
    args = p.parse_args()
    run_sweep(Path(args.run_dir), n_boot=args.n_boot)


if __name__ == "__main__":
    main()
