#!/usr/bin/env python3
"""Biological Scenario: Tumor Margin Invasion

Scientific Question:
    Does ecDNA-driven EGFR amplification causally drive invasive behavior
    at the tumor margin, or is the association confounded by the hypoxic
    microenvironment that both selects for high-ecDNA cells and promotes migration?

Hypothesis:
    If ecDNA causally drives invasion:
    - IV estimate of migration effect (delta) should be positive
    - Effect should persist after controlling for hypoxia
    - Margin cells should show causal effect, not just correlation

Validation:
    - Compare IV vs OLS estimates (OLS biased by hypoxia confounding)
    - Stratify by tumor region (core vs margin vs infiltrating)
    - Test whether effect varies with distance from tumor center

Usage:
    python scripts/scenarios/margin_invasion.py --output-dir results/margin_invasion
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np


def run_margin_invasion_study(
    base_config: Path,
    output_dir: Path,
    n_replicates: int = 10,
) -> dict:
    """Run margin invasion study.

    Simulates tumors with varying migration effects and analyzes
    whether IV can identify the causal effect of ecDNA on invasion.
    """
    from causanta.simulate.sweep import SweepConfig, run_sweep
    from causanta.analyze.heterogeneity import stratified_iv_by_region
    from causanta.analyze.bootstrap import bootstrap_iv_estimate

    print("=" * 60)
    print("BIOLOGICAL SCENARIO: Tumor Margin Invasion")
    print("=" * 60)
    print("\nQuestion: Does ecDNA causally drive invasion at tumor margins?")
    print("Hypothesis: IV should detect causal effect, OLS will be biased by hypoxia")

    # Vary the migration effect parameter
    sweep_config = SweepConfig(
        name="margin_invasion",
        base_config_path=base_config,
        output_dir=output_dir,
        n_replicates=n_replicates,
        vary_params={
            "cell_types.6.ecDNA_effect_on_migration": [0.0, 0.025, 0.05, 0.075, 0.10],
        },
    )

    print(f"\nRunning {sweep_config.total_runs()} simulations...")

    def progress(i, total, result):
        status = "OK" if result["status"] == "success" else "FAIL"
        print(f"  [{i}/{total}] {result['name']}: {status}")

    results = run_sweep(sweep_config, n_workers=4, progress_callback=progress)

    # Analyze results
    print("\n" + "=" * 60)
    print("ANALYSIS")
    print("=" * 60)

    analysis_results = {
        "scenario": "margin_invasion",
        "question": "Does ecDNA causally drive invasion?",
        "conditions": [],
    }

    for result in results:
        if result["status"] != "success":
            continue

        run_dir = Path(result["output_dir"])
        cells_file = run_dir / "cells_final.json"

        if not cells_file.exists():
            continue

        with open(cells_file) as f:
            cells_data = json.load(f)

        tumor_cells = [c for c in cells_data if c.get("cell_type") == 6]

        if len(tumor_cells) < 50:
            continue

        # Extract data
        Z = np.array([c["ecDNA_count"] for c in tumor_cells], dtype=float)
        D = np.array([c.get("egfr_expression", c["ecDNA_count"]) for c in tumor_cells], dtype=float)
        Y = np.array([c.get("migration_rate", 0) for c in tumor_cells], dtype=float)

        # Bootstrap IV estimate
        try:
            bootstrap_result = bootstrap_iv_estimate(Z, D, Y, n_bootstrap=500)

            true_effect = result["params"].get("cell_types.6.ecDNA_effect_on_migration", 0.05)

            condition_result = {
                "true_delta": true_effect,
                "iv_estimate": bootstrap_result.point_estimate,
                "iv_ci_lower": bootstrap_result.ci_lower,
                "iv_ci_upper": bootstrap_result.ci_upper,
                "ols_estimate": bootstrap_result.ols_estimate,
                "bias_iv": abs(bootstrap_result.point_estimate - true_effect),
                "bias_ols": abs(bootstrap_result.ols_estimate - true_effect),
                "covers_truth": bootstrap_result.ci_lower <= true_effect <= bootstrap_result.ci_upper,
                "n_cells": len(tumor_cells),
            }

            analysis_results["conditions"].append(condition_result)

            print(f"\nTrue delta = {true_effect:.3f}")
            print(f"  IV:  {bootstrap_result.point_estimate:.4f} [{bootstrap_result.ci_lower:.4f}, {bootstrap_result.ci_upper:.4f}]")
            print(f"  OLS: {bootstrap_result.ols_estimate:.4f}")
            print(f"  IV bias:  {condition_result['bias_iv']:.4f}")
            print(f"  OLS bias: {condition_result['bias_ols']:.4f}")

        except Exception as e:
            print(f"  Analysis failed: {e}")

    # Stratified analysis by region
    print("\n" + "-" * 40)
    print("Stratified Analysis by Tumor Region")
    print("-" * 40)

    # Use last successful run for stratified analysis
    for result in reversed(results):
        if result["status"] == "success":
            run_dir = Path(result["output_dir"])
            cells_file = run_dir / "cells_final.json"
            if cells_file.exists():
                with open(cells_file) as f:
                    cells_data = json.load(f)

                try:
                    het_results = stratified_iv_by_region(
                        cells_data,
                        outcome_name="migration_rate",
                        center_x=500,
                        center_y=500,
                    )

                    print("\nRegion-stratified IV estimates:")
                    for region in ["core", "margin", "infiltrating"]:
                        if region in het_results.get("by_region", {}):
                            r = het_results["by_region"][region]
                            print(f"  {region}: {r.get('iv_estimate', 0):.4f} (n={r.get('n_cells', 0)})")

                    analysis_results["heterogeneity"] = het_results

                except Exception as e:
                    print(f"  Stratified analysis failed: {e}")

                break

    # Summary
    print("\n" + "=" * 60)
    print("CONCLUSIONS")
    print("=" * 60)

    if analysis_results["conditions"]:
        mean_iv_bias = np.mean([c["bias_iv"] for c in analysis_results["conditions"]])
        mean_ols_bias = np.mean([c["bias_ols"] for c in analysis_results["conditions"]])
        coverage = np.mean([c["covers_truth"] for c in analysis_results["conditions"]])

        print(f"\n1. Mean IV bias:  {mean_iv_bias:.4f}")
        print(f"2. Mean OLS bias: {mean_ols_bias:.4f}")
        print(f"3. CI coverage:   {coverage*100:.1f}%")

        if mean_iv_bias < mean_ols_bias * 0.5:
            print("\nResult: IV successfully identifies causal effect of ecDNA on invasion")
            print("        while OLS is biased by hypoxia confounding.")
            analysis_results["conclusion"] = "IV_SUPERIOR"
        else:
            print("\nResult: Both methods show similar bias; may need larger sample size")
            analysis_results["conclusion"] = "INCONCLUSIVE"

    # Save results
    results_file = output_dir / "margin_invasion_results.json"
    with open(results_file, "w") as f:
        json.dump(analysis_results, f, indent=2, default=str)

    print(f"\nResults saved to: {results_file}")

    return analysis_results


def main():
    parser = argparse.ArgumentParser(description="Tumor margin invasion scenario")
    parser.add_argument("--config", type=Path, default=Path("causanta/simulate/params/default.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/margin_invasion"))
    parser.add_argument("--n-replicates", type=int, default=10)

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    run_margin_invasion_study(args.config, args.output_dir, args.n_replicates)


if __name__ == "__main__":
    main()
