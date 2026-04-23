#!/usr/bin/env python3
"""Biological Scenario: Driver vs Passenger ecDNA

Scientific Question:
    Can IV distinguish between ecDNA that causally drives tumor phenotypes
    (driver) and ecDNA that is passively amplified without functional
    consequences (passenger)?

Hypothesis:
    - Driver ecDNA should show positive IV estimates for phenotypic outcomes
    - Passenger ecDNA (effect = 0) should show null IV estimates
    - OLS may show spurious associations for passengers due to confounding

Simulation Design:
    1. Simulate with varying effect sizes (0.0 = passenger, >0 = driver)
    2. Use IV to estimate effects
    3. Compare detection of true nulls vs true effects

This addresses a key clinical question: not all ecDNA is functionally
relevant, and distinguishing drivers from passengers is crucial for
precision oncology.

Usage:
    python scripts/scenarios/driver_vs_passenger.py --output-dir results/driver_passenger
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np


def run_driver_passenger_study(
    base_config: Path,
    output_dir: Path,
    n_replicates: int = 10,
) -> dict:
    """Run driver vs passenger discrimination study.

    Tests whether IV can correctly identify:
    1. True drivers (effect > 0): positive significant estimate
    2. True passengers (effect = 0): null estimate
    """
    from causanta.simulate.sweep import SweepConfig, run_sweep
    from causanta.analyze.bootstrap import bootstrap_iv_estimate

    print("=" * 60)
    print("BIOLOGICAL SCENARIO: Driver vs Passenger ecDNA")
    print("=" * 60)
    print("\nQuestion: Can IV distinguish driver from passenger ecDNA?")
    print("Key: Passengers have effect = 0, drivers have effect > 0")

    # Define effect sizes to test
    # 0.0 = passenger (no functional effect)
    # 0.05-0.20 = driver (functional effect)
    effect_sizes = [0.0, 0.02, 0.05, 0.10, 0.20]

    all_results = []

    for effect in effect_sizes:
        is_driver = effect > 0
        label = "DRIVER" if is_driver else "PASSENGER"

        print(f"\n{'='*40}")
        print(f"Testing: {label} (effect = {effect})")
        print("=" * 40)

        sweep_config = SweepConfig(
            name=f"effect_{effect}",
            base_config_path=base_config,
            output_dir=output_dir / f"effect_{effect}",
            n_replicates=n_replicates,
            fixed_params={
                "cell_types.6.ecDNA_effect_on_VEGF": effect,
            },
        )

        def progress(i, total, result):
            status = "OK" if result["status"] == "success" else "FAIL"
            print(f"  [{i}/{total}] {result['name']}: {status}")

        results = run_sweep(sweep_config, n_workers=4, progress_callback=progress)

        # Analyze each replicate
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
            Y = np.array([c.get("VEGF_secretion", 0) for c in tumor_cells], dtype=float)

            # Bootstrap IV estimate
            try:
                bootstrap_result = bootstrap_iv_estimate(Z, D, Y, n_bootstrap=500)

                # Determine if we detected an effect
                # Effect detected if CI excludes 0
                ci_excludes_zero = (bootstrap_result.ci_lower > 0) or (bootstrap_result.ci_upper < 0)

                # Correct classification
                # True positive: driver detected as driver
                # True negative: passenger detected as passenger
                # False positive: passenger detected as driver (Type I error)
                # False negative: driver detected as passenger (Type II error)

                if is_driver and ci_excludes_zero:
                    classification = "TRUE_POSITIVE"
                elif is_driver and not ci_excludes_zero:
                    classification = "FALSE_NEGATIVE"
                elif not is_driver and not ci_excludes_zero:
                    classification = "TRUE_NEGATIVE"
                else:  # passenger but CI excludes zero
                    classification = "FALSE_POSITIVE"

                all_results.append({
                    "true_effect": effect,
                    "is_driver": is_driver,
                    "iv_estimate": bootstrap_result.point_estimate,
                    "iv_ci_lower": bootstrap_result.ci_lower,
                    "iv_ci_upper": bootstrap_result.ci_upper,
                    "ols_estimate": bootstrap_result.ols_estimate,
                    "ci_excludes_zero": ci_excludes_zero,
                    "classification": classification,
                    "n_cells": len(tumor_cells),
                })

            except Exception as e:
                print(f"    Analysis failed: {e}")

    # Summary analysis
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)

    analysis_results = {
        "scenario": "driver_vs_passenger",
        "question": "Can IV distinguish driver from passenger ecDNA?",
        "replicates": all_results,
    }

    # Compute classification metrics
    classifications = [r["classification"] for r in all_results]
    tp = classifications.count("TRUE_POSITIVE")
    tn = classifications.count("TRUE_NEGATIVE")
    fp = classifications.count("FALSE_POSITIVE")
    fn = classifications.count("FALSE_NEGATIVE")

    total = len(classifications)
    drivers = sum(1 for r in all_results if r["is_driver"])
    passengers = total - drivers

    print(f"\nClassification Results (n = {total}):")
    print(f"  True Positives (drivers detected):     {tp}")
    print(f"  True Negatives (passengers detected):  {tn}")
    print(f"  False Positives (Type I error):        {fp}")
    print(f"  False Negatives (Type II error):       {fn}")

    # Metrics
    sensitivity = tp / drivers if drivers > 0 else 0  # Power
    specificity = tn / passengers if passengers > 0 else 0  # 1 - Type I error rate
    fpr = fp / passengers if passengers > 0 else 0  # Type I error rate

    print(f"\nMetrics:")
    print(f"  Sensitivity (Power):  {sensitivity*100:.1f}%")
    print(f"  Specificity:          {specificity*100:.1f}%")
    print(f"  False Positive Rate:  {fpr*100:.1f}%")

    analysis_results["metrics"] = {
        "sensitivity": sensitivity,
        "specificity": specificity,
        "false_positive_rate": fpr,
        "true_positives": tp,
        "true_negatives": tn,
        "false_positives": fp,
        "false_negatives": fn,
    }

    # By effect size
    print("\n" + "-" * 40)
    print("Detection Rate by Effect Size")
    print("-" * 40)

    for effect in effect_sizes:
        effect_results = [r for r in all_results if r["true_effect"] == effect]
        if effect_results:
            detected = sum(1 for r in effect_results if r["ci_excludes_zero"])
            rate = detected / len(effect_results)
            mean_iv = np.mean([r["iv_estimate"] for r in effect_results])
            mean_ols = np.mean([r["ols_estimate"] for r in effect_results])
            label = "DRIVER" if effect > 0 else "PASSENGER"
            print(f"  {label} (β={effect}): {rate*100:5.1f}% detected | IV={mean_iv:.4f} OLS={mean_ols:.4f}")

    # Conclusions
    print("\n" + "=" * 60)
    print("CONCLUSIONS")
    print("=" * 60)

    if specificity >= 0.90 and sensitivity >= 0.70:
        print("\nResult: IV successfully discriminates drivers from passengers")
        print("        - Low false positive rate (good specificity)")
        print("        - Reasonable power to detect true effects")
        analysis_results["conclusion"] = "SUCCESSFUL_DISCRIMINATION"
    elif specificity >= 0.90:
        print("\nResult: IV has good specificity but limited power")
        print("        - Low false positive rate (passengers not misclassified)")
        print("        - May need larger sample size to detect weak drivers")
        analysis_results["conclusion"] = "GOOD_SPECIFICITY_LOW_POWER"
    elif sensitivity >= 0.80:
        print("\nResult: IV has good power but elevated false positives")
        print("        - Detects most drivers")
        print("        - But may incorrectly flag some passengers")
        analysis_results["conclusion"] = "GOOD_POWER_HIGH_FPR"
    else:
        print("\nResult: IV struggles to discriminate")
        print("        - May need larger sample size or stronger effects")
        analysis_results["conclusion"] = "INCONCLUSIVE"

    # Save results
    results_file = output_dir / "driver_passenger_results.json"
    with open(results_file, "w") as f:
        json.dump(analysis_results, f, indent=2, default=str)

    print(f"\nResults saved to: {results_file}")

    return analysis_results


def main():
    parser = argparse.ArgumentParser(description="Driver vs passenger ecDNA scenario")
    parser.add_argument("--config", type=Path, default=Path("causanta/simulate/params/default.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/driver_passenger"))
    parser.add_argument("--n-replicates", type=int, default=10)

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    run_driver_passenger_study(args.config, args.output_dir, args.n_replicates)


if __name__ == "__main__":
    main()
