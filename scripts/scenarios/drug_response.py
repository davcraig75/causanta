#!/usr/bin/env python3
"""Biological Scenario: Drug Response Heterogeneity

Scientific Question:
    Can we predict drug response heterogeneity from ecDNA copy number?
    Does ecDNA-driven EGFR overexpression predict response to EGFR inhibitors?

Hypothesis:
    If ecDNA causally drives EGFR-dependent phenotypes:
    - Cells with more ecDNA should be more dependent on EGFR signaling
    - Simulated EGFR inhibition should have larger effect on high-ecDNA cells
    - IV can estimate the dose-response relationship

Simulation Design:
    1. Simulate baseline tumor growth
    2. Apply simulated drug effect (reduce EGFR effect parameters)
    3. Compare tumor composition before/after treatment
    4. Use IV to estimate ecDNA → drug sensitivity relationship

Usage:
    python scripts/scenarios/drug_response.py --output-dir results/drug_response
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np


def simulate_drug_effect(
    base_params: dict,
    drug_efficacy: float,
) -> dict:
    """Modify parameters to simulate EGFR inhibitor effect.

    Args:
        base_params: Original simulation parameters
        drug_efficacy: Fraction of EGFR effect blocked (0.0 to 1.0)

    Returns:
        Modified parameters with reduced EGFR effects
    """
    import copy
    params = copy.deepcopy(base_params)

    # Reduce ecDNA effect parameters proportionally
    tumor_params = params["cell_types"][6]

    for key in ["ecDNA_effect_on_division", "ecDNA_effect_on_VEGF",
                "ecDNA_effect_on_migration", "ecDNA_effect_on_survival"]:
        if key in tumor_params:
            tumor_params[key] *= (1 - drug_efficacy)

    return params


def run_drug_response_study(
    base_config: Path,
    output_dir: Path,
    n_replicates: int = 10,
) -> dict:
    """Run drug response heterogeneity study.

    Simulates tumors under varying drug conditions and analyzes
    the relationship between ecDNA and drug sensitivity.
    """
    from causanta.simulate.sweep import SweepConfig, run_sweep

    print("=" * 60)
    print("BIOLOGICAL SCENARIO: Drug Response Heterogeneity")
    print("=" * 60)
    print("\nQuestion: Does ecDNA predict response to EGFR inhibitors?")
    print("Hypothesis: High-ecDNA cells should be more drug-sensitive")

    # Load base config
    with open(base_config) as f:
        base_params = json.load(f)

    # Create configs for different drug doses
    drug_doses = [0.0, 0.25, 0.50, 0.75, 1.0]  # Fraction of effect blocked

    all_results = []

    for dose in drug_doses:
        print(f"\n{'='*40}")
        print(f"Drug efficacy: {dose*100:.0f}%")
        print("=" * 40)

        # Create modified config
        modified_params = simulate_drug_effect(base_params, dose)

        # Save to temp file
        dose_config = output_dir / f"params_dose{int(dose*100)}.json"
        with open(dose_config, "w") as f:
            json.dump(modified_params, f, indent=2)

        # Run sweep
        sweep_config = SweepConfig(
            name=f"drug_dose{int(dose*100)}",
            base_config_path=dose_config,
            output_dir=output_dir / f"dose{int(dose*100)}",
            n_replicates=n_replicates,
        )

        def progress(i, total, result):
            status = "OK" if result["status"] == "success" else "FAIL"
            print(f"  [{i}/{total}] {result['name']}: {status}")

        results = run_sweep(sweep_config, n_workers=4, progress_callback=progress)

        # Analyze results
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

            if len(tumor_cells) < 30:
                continue

            # Compute summary statistics
            ecdna_counts = [c["ecDNA_count"] for c in tumor_cells]

            all_results.append({
                "drug_dose": dose,
                "n_tumor_cells": len(tumor_cells),
                "mean_ecdna": np.mean(ecdna_counts),
                "std_ecdna": np.std(ecdna_counts),
                "ecdna_range": [min(ecdna_counts), max(ecdna_counts)],
            })

    # Analysis
    print("\n" + "=" * 60)
    print("ANALYSIS: ecDNA vs Drug Response")
    print("=" * 60)

    analysis_results = {
        "scenario": "drug_response",
        "question": "Does ecDNA predict EGFR inhibitor response?",
        "conditions": all_results,
    }

    # Summarize by drug dose
    print("\nTumor cell counts by drug dose:")
    for dose in drug_doses:
        dose_results = [r for r in all_results if r["drug_dose"] == dose]
        if dose_results:
            mean_cells = np.mean([r["n_tumor_cells"] for r in dose_results])
            mean_ecdna = np.mean([r["mean_ecdna"] for r in dose_results])
            print(f"  {dose*100:3.0f}% efficacy: {mean_cells:6.0f} cells, mean ecDNA = {mean_ecdna:.1f}")

    # Test dose-response relationship
    if len(all_results) > 3:
        doses = [r["drug_dose"] for r in all_results]
        cells = [r["n_tumor_cells"] for r in all_results]

        # Simple correlation
        from scipy import stats
        corr, pval = stats.spearmanr(doses, cells)

        print(f"\nDose-response correlation: r = {corr:.3f}, p = {pval:.4f}")

        if corr < -0.3 and pval < 0.05:
            print("Result: Significant negative dose-response relationship detected")
            print("        Higher drug doses reduce tumor cell counts")
            analysis_results["dose_response"] = "SIGNIFICANT"
        else:
            print("Result: No significant dose-response relationship")
            analysis_results["dose_response"] = "NOT_SIGNIFICANT"

    # ecDNA stratification analysis
    print("\n" + "-" * 40)
    print("ecDNA Stratification Analysis")
    print("-" * 40)

    # Use untreated condition for ecDNA analysis
    untreated_results = [r for r in all_results if r["drug_dose"] == 0.0]
    if untreated_results:
        print(f"\nUntreated condition:")
        print(f"  Mean ecDNA: {np.mean([r['mean_ecdna'] for r in untreated_results]):.1f}")

        # Compare to fully treated
        treated_results = [r for r in all_results if r["drug_dose"] == 1.0]
        if treated_results:
            print(f"\nFully treated condition:")
            print(f"  Mean ecDNA: {np.mean([r['mean_ecdna'] for r in treated_results]):.1f}")

            # Selection effect
            untreated_ecdna = np.mean([r["mean_ecdna"] for r in untreated_results])
            treated_ecdna = np.mean([r["mean_ecdna"] for r in treated_results])

            if treated_ecdna < untreated_ecdna:
                print("\nResult: Drug treatment selects against high-ecDNA cells")
                print("        (mean ecDNA decreased under treatment)")
                analysis_results["selection"] = "AGAINST_HIGH_ECDNA"
            else:
                print("\nResult: Drug treatment selects for high-ecDNA cells")
                print("        (high-ecDNA cells may have acquired resistance)")
                analysis_results["selection"] = "FOR_HIGH_ECDNA"

    # Save results
    results_file = output_dir / "drug_response_results.json"
    with open(results_file, "w") as f:
        json.dump(analysis_results, f, indent=2, default=str)

    print(f"\nResults saved to: {results_file}")

    return analysis_results


def main():
    parser = argparse.ArgumentParser(description="Drug response heterogeneity scenario")
    parser.add_argument("--config", type=Path, default=Path("causanta/simulate/params/default.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/drug_response"))
    parser.add_argument("--n-replicates", type=int, default=5)

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    run_drug_response_study(args.config, args.output_dir, args.n_replicates)


if __name__ == "__main__":
    main()
