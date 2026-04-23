#!/usr/bin/env python3
"""Biological Scenario: Immune Selection on ecDNA

Scientific Question:
    Does immune pressure select for or against high-ecDNA tumor cells?
    What are the causal relationships between ecDNA, immune visibility,
    and survival under immune attack?

Hypotheses:
    H1: High ecDNA → More EGFR → More visible to immune system → Selected against
    H2: High ecDNA → Faster division → Can outpace immune killing → Selected for
    H3: ecDNA effect varies with immune pressure (interaction)

Simulation Design:
    1. Vary immune recruitment strength
    2. Compare ecDNA distribution in survivors vs all cells
    3. Use IV to estimate ecDNA → immune susceptibility relationship

Clinical Relevance:
    Understanding whether ecDNA-high tumors are immune hot or cold
    can inform immunotherapy patient selection.

Usage:
    python scripts/scenarios/immune_selection.py --output-dir results/immune_selection
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import numpy as np


def run_immune_selection_study(
    base_config: Path,
    output_dir: Path,
    n_replicates: int = 10,
) -> dict:
    """Run immune selection study.

    Tests how immune pressure affects ecDNA distribution in tumors
    and whether IV can identify causal relationships.
    """
    from causanta.simulate.sweep import SweepConfig, run_sweep
    from causanta.analyze.bootstrap import bootstrap_iv_estimate

    print("=" * 60)
    print("BIOLOGICAL SCENARIO: Immune Selection on ecDNA")
    print("=" * 60)
    print("\nQuestion: Does immune pressure select for or against high-ecDNA cells?")
    print("Key: Compare ecDNA distribution under varying immune strength")

    # Vary immune recruitment strength
    immune_strengths = [0.0, 0.001, 0.005, 0.01, 0.02]

    all_results = []

    for strength in immune_strengths:
        label = "NO_IMMUNE" if strength == 0 else f"IMMUNE_{strength}"

        print(f"\n{'='*40}")
        print(f"Testing: {label}")
        print("=" * 40)

        sweep_config = SweepConfig(
            name=f"immune_{strength}",
            base_config_path=base_config,
            output_dir=output_dir / f"immune_{strength}",
            n_replicates=n_replicates,
            fixed_params={
                "immune_recruitment.base_rate_per_tumor_per_hr": strength,
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
            lineage_file = run_dir / "lineage.json"

            if not cells_file.exists():
                continue

            with open(cells_file) as f:
                cells_data = json.load(f)

            tumor_cells = [c for c in cells_data if c.get("cell_type") == 6]
            immune_cells = [c for c in cells_data if c.get("cell_type") == 7]

            if len(tumor_cells) < 10:
                continue

            # Compute statistics
            ecdna_counts = [c["ecDNA_count"] for c in tumor_cells]

            # ecDNA distribution statistics
            mean_ecdna = np.mean(ecdna_counts)
            std_ecdna = np.std(ecdna_counts)
            high_ecdna_frac = np.mean([1 if e > 30 else 0 for e in ecdna_counts])

            all_results.append({
                "immune_strength": strength,
                "n_tumor_cells": len(tumor_cells),
                "n_immune_cells": len(immune_cells),
                "mean_ecdna": mean_ecdna,
                "std_ecdna": std_ecdna,
                "high_ecdna_fraction": high_ecdna_frac,
                "ecdna_range": [min(ecdna_counts), max(ecdna_counts)],
            })

    # Analysis
    print("\n" + "=" * 60)
    print("ANALYSIS: Immune Selection Effects")
    print("=" * 60)

    analysis_results = {
        "scenario": "immune_selection",
        "question": "Does immune pressure select for/against high-ecDNA cells?",
        "conditions": all_results,
    }

    # Compare ecDNA distribution across immune strengths
    print("\nTumor survival and ecDNA by immune strength:")
    print("-" * 60)

    for strength in immune_strengths:
        strength_results = [r for r in all_results if r["immune_strength"] == strength]
        if strength_results:
            mean_tumor = np.mean([r["n_tumor_cells"] for r in strength_results])
            mean_immune = np.mean([r["n_immune_cells"] for r in strength_results])
            mean_ecdna = np.mean([r["mean_ecdna"] for r in strength_results])
            high_frac = np.mean([r["high_ecdna_fraction"] for r in strength_results])

            label = "NONE" if strength == 0 else f"{strength:.3f}"
            print(f"  Immune={label}: {mean_tumor:6.0f} tumor, {mean_immune:5.0f} immune, "
                  f"mean_ecDNA={mean_ecdna:.1f}, high_frac={high_frac:.2f}")

    # Statistical test: does ecDNA change with immune pressure?
    print("\n" + "-" * 40)
    print("Selection Analysis")
    print("-" * 40)

    if len(all_results) > 3:
        from scipy import stats

        immune_vals = [r["immune_strength"] for r in all_results]
        ecdna_vals = [r["mean_ecdna"] for r in all_results]

        corr, pval = stats.spearmanr(immune_vals, ecdna_vals)

        print(f"\nCorrelation: immune strength vs mean ecDNA")
        print(f"  Spearman r = {corr:.3f}, p = {pval:.4f}")

        if pval < 0.05:
            if corr > 0:
                print("\nResult: Immune pressure SELECTS FOR high-ecDNA cells")
                print("        Interpretation: High-ecDNA cells escape immune killing")
                print("        (possibly through faster division outpacing killing)")
                analysis_results["selection_direction"] = "FOR_HIGH_ECDNA"
            else:
                print("\nResult: Immune pressure SELECTS AGAINST high-ecDNA cells")
                print("        Interpretation: High-ecDNA cells are more immune visible")
                print("        (possibly due to higher EGFR/neoantigen presentation)")
                analysis_results["selection_direction"] = "AGAINST_HIGH_ECDNA"
        else:
            print("\nResult: No significant selection effect detected")
            analysis_results["selection_direction"] = "NEUTRAL"

        analysis_results["correlation"] = {
            "spearman_r": corr,
            "p_value": pval,
        }

    # IV analysis for survival
    print("\n" + "-" * 40)
    print("IV Analysis: ecDNA → Survival")
    print("-" * 40)

    # Use high-immune condition for IV analysis
    high_immune_results = [r for r in all_results if r["immune_strength"] >= 0.01]
    if high_immune_results:
        # Find a run with enough cells
        for result in all_results:
            if result["immune_strength"] >= 0.01 and result["n_tumor_cells"] >= 50:
                run_dir = output_dir / f"immune_{result['immune_strength']}"

                # Find actual run directory
                for subdir in run_dir.iterdir():
                    cells_file = subdir / "cells_final.json"
                    if cells_file.exists():
                        with open(cells_file) as f:
                            cells_data = json.load(f)

                        tumor_cells = [c for c in cells_data if c.get("cell_type") == 6]

                        if len(tumor_cells) >= 50:
                            Z = np.array([c["ecDNA_count"] for c in tumor_cells], dtype=float)
                            D = np.array([c.get("egfr_expression", c["ecDNA_count"]) for c in tumor_cells], dtype=float)
                            # Use survival proxy (cells that exist survived)
                            Y = np.ones(len(tumor_cells))  # All survivors

                            # For meaningful analysis, we'd need death events
                            # Here we analyze proliferation rate as proxy
                            if "proliferation_rate" in tumor_cells[0]:
                                Y = np.array([c["proliferation_rate"] for c in tumor_cells], dtype=float)

                                try:
                                    bootstrap_result = bootstrap_iv_estimate(Z, D, Y, n_bootstrap=500)

                                    print(f"\nUnder immune pressure (strength={result['immune_strength']}):")
                                    print(f"  IV estimate:  {bootstrap_result.point_estimate:.4f}")
                                    print(f"  95% CI:       [{bootstrap_result.ci_lower:.4f}, {bootstrap_result.ci_upper:.4f}]")
                                    print(f"  OLS estimate: {bootstrap_result.ols_estimate:.4f}")

                                    analysis_results["iv_under_immune"] = {
                                        "iv_estimate": bootstrap_result.point_estimate,
                                        "ci_lower": bootstrap_result.ci_lower,
                                        "ci_upper": bootstrap_result.ci_upper,
                                        "ols_estimate": bootstrap_result.ols_estimate,
                                    }

                                except Exception as e:
                                    print(f"  IV analysis failed: {e}")

                            break
                    if "iv_under_immune" in analysis_results:
                        break
                break

    # Conclusions
    print("\n" + "=" * 60)
    print("CONCLUSIONS")
    print("=" * 60)

    selection = analysis_results.get("selection_direction", "UNKNOWN")

    if selection == "FOR_HIGH_ECDNA":
        print("\nImplications for immunotherapy:")
        print("  - High-ecDNA tumors may be resistant to immune checkpoint blockade")
        print("  - Consider combination with EGFR inhibitors to slow proliferation")
        print("  - ecDNA could be a biomarker for immunotherapy resistance")
    elif selection == "AGAINST_HIGH_ECDNA":
        print("\nImplications for immunotherapy:")
        print("  - High-ecDNA tumors may be 'hot' (immune visible)")
        print("  - Checkpoint blockade may be more effective in ecDNA-high tumors")
        print("  - ecDNA could be a biomarker for immunotherapy response")
    else:
        print("\nImplications:")
        print("  - No clear selection effect detected")
        print("  - ecDNA status may not predict immunotherapy response")
        print("  - Consider other biomarkers for patient selection")

    # Save results
    results_file = output_dir / "immune_selection_results.json"
    with open(results_file, "w") as f:
        json.dump(analysis_results, f, indent=2, default=str)

    print(f"\nResults saved to: {results_file}")

    return analysis_results


def main():
    parser = argparse.ArgumentParser(description="Immune selection scenario")
    parser.add_argument("--config", type=Path, default=Path("causanta/simulate/params/default.json"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/immune_selection"))
    parser.add_argument("--n-replicates", type=int, default=10)

    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    run_immune_selection_study(args.config, args.output_dir, args.n_replicates)


if __name__ == "__main__":
    main()
