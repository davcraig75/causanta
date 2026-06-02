#!/usr/bin/env python3
"""Simulation study runner for CAUSANTA validation.

Orchestrates:
1. Baseline validation (multiple replicates)
2. Effect size sweep
3. Sample size study
4. Confounding strength study
5. Comprehensive analysis of all results
6. Figure generation

Usage:
    python scripts/run_simulation_study.py --study baseline --n-replicates 20
    python scripts/run_simulation_study.py --study full --output-dir results/full_study
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np


def run_baseline_study(
    base_config: Path,
    output_dir: Path,
    n_replicates: int = 20,
    n_workers: int = 4,
) -> dict:
    """Run baseline validation study.

    Multiple replicates with default parameters to establish
    baseline behavior and variance estimates.
    """
    from causanta.simulate.sweep import baseline_validation_sweep, run_sweep

    print(f"\n{'='*60}")
    print("BASELINE VALIDATION STUDY")
    print(f"{'='*60}")
    print(f"Replicates: {n_replicates}")
    print(f"Output: {output_dir}")

    sweep_config = baseline_validation_sweep(
        base_config=base_config,
        output_dir=output_dir / "baseline",
        n_replicates=n_replicates,
    )

    def progress(i, total, result):
        status = "OK" if result["status"] == "success" else "FAIL"
        print(f"  [{i}/{total}] {result['name']}: {status}")

    results = run_sweep(sweep_config, n_workers=n_workers, progress_callback=progress)

    success = sum(1 for r in results if r["status"] == "success")
    print(f"\nCompleted: {success}/{len(results)} successful")

    return {
        "study": "baseline",
        "n_replicates": n_replicates,
        "success_rate": success / len(results),
        "results": results,
    }


def run_effect_size_study(
    base_config: Path,
    output_dir: Path,
    n_replicates: int = 10,
    n_workers: int = 4,
) -> dict:
    """Run effect size sweep study.

    Varies causal effect parameters to study detection power
    across different effect magnitudes.
    """
    from causanta.simulate.sweep import SweepConfig, run_sweep

    print(f"\n{'='*60}")
    print("EFFECT SIZE SWEEP STUDY")
    print(f"{'='*60}")

    # Define effect size variations
    # Focus on each parameter individually to avoid combinatorial explosion
    effect_variations = {
        "ecDNA_effect_on_division": [0.15, 0.30, 0.45],
        "ecDNA_effect_on_VEGF": [0.05, 0.10, 0.15],
        "ecDNA_effect_on_migration": [0.025, 0.05, 0.075],
        "ecDNA_effect_on_survival": [0.10, 0.20, 0.30],
    }

    all_results = []

    for param_name, values in effect_variations.items():
        print(f"\n  Varying: {param_name}")

        sweep_config = SweepConfig(
            name=f"effect_{param_name}",
            base_config_path=base_config,
            output_dir=output_dir / f"effect_{param_name}",
            n_replicates=n_replicates,
            vary_params={f"cell_types.6.{param_name}": values},
        )

        def progress(i, total, result):
            status = "OK" if result["status"] == "success" else "FAIL"
            print(f"    [{i}/{total}] {result['name']}: {status}")

        results = run_sweep(sweep_config, n_workers=n_workers, progress_callback=progress)
        all_results.extend(results)

    success = sum(1 for r in all_results if r["status"] == "success")
    print(f"\nTotal completed: {success}/{len(all_results)} successful")

    return {
        "study": "effect_size",
        "n_replicates": n_replicates,
        "variations": effect_variations,
        "success_rate": success / len(all_results),
        "results": all_results,
    }


def run_sample_size_study(
    base_config: Path,
    output_dir: Path,
    n_replicates: int = 10,
    n_workers: int = 4,
) -> dict:
    """Run sample size study.

    Varies simulation duration to generate different sample sizes
    for power analysis.
    """
    from causanta.simulate.sweep import sample_size_sweep, run_sweep

    print(f"\n{'='*60}")
    print("SAMPLE SIZE STUDY")
    print(f"{'='*60}")

    sweep_config = sample_size_sweep(
        base_config=base_config,
        output_dir=output_dir / "sample_size",
        n_replicates=n_replicates,
    )

    def progress(i, total, result):
        status = "OK" if result["status"] == "success" else "FAIL"
        print(f"  [{i}/{total}] {result['name']}: {status}")

    results = run_sweep(sweep_config, n_workers=n_workers, progress_callback=progress)

    success = sum(1 for r in results if r["status"] == "success")
    print(f"\nCompleted: {success}/{len(results)} successful")

    return {
        "study": "sample_size",
        "n_replicates": n_replicates,
        "success_rate": success / len(results),
        "results": results,
    }


def run_confounding_study(
    base_config: Path,
    output_dir: Path,
    n_replicates: int = 10,
    n_workers: int = 4,
) -> dict:
    """Run confounding strength study.

    Varies hypoxia parameters to study IV performance under
    different confounding scenarios.
    """
    from causanta.simulate.sweep import confounding_strength_sweep, run_sweep

    print(f"\n{'='*60}")
    print("CONFOUNDING STRENGTH STUDY")
    print(f"{'='*60}")

    sweep_config = confounding_strength_sweep(
        base_config=base_config,
        output_dir=output_dir / "confounding",
        n_replicates=n_replicates,
    )

    def progress(i, total, result):
        status = "OK" if result["status"] == "success" else "FAIL"
        print(f"  [{i}/{total}] {result['name']}: {status}")

    results = run_sweep(sweep_config, n_workers=n_workers, progress_callback=progress)

    success = sum(1 for r in results if r["status"] == "success")
    print(f"\nCompleted: {success}/{len(results)} successful")

    return {
        "study": "confounding",
        "n_replicates": n_replicates,
        "success_rate": success / len(results),
        "results": results,
    }


def analyze_study_results(study_dir: Path) -> dict:
    """Run comprehensive analysis on study results."""
    from causanta.analyze.bootstrap import bootstrap_iv_estimate
    from causanta.analyze.power import comprehensive_power_analysis
    from causanta.analyze.heterogeneity import stratified_iv_by_region

    print(f"\n{'='*60}")
    print("ANALYZING RESULTS")
    print(f"{'='*60}")

    analysis_results = {
        "iv_estimates": [],
        "bootstrap_cis": [],
        "power_analyses": [],
        "heterogeneity": [],
    }

    # Find all successful run directories
    for sweep_results_file in study_dir.rglob("sweep_results.json"):
        with open(sweep_results_file) as f:
            results = json.load(f)

        for result in results:
            if result["status"] != "success":
                continue

            run_dir = Path(result["output_dir"])
            cells_file = run_dir / "cells_final.json"

            if not cells_file.exists():
                continue

            print(f"  Analyzing: {result['name']}")

            try:
                with open(cells_file) as f:
                    cells_data = json.load(f)

                # Filter tumor cells
                tumor_cells = [c for c in cells_data if c.get("cell_type") == 6]

                if len(tumor_cells) < 50:
                    continue

                # Extract data
                Z = np.array([c["ecDNA_count"] for c in tumor_cells], dtype=float)
                D = np.array([c.get("egfr_expression", c["ecDNA_count"]) for c in tumor_cells], dtype=float)
                Y = np.array([c.get("VEGF_secretion", 0) for c in tumor_cells], dtype=float)

                # Bootstrap IV
                bootstrap_result = bootstrap_iv_estimate(Z, D, Y, n_bootstrap=1000)

                analysis_results["iv_estimates"].append({
                    "name": result["name"],
                    "params": result.get("params", {}),
                    "iv_estimate": bootstrap_result.point_estimate,
                    "ci_lower": bootstrap_result.ci_lower,
                    "ci_upper": bootstrap_result.ci_upper,
                    "n_cells": len(tumor_cells),
                })

            except Exception as e:
                print(f"    Error: {e}")

    return analysis_results


def generate_figures(study_dir: Path, output_dir: Path) -> None:
    """Generate publication figures from study results."""
    from causanta.visualize import (
        apply_nature_style,
        create_main_figure_1,
        create_main_figure_2,
        create_main_figure_3,
        create_main_figure_4,
        create_main_figure_5,
        create_main_figure_6,
    )

    print(f"\n{'='*60}")
    print("GENERATING FIGURES")
    print(f"{'='*60}")

    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    apply_nature_style()

    # Figure 1: Framework
    print("  Figure 1: Conceptual framework...")
    create_main_figure_1(figures_dir)

    # Figure 2: Validation
    print("  Figure 2: Instrument validation...")
    create_main_figure_2(figures_dir)

    # Figure 3: Estimation
    print("  Figure 3: Causal estimation...")
    create_main_figure_3(figures_dir)

    # Figure 4: Power
    print("  Figure 4: Power analysis...")
    create_main_figure_4(figures_dir)

    # Figure 5: Heterogeneity
    print("  Figure 5: Effect heterogeneity...")
    create_main_figure_5(figures_dir)

    # Figure 6: Robustness
    print("  Figure 6: Robustness...")
    create_main_figure_6(figures_dir)

    print(f"\nFigures saved to: {figures_dir}")


def evaluate_success_criteria(analysis_results: dict) -> dict:
    """Evaluate results against publication success criteria."""
    print(f"\n{'='*60}")
    print("SUCCESS CRITERIA EVALUATION")
    print(f"{'='*60}")

    criteria = {}

    iv_estimates = analysis_results.get("iv_estimates", [])

    if iv_estimates:
        # Criterion 1: IV accuracy (<15% bias)
        true_effect = 0.10  # Ground truth beta
        biases = []
        for est in iv_estimates:
            bias = abs(est["iv_estimate"] - true_effect) / true_effect
            biases.append(bias)
        mean_bias = np.mean(biases)
        criteria["iv_accuracy"] = {
            "target": "< 15% bias",
            "achieved": f"{mean_bias*100:.1f}% bias",
            "passed": mean_bias < 0.15,
        }

        # Criterion 2: CI coverage (>= 90%)
        coverage_count = sum(
            1 for est in iv_estimates
            if est["ci_lower"] <= true_effect <= est["ci_upper"]
        )
        coverage = coverage_count / len(iv_estimates) if iv_estimates else 0
        criteria["ci_coverage"] = {
            "target": ">= 90%",
            "achieved": f"{coverage*100:.1f}%",
            "passed": coverage >= 0.90,
        }

    # Print results
    for name, result in criteria.items():
        status = "PASS" if result["passed"] else "FAIL"
        print(f"  {name}: {result['achieved']} (target: {result['target']}) [{status}]")

    return criteria


def main():
    parser = argparse.ArgumentParser(description="Run CAUSANTA simulation study")
    parser.add_argument("--study", choices=["baseline", "effect", "sample", "confound", "full"],
                        default="baseline", help="Study type to run")
    parser.add_argument("--config", type=Path, default=Path("causanta/simulate/params/default.json"),
                        help="Base configuration file")
    parser.add_argument("--output-dir", type=Path, default=Path("results/simulation_study"),
                        help="Output directory")
    parser.add_argument("--n-replicates", type=int, default=10,
                        help="Number of replicates per condition")
    parser.add_argument("--n-workers", type=int, default=4,
                        help="Number of parallel workers")
    parser.add_argument("--analyze-only", action="store_true",
                        help="Only run analysis on existing results")
    parser.add_argument("--figures-only", action="store_true",
                        help="Only generate figures from existing results")

    args = parser.parse_args()

    # Setup
    args.output_dir.mkdir(parents=True, exist_ok=True)

    start_time = time.time()
    print(f"\nCAUSANTA Simulation Study")
    print(f"Started: {datetime.now().isoformat()}")

    all_results = {}

    if not args.analyze_only and not args.figures_only:
        # Run simulations
        if args.study in ["baseline", "full"]:
            all_results["baseline"] = run_baseline_study(
                args.config, args.output_dir, args.n_replicates, args.n_workers
            )

        if args.study in ["effect", "full"]:
            all_results["effect_size"] = run_effect_size_study(
                args.config, args.output_dir, args.n_replicates, args.n_workers
            )

        if args.study in ["sample", "full"]:
            all_results["sample_size"] = run_sample_size_study(
                args.config, args.output_dir, args.n_replicates, args.n_workers
            )

        if args.study in ["confound", "full"]:
            all_results["confounding"] = run_confounding_study(
                args.config, args.output_dir, args.n_replicates, args.n_workers
            )

    if not args.figures_only:
        # Analyze results
        analysis_results = analyze_study_results(args.output_dir)

        # Save analysis
        analysis_file = args.output_dir / "analysis_results.json"
        with open(analysis_file, "w") as f:
            json.dump(analysis_results, f, indent=2, default=str)
        print(f"\nAnalysis saved to: {analysis_file}")

        # Evaluate success criteria
        criteria = evaluate_success_criteria(analysis_results)

        # Save criteria evaluation
        criteria_file = args.output_dir / "success_criteria.json"
        with open(criteria_file, "w") as f:
            json.dump(criteria, f, indent=2)

    # Generate figures
    generate_figures(args.output_dir, args.output_dir)

    # Summary
    elapsed = time.time() - start_time
    print(f"\n{'='*60}")
    print("STUDY COMPLETE")
    print(f"{'='*60}")
    print(f"Elapsed time: {elapsed/60:.1f} minutes")
    print(f"Output directory: {args.output_dir}")


if __name__ == "__main__":
    main()
