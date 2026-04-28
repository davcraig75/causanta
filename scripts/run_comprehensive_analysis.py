#!/usr/bin/env python3
"""Comprehensive analysis runner for CAUSANTA simulations.

Runs all analysis methods on simulation output and generates
a complete results package suitable for publication.

Usage:
    python scripts/run_comprehensive_analysis.py output/run_*/data/
    python scripts/run_comprehensive_analysis.py output/run_*/data/ --output results/
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from causanta.analyze.loader import SimulationData, load_simulation_output
from causanta.analyze.effects import run_causal_analysis, analyze_ecDNA_segregation
from causanta.analyze.iv import estimate_iv_effects, estimate_segregation_iv
from causanta.analyze.discovery import discover_causal_structure, compare_to_ground_truth
from causanta.analyze.sensitivity import rosenbaum_bounds, e_value_analysis, omitted_variable_bias, placebo_test
from causanta.analyze.regression import ols_with_inference, polynomial_regression
from causanta.analyze.matching import propensity_score_matching
from causanta.analyze.bootstrap import bootstrap_causal_effects, bootstrap_segregation_test, bootstrap_effect_comparison
from causanta.analyze.power import analyze_power_from_data, comprehensive_power_analysis
from causanta.analyze.heterogeneity import (
    stratified_iv_by_region,
    stratified_iv_by_hypoxia,
    stratified_iv_by_ecDNA_level,
    radial_effect_profile,
    comprehensive_heterogeneity_analysis,
)


def run_segregation_analysis(data: SimulationData) -> dict[str, Any]:
    """Analyze ecDNA segregation statistics."""
    print("  Running segregation analysis...")

    lineage = data.lineage
    if not lineage:
        return {"error": "No lineage data"}

    # Basic segregation stats
    seg_stats = analyze_ecDNA_segregation(lineage)

    # Compute daughter fractions for bootstrap test
    daughter_fractions = []
    for rec in lineage:
        total = rec["parent_ecDNA_after"] + rec["daughter_ecDNA"]
        if total > 0:
            daughter_fractions.append(rec["daughter_ecDNA"] / total)

    # Bootstrap test
    if len(daughter_fractions) >= 20:
        boot_test = bootstrap_segregation_test(
            np.array(daughter_fractions),
            null_value=0.5,
            n_bootstrap=1000,
        )
    else:
        boot_test = {"error": "Insufficient data"}

    return {
        "basic_stats": seg_stats,
        "bootstrap_test": boot_test,
        "n_divisions": len(lineage),
    }


def run_first_stage_analysis(data: SimulationData) -> dict[str, Any]:
    """Analyze first-stage regression (ecDNA → EGFR)."""
    print("  Running first-stage analysis...")

    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 30:
        return {"error": "Insufficient tumor cells"}

    Z = np.array([c["ecDNA_count"] for c in tumor_cells], dtype=float)
    D = np.array([c.get("egfr_expression", 0) for c in tumor_cells], dtype=float)

    if D.sum() == 0:
        D = Z  # Fallback if EGFR not tracked

    n = len(Z)
    X = np.column_stack([np.ones(n), Z])

    try:
        beta = np.linalg.lstsq(X, D, rcond=None)[0]
        D_hat = X @ beta

        # R-squared
        ss_res = np.sum((D - D_hat) ** 2)
        ss_tot = np.sum((D - np.mean(D)) ** 2)
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        # F-statistic
        k = 1  # number of instruments
        f_stat = (r_squared / k) / ((1 - r_squared) / (n - k - 1)) if r_squared < 1 else np.inf

        # Standard error
        mse = ss_res / (n - 2)
        XtX_inv = np.linalg.pinv(X.T @ X)
        se = np.sqrt(mse * XtX_inv[1, 1])

        return {
            "intercept": float(beta[0]),
            "coefficient": float(beta[1]),
            "se": float(se),
            "r_squared": float(r_squared),
            "f_statistic": float(f_stat),
            "n_observations": n,
            "weak_instrument": f_stat < 10,
            "interpretation": (
                f"Each ecDNA copy adds {beta[1]:.3f} units of EGFR expression. "
                f"F-statistic = {f_stat:.1f} {'(strong instrument)' if f_stat > 10 else '(weak instrument!)'}"
            ),
        }
    except Exception as e:
        return {"error": str(e)}


def run_iv_analysis(data: SimulationData) -> dict[str, Any]:
    """Run IV/2SLS estimation for all outcomes."""
    print("  Running IV analysis...")

    outcomes = ["VEGF_secretion", "migration_rate"]
    iv_results = estimate_iv_effects(data, outcomes)

    # Also run segregation IV
    seg_iv = estimate_segregation_iv(data)

    # Bootstrap CIs
    tumor_cells = data.tumor_cells
    if len(tumor_cells) >= 30:
        boot_results = bootstrap_causal_effects(
            tumor_cells, outcomes, n_bootstrap=500
        )
    else:
        boot_results = {}

    # Effect comparison (OLS vs IV)
    comparisons = {}
    if len(tumor_cells) >= 30:
        Z = np.array([c["ecDNA_count"] for c in tumor_cells], dtype=float)
        D = np.array([c.get("egfr_expression", c["ecDNA_count"]) for c in tumor_cells], dtype=float)
        O2 = np.array([c["O2_local"] for c in tumor_cells], dtype=float)
        glucose = np.array([c["glucose_local"] for c in tumor_cells], dtype=float)
        covariates = np.column_stack([O2, glucose])

        for outcome in outcomes:
            if outcome in tumor_cells[0]:
                Y = np.array([c[outcome] for c in tumor_cells], dtype=float)
                comp = bootstrap_effect_comparison(Z, D, Y, covariates, n_bootstrap=500)
                comparisons[outcome] = comp

    return {
        "iv_estimates": [r.to_dict() for r in iv_results],
        "segregation_iv": seg_iv.to_dict() if seg_iv else None,
        "bootstrap_cis": {k: v.to_dict() if hasattr(v, 'to_dict') else v
                         for k, v in boot_results.items()},
        "ols_iv_comparison": comparisons,
    }


def run_causal_discovery(data: SimulationData) -> dict[str, Any]:
    """Run causal discovery algorithms."""
    print("  Running causal discovery...")

    # PC algorithm
    pc_result = discover_causal_structure(data, algorithm="pc", alpha=0.05)

    # GES algorithm
    ges_result = discover_causal_structure(data, algorithm="ges")

    # Compare to ground truth
    ground_truth_edges = [
        ("ecDNA_count", "egfr_expression"),
        ("egfr_expression", "VEGF_secretion"),
        ("egfr_expression", "migration_rate"),
        ("O2_local", "VEGF_secretion"),
        ("O2_local", "migration_rate"),
    ]

    pc_comparison = compare_to_ground_truth(pc_result, ground_truth_edges)
    ges_comparison = compare_to_ground_truth(ges_result, ground_truth_edges)

    return {
        "pc_algorithm": pc_result.to_dict(),
        "ges_algorithm": ges_result.to_dict(),
        "pc_vs_truth": pc_comparison,
        "ges_vs_truth": ges_comparison,
    }


def run_sensitivity_analysis(data: SimulationData) -> dict[str, Any]:
    """Run sensitivity analyses."""
    print("  Running sensitivity analysis...")

    results = {}

    # Rosenbaum bounds
    for outcome in ["VEGF_secretion", "migration_rate"]:
        try:
            rb = rosenbaum_bounds(data, outcome_name=outcome)
            results[f"rosenbaum_{outcome}"] = rb.to_dict()
        except Exception as e:
            results[f"rosenbaum_{outcome}"] = {"error": str(e)}

    # E-values (need an IV estimate first)
    iv_results = estimate_iv_effects(data, ["migration_rate"])
    if iv_results:
        iv_est = iv_results[0]
        ev = e_value_analysis(iv_est.second_stage_coef, iv_est.second_stage_se)
        results["e_value_migration"] = ev.to_dict()

    # Omitted variable bias
    try:
        ovb = omitted_variable_bias(data)
        results["omitted_variable_bias"] = ovb.to_dict()
    except Exception as e:
        results["omitted_variable_bias"] = {"error": str(e)}

    # Placebo tests
    try:
        placebo = placebo_test(data)
        results["placebo_tests"] = placebo
    except Exception as e:
        results["placebo_tests"] = {"error": str(e)}

    return results


def run_heterogeneity_analysis(data: SimulationData) -> dict[str, Any]:
    """Run effect heterogeneity analyses."""
    print("  Running heterogeneity analysis...")

    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 50:
        return {"error": "Insufficient tumor cells"}

    results = comprehensive_heterogeneity_analysis(tumor_cells)
    return results


def run_power_analysis(data: SimulationData) -> dict[str, Any]:
    """Run power analysis."""
    print("  Running power analysis...")

    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 30:
        return {"error": "Insufficient tumor cells"}

    # Data-based power analysis
    power_result = analyze_power_from_data(tumor_cells)

    # Comprehensive power curves
    power_curves = comprehensive_power_analysis(
        effect_sizes=[0.01, 0.02, 0.05, 0.1],
        sample_sizes=[100, 200, 500, 1000, 2000, 5000],
        first_stage_r2=power_result.first_stage_r2,
    )

    return {
        "from_data": power_result.to_dict(),
        "power_curves": power_curves,
    }


def run_ground_truth_comparison(data: SimulationData, config: dict) -> dict[str, Any]:
    """Compare estimates to ground truth from config."""
    print("  Running ground truth comparison...")

    # Extract ground truth from config
    tumor_config = config.get("cell_types", {}).get("6", {})
    ground_truth = {
        "alpha": tumor_config.get("ecDNA_effect_on_division", 0.3),
        "beta": tumor_config.get("ecDNA_effect_on_VEGF", 0.1),
        "delta": tumor_config.get("ecDNA_effect_on_migration", 0.05),
        "gamma": tumor_config.get("ecDNA_effect_on_survival", 0.5),
    }

    # Get estimates
    iv_results = estimate_iv_effects(data, ["VEGF_secretion", "migration_rate"])

    comparisons = {}
    for iv_est in iv_results:
        if iv_est.outcome_name == "VEGF_secretion":
            true_val = ground_truth["beta"]
            param = "beta"
        elif iv_est.outcome_name == "migration_rate":
            true_val = ground_truth["delta"]
            param = "delta"
        else:
            continue

        iv_bias = (iv_est.second_stage_coef - true_val) / true_val if true_val != 0 else 0
        ols_bias = (iv_est.ols_coef - true_val) / true_val if true_val != 0 else 0

        comparisons[param] = {
            "ground_truth": true_val,
            "iv_estimate": iv_est.second_stage_coef,
            "iv_se": iv_est.second_stage_se,
            "iv_relative_bias": iv_bias,
            "ols_estimate": iv_est.ols_coef,
            "ols_relative_bias": ols_bias,
            "iv_within_2se": abs(iv_est.second_stage_coef - true_val) <= 2 * iv_est.second_stage_se,
        }

    return {
        "ground_truth": ground_truth,
        "comparisons": comparisons,
    }


def generate_summary(results: dict[str, Any]) -> str:
    """Generate human-readable summary of results."""
    lines = [
        "=" * 70,
        "CAUSANTA Comprehensive Analysis Summary",
        "=" * 70,
        "",
    ]

    # Segregation
    if "segregation" in results:
        seg = results["segregation"]
        if "basic_stats" in seg:
            stats = seg["basic_stats"]
            lines.append("ecDNA SEGREGATION:")
            lines.append(f"  Divisions analyzed: {stats.get('n_divisions', 'N/A')}")
            lines.append(f"  Mean daughter fraction: {stats.get('mean_daughter_fraction', 'N/A'):.3f} (expected: 0.5)")
            lines.append(f"  Variance: {stats.get('segregation_variance', 'N/A'):.4f}")
            lines.append("")

    # First stage
    if "first_stage" in results:
        fs = results["first_stage"]
        if "f_statistic" in fs:
            lines.append("FIRST-STAGE REGRESSION (ecDNA → EGFR):")
            lines.append(f"  F-statistic: {fs['f_statistic']:.1f} {'(STRONG)' if fs['f_statistic'] > 10 else '(WEAK!)'}")
            lines.append(f"  R²: {fs['r_squared']:.3f}")
            lines.append(f"  Coefficient: {fs['coefficient']:.3f} (SE: {fs['se']:.3f})")
            lines.append("")

    # IV estimates
    if "iv_analysis" in results:
        iv = results["iv_analysis"]
        if "iv_estimates" in iv:
            lines.append("IV ESTIMATES:")
            for est in iv["iv_estimates"]:
                lines.append(f"  {est['outcome']}:")
                lines.append(f"    IV:  {est['second_stage_coef']:.4f} (SE: {est['second_stage_se']:.4f})")
                lines.append(f"    OLS: {est['ols_coef']:.4f}")
                lines.append(f"    Bias ratio (OLS/IV): {est['bias_ratio']:.2f}")
            lines.append("")

    # Ground truth comparison
    if "ground_truth_comparison" in results:
        gt = results["ground_truth_comparison"]
        if "comparisons" in gt:
            lines.append("GROUND TRUTH COMPARISON:")
            for param, comp in gt["comparisons"].items():
                lines.append(f"  {param}:")
                lines.append(f"    True: {comp['ground_truth']:.4f}")
                lines.append(f"    IV:   {comp['iv_estimate']:.4f} (bias: {comp['iv_relative_bias']*100:.1f}%)")
                lines.append(f"    OLS:  {comp['ols_estimate']:.4f} (bias: {comp['ols_relative_bias']*100:.1f}%)")
                lines.append(f"    Within 2SE: {'Yes' if comp['iv_within_2se'] else 'No'}")
            lines.append("")

    # Sensitivity
    if "sensitivity" in results:
        sens = results["sensitivity"]
        lines.append("SENSITIVITY ANALYSIS:")
        for key, val in sens.items():
            if isinstance(val, dict) and "robustness_value" in val:
                lines.append(f"  {key}: Γ* = {val.get('breakdown_point', 'N/A'):.2f}")
        lines.append("")

    # Power
    if "power_analysis" in results:
        pa = results["power_analysis"]
        if "from_data" in pa:
            fd = pa["from_data"]
            lines.append("POWER ANALYSIS:")
            lines.append(f"  First-stage R²: {fd.get('first_stage_r2', 'N/A'):.3f}")
            lines.append(f"  Current power (δ=0.05): {fd.get('achieved_power', 'N/A'):.2%}")
            lines.append(f"  Required N for 80% power: {fd.get('required_n', 'N/A')}")
            lines.append("")

    lines.append("=" * 70)
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Run comprehensive CAUSANTA analysis"
    )
    parser.add_argument(
        "data_dir",
        type=Path,
        help="Path to simulation data directory",
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        default=None,
        help="Output directory for results (default: data_dir/../analysis)",
    )
    parser.add_argument(
        "--skip-bootstrap",
        action="store_true",
        help="Skip bootstrap analyses (faster)",
    )
    args = parser.parse_args()

    data_dir = args.data_dir
    if not data_dir.exists():
        print(f"Error: Data directory not found: {data_dir}")
        sys.exit(1)

    # Determine output directory
    if args.output:
        output_dir = args.output
    else:
        output_dir = data_dir.parent / "analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"CAUSANTA Comprehensive Analysis")
    print(f"  Data: {data_dir}")
    print(f"  Output: {output_dir}")
    print()

    # Load data
    print("Loading simulation data...")
    data = load_simulation_output(data_dir)

    # Load config for ground truth
    config_path = data_dir / "params.json"
    if config_path.exists():
        with open(config_path) as f:
            config = json.load(f)
    else:
        config = {}

    # Run all analyses
    results = {}

    print("\nRunning analyses:")

    results["segregation"] = run_segregation_analysis(data)
    results["first_stage"] = run_first_stage_analysis(data)
    results["iv_analysis"] = run_iv_analysis(data)
    results["causal_discovery"] = run_causal_discovery(data)
    results["sensitivity"] = run_sensitivity_analysis(data)
    results["heterogeneity"] = run_heterogeneity_analysis(data)
    results["power_analysis"] = run_power_analysis(data)
    results["ground_truth_comparison"] = run_ground_truth_comparison(data, config)

    # Generate summary
    summary = generate_summary(results)
    print("\n" + summary)

    # Save results
    results_path = output_dir / "comprehensive_results.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults saved to: {results_path}")

    summary_path = output_dir / "analysis_summary.txt"
    with open(summary_path, "w") as f:
        f.write(summary)
    print(f"Summary saved to: {summary_path}")


if __name__ == "__main__":
    main()
