#!/usr/bin/env python3
"""Example: Running causal analysis on simulation output.

This script demonstrates how to use the causanta.analyze module
to perform comprehensive causal inference on simulation data.

Usage:
    python examples/analyze_simulation.py output/run_YYYYMMDD_HHMMSS/
"""

from __future__ import annotations

import sys
from pathlib import Path

from causanta.analyze import (
    # Data loading
    load_simulation_output,
    # IV estimation
    estimate_iv_effects,
    estimate_segregation_iv,
    # Causal discovery
    discover_causal_structure,
    compare_to_ground_truth,
    # Propensity score matching
    propensity_score_matching,
    inverse_propensity_weighting,
    # Regression adjustment
    linear_regression_adjustment,
    doubly_robust_estimation,
    # SEM
    mediation_analysis,
    full_sem,
    # Sensitivity
    rosenbaum_bounds,
    omitted_variable_bias,
    placebo_test,
    # Report
    generate_analysis_report,
)


def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_simulation.py <output_dir>")
        sys.exit(1)

    output_dir = Path(sys.argv[1])
    if not output_dir.exists():
        print(f"Error: Directory not found: {output_dir}")
        sys.exit(1)

    print(f"Loading simulation data from: {output_dir}")
    data = load_simulation_output(output_dir, load_all_timesteps=True)

    print(f"\nData Summary:")
    print(f"  Timesteps: {data.n_timesteps}")
    print(f"  Tumor cells: {len(data.tumor_cells)}")
    print(f"  Division events: {data.n_divisions}")

    # 1. IV Estimation
    print("\n" + "=" * 60)
    print("1. Instrumental Variable Estimation")
    print("=" * 60)
    iv_results = estimate_iv_effects(data)
    for result in iv_results:
        print(f"\n  {result.outcome_name}:")
        print(f"    2SLS estimate: {result.second_stage_coef:.4f} ({result.second_stage_se:.4f})")
        print(f"    First-stage F: {result.first_stage_f_stat:.1f}")
        print(f"    Weak IV: {result.weak_instrument}")

    # Segregation as IV
    seg_iv = estimate_segregation_iv(data)
    if seg_iv:
        print(f"\n  Segregation IV:")
        print(f"    Effect: {seg_iv.second_stage_coef:.4f}")
    else:
        print(f"\n  Segregation IV: Insufficient lineage data")

    # 2. Causal Discovery
    print("\n" + "=" * 60)
    print("2. Causal Discovery")
    print("=" * 60)
    pc_dag = discover_causal_structure(data, algorithm="pc")
    print(f"\n  PC Algorithm:")
    print(f"    Edges discovered: {len(pc_dag.edges)}")
    for edge in pc_dag.edges[:5]:
        print(f"      {edge}")

    ges_dag = discover_causal_structure(data, algorithm="ges")
    print(f"\n  GES Algorithm:")
    print(f"    Edges: {len(ges_dag.edges)}")

    # 3. Propensity Score Matching
    print("\n" + "=" * 60)
    print("3. Propensity Score Matching")
    print("=" * 60)
    matching_result = propensity_score_matching(data)
    print(f"\n  Matching:")
    print(f"    N treated: {matching_result.n_treated}")
    print(f"    N matched: {matching_result.n_matched}")
    print(f"    ATT: {matching_result.att:.4f} ({matching_result.att_se:.4f})")

    ipw_result = inverse_propensity_weighting(data)
    print(f"\n  IPW:")
    print(f"    ATE: {ipw_result['ate']:.4f} ({ipw_result['ate_se']:.4f})")

    # 4. Regression Adjustment
    print("\n" + "=" * 60)
    print("4. Regression Adjustment")
    print("=" * 60)
    reg_result = linear_regression_adjustment(data)
    print(f"\n  Linear Regression:")
    print(f"    Effect: {reg_result.treatment_coef:.4f} ({reg_result.treatment_se:.4f})")
    print(f"    R-squared: {reg_result.r_squared:.3f}")

    dr_result = doubly_robust_estimation(data)
    print(f"\n  Doubly Robust:")
    print(f"    ATE: {dr_result['ate']:.4f} ({dr_result['ate_se']:.4f})")

    # 5. Structural Equation Modeling
    print("\n" + "=" * 60)
    print("5. Structural Equation Modeling")
    print("=" * 60)
    mediation = mediation_analysis(data)
    print(f"\n  Mediation Analysis:")
    print(f"    Direct effect: {list(mediation.direct_effects.values())[0] if mediation.direct_effects else 0:.4f}")
    print(f"    Indirect effect: {mediation.indirect_effects.get('total_indirect', 0):.4f}")

    sem_result = full_sem(data)
    print(f"\n  Full SEM:")
    print(f"    Paths estimated: {len(sem_result.paths)}")
    for path in sem_result.paths[:3]:
        print(f"      {path.from_var} -> {path.to_var}: {path.coefficient:.4f}")

    # 6. Sensitivity Analysis
    print("\n" + "=" * 60)
    print("6. Sensitivity Analysis")
    print("=" * 60)
    rosenbaum = rosenbaum_bounds(data)
    print(f"\n  Rosenbaum Bounds:")
    print(f"    Original effect: {rosenbaum.original_estimate:.4f}")
    print(f"    Breakdown point (gamma): {rosenbaum.breakdown_point:.2f}")

    ovb = omitted_variable_bias(data)
    print(f"\n  Omitted Variable Bias:")
    print(f"    Robustness value: {ovb.robustness_value:.3f}")

    placebo = placebo_test(data)
    print(f"\n  Placebo Test:")
    print(f"    {placebo.get('interpretation', 'N/A')}")

    # Generate full HTML report
    print("\n" + "=" * 60)
    print("Generating HTML Report...")
    print("=" * 60)
    report_path = generate_analysis_report(
        data=data,
        output_path=output_dir / "analysis_report.html",
        methods=["iv", "discovery", "matching", "regression", "sem", "sensitivity"],
    )
    print(f"\nReport generated: {report_path}")


if __name__ == "__main__":
    main()
