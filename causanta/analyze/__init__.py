"""CAUSANTA Causal Analysis Framework.

Provides comprehensive causal inference methods for analyzing
simulation output or real spatial tissue data:

- Instrumental Variable (IV) / 2SLS estimation
- Causal discovery algorithms (PC, GES)
- Propensity score matching
- Regression adjustment with covariate control
- Structural Equation Modeling (SEM)
- Sensitivity analysis for unmeasured confounding
"""

# Data loading (single source of truth)
from .loader import (
    SimulationData,
    load_simulation_output,
    load_cells_tsv,
    load_lineage_tsv,
    load_summary_log,
    cells_to_array,
)

# Ground truth effect estimation
from .effects import (
    EffectEstimate,
    CausalEffectAnalysis,
    run_causal_analysis,
    analyze_ecDNA_segregation,
)

# IV estimation
from .iv import (
    IVEstimate,
    two_stage_least_squares,
    estimate_iv_effects,
    estimate_segregation_iv,
)

# Causal discovery
from .discovery import (
    DiscoveredDAG,
    pc_algorithm,
    ges_algorithm,
    discover_causal_structure,
    compare_to_ground_truth,
)

# Propensity score matching
from .matching import (
    MatchingResult,
    propensity_score_matching,
    inverse_propensity_weighting,
)

# Regression adjustment
from .regression import (
    RegressionResult,
    linear_regression_adjustment,
    polynomial_regression,
    interaction_regression,
    doubly_robust_estimation,
)

# Structural equation modeling
from .sem import (
    PathCoefficient,
    SEMResult,
    mediation_analysis,
    full_sem,
    sobel_test,
)

# Sensitivity analysis
from .sensitivity import (
    SensitivityResult,
    rosenbaum_bounds,
    e_value_analysis,
    omitted_variable_bias,
    placebo_test,
)

# Report generation
from .report import generate_analysis_report

__all__ = [
    # Data loading
    "SimulationData",
    "load_simulation_output",
    "load_cells_tsv",
    "load_lineage_tsv",
    "load_summary_log",
    "cells_to_array",
    # Effects
    "EffectEstimate",
    "CausalEffectAnalysis",
    "run_causal_analysis",
    "analyze_ecDNA_segregation",
    # IV
    "IVEstimate",
    "two_stage_least_squares",
    "estimate_iv_effects",
    "estimate_segregation_iv",
    # Discovery
    "DiscoveredDAG",
    "pc_algorithm",
    "ges_algorithm",
    "discover_causal_structure",
    "compare_to_ground_truth",
    # Matching
    "MatchingResult",
    "propensity_score_matching",
    "inverse_propensity_weighting",
    # Regression
    "RegressionResult",
    "linear_regression_adjustment",
    "polynomial_regression",
    "interaction_regression",
    "doubly_robust_estimation",
    # SEM
    "PathCoefficient",
    "SEMResult",
    "mediation_analysis",
    "full_sem",
    "sobel_test",
    # Sensitivity
    "SensitivityResult",
    "rosenbaum_bounds",
    "e_value_analysis",
    "omitted_variable_bias",
    "placebo_test",
    # Report
    "generate_analysis_report",
]
