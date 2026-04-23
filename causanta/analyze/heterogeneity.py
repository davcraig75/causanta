"""Effect heterogeneity analysis for spatial causal inference.

Implements methods to detect and quantify how causal effects vary:
- By tumor region (core, margin, infiltrating)
- Over simulation time (temporal dynamics)
- By microenvironmental context (hypoxia, nutrient status)
- By cell characteristics (generation, ecDNA level)

Key insight: Spatial causal inference can reveal WHERE effects are strongest,
guiding targeted therapeutic interventions.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import stats

from .iv import two_stage_least_squares


@dataclass
class StratifiedEffect:
    """Effect estimate for a specific stratum."""

    stratum_name: str
    n_observations: int
    iv_estimate: float
    iv_se: float
    iv_ci_lower: float
    iv_ci_upper: float
    ols_estimate: float
    ols_se: float
    first_stage_f: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "stratum": self.stratum_name,
            "n": self.n_observations,
            "iv_estimate": self.iv_estimate,
            "iv_se": self.iv_se,
            "iv_ci_lower": self.iv_ci_lower,
            "iv_ci_upper": self.iv_ci_upper,
            "ols_estimate": self.ols_estimate,
            "ols_se": self.ols_se,
            "first_stage_f": self.first_stage_f,
        }


@dataclass
class HeterogeneityAnalysis:
    """Result of heterogeneity analysis."""

    analysis_type: str
    outcome_name: str
    overall_effect: float
    stratified_effects: list[StratifiedEffect]
    heterogeneity_test: dict[str, float]
    interpretation: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "analysis_type": self.analysis_type,
            "outcome": self.outcome_name,
            "overall_effect": self.overall_effect,
            "stratified_effects": [s.to_dict() for s in self.stratified_effects],
            "heterogeneity_test": self.heterogeneity_test,
            "interpretation": self.interpretation,
        }


def _compute_tumor_center(cells: list[dict]) -> tuple[float, float]:
    """Compute centroid of tumor cells."""
    tumor_cells = [c for c in cells if c.get("cell_type") == 6]
    if not tumor_cells:
        return 0.0, 0.0
    x_mean = np.mean([c["x"] for c in tumor_cells])
    y_mean = np.mean([c["y"] for c in tumor_cells])
    return float(x_mean), float(y_mean)


def _compute_tumor_radius(cells: list[dict], center: tuple[float, float]) -> float:
    """Compute approximate tumor radius (95th percentile distance)."""
    tumor_cells = [c for c in cells if c.get("cell_type") == 6]
    if len(tumor_cells) < 5:
        return 100.0
    distances = [
        np.sqrt((c["x"] - center[0]) ** 2 + (c["y"] - center[1]) ** 2)
        for c in tumor_cells
    ]
    return float(np.percentile(distances, 95))


def _classify_tumor_region(
    cell: dict,
    center: tuple[float, float],
    radius: float,
    margin_width: float = 50.0,
) -> str:
    """Classify cell into tumor region.

    Regions:
    - core: > margin_width from tumor edge
    - margin: within margin_width of tumor edge
    - infiltrating: beyond tumor radius (if cell exists there)
    """
    dist = np.sqrt((cell["x"] - center[0]) ** 2 + (cell["y"] - center[1]) ** 2)

    if dist > radius:
        return "infiltrating"
    elif dist > radius - margin_width:
        return "margin"
    else:
        return "core"


def _estimate_iv_for_stratum(
    cells: list[dict],
    outcome_name: str,
    alpha: float = 0.05,
) -> StratifiedEffect | None:
    """Estimate IV effect for a stratum of cells."""
    if len(cells) < 20:
        return None

    # Extract data
    Z = np.array([c["ecDNA_count"] for c in cells], dtype=float)
    D = np.array([c.get("egfr_expression", c["ecDNA_count"]) for c in cells], dtype=float)
    Y = np.array([c[outcome_name] for c in cells], dtype=float)
    O2 = np.array([c.get("O2_local", 20.0) for c in cells], dtype=float)
    glucose = np.array([c.get("glucose_local", 5.0) for c in cells], dtype=float)

    covariates = np.column_stack([O2, glucose])

    try:
        (first_coef, first_se, f_stat,
         second_coef, second_se, ols_coef) = two_stage_least_squares(
            Z, D, Y, covariates
        )

        # Confidence interval
        z_crit = stats.norm.ppf(1 - alpha / 2)
        ci_lower = second_coef - z_crit * second_se
        ci_upper = second_coef + z_crit * second_se

        # OLS SE
        n = len(cells)
        X_ols = np.column_stack([np.ones(n), D, covariates])
        _, se_ols, _ = np.linalg.lstsq(X_ols.T @ X_ols, X_ols.T @ Y, rcond=None)
        resid_ols = Y - X_ols @ np.linalg.lstsq(X_ols, Y, rcond=None)[0]
        mse_ols = np.sum(resid_ols ** 2) / (n - X_ols.shape[1])
        XtX_inv = np.linalg.pinv(X_ols.T @ X_ols)
        ols_se = np.sqrt(mse_ols * XtX_inv[1, 1])

        return StratifiedEffect(
            stratum_name="",  # Set by caller
            n_observations=n,
            iv_estimate=second_coef,
            iv_se=second_se,
            iv_ci_lower=ci_lower,
            iv_ci_upper=ci_upper,
            ols_estimate=ols_coef,
            ols_se=ols_se,
            first_stage_f=f_stat,
        )
    except Exception:
        return None


def stratified_iv_by_region(
    cells_data: list[dict],
    outcome_name: str = "migration_rate",
    margin_width: float = 50.0,
) -> HeterogeneityAnalysis:
    """Estimate causal effects stratified by tumor region.

    Divides tumor into:
    - Core: cells well inside tumor bulk
    - Margin: cells at invasive front
    - Infiltrating: cells beyond main tumor mass

    This analysis reveals whether EGFR effects on migration/invasion
    vary by spatial context, with therapeutic implications.

    Args:
        cells_data: Cell data from simulation
        outcome_name: Outcome variable
        margin_width: Width of margin region in um

    Returns:
        HeterogeneityAnalysis with regional effects
    """
    # Filter to tumor cells
    tumor_cells = [c for c in cells_data if c.get("cell_type") == 6]

    if len(tumor_cells) < 50:
        return HeterogeneityAnalysis(
            analysis_type="spatial_region",
            outcome_name=outcome_name,
            overall_effect=0.0,
            stratified_effects=[],
            heterogeneity_test={},
            interpretation="Insufficient tumor cells for regional analysis",
        )

    # Compute tumor geometry
    center = _compute_tumor_center(cells_data)
    radius = _compute_tumor_radius(cells_data, center)

    # Classify cells
    regions = {"core": [], "margin": [], "infiltrating": []}
    for cell in tumor_cells:
        region = _classify_tumor_region(cell, center, radius, margin_width)
        regions[region].append(cell)

    # Overall effect
    overall_result = _estimate_iv_for_stratum(tumor_cells, outcome_name)
    overall_effect = overall_result.iv_estimate if overall_result else 0.0

    # Stratified effects
    stratified_effects = []
    for region_name, region_cells in regions.items():
        if len(region_cells) < 20:
            continue
        result = _estimate_iv_for_stratum(region_cells, outcome_name)
        if result:
            result = StratifiedEffect(
                stratum_name=region_name,
                n_observations=result.n_observations,
                iv_estimate=result.iv_estimate,
                iv_se=result.iv_se,
                iv_ci_lower=result.iv_ci_lower,
                iv_ci_upper=result.iv_ci_upper,
                ols_estimate=result.ols_estimate,
                ols_se=result.ols_se,
                first_stage_f=result.first_stage_f,
            )
            stratified_effects.append(result)

    # Test for heterogeneity (Cochran's Q test)
    if len(stratified_effects) >= 2:
        effects = [s.iv_estimate for s in stratified_effects]
        ses = [s.iv_se for s in stratified_effects]
        weights = [1 / (se ** 2) if se > 0 else 0 for se in ses]
        weighted_mean = sum(w * e for w, e in zip(weights, effects)) / sum(weights) if sum(weights) > 0 else 0

        Q = sum(w * (e - weighted_mean) ** 2 for w, e in zip(weights, effects))
        df = len(effects) - 1
        p_value = 1 - stats.chi2.cdf(Q, df) if df > 0 else 1.0

        heterogeneity_test = {
            "Q_statistic": Q,
            "df": df,
            "p_value": p_value,
            "significant": p_value < 0.05,
        }
    else:
        heterogeneity_test = {}

    # Interpretation
    if stratified_effects:
        max_effect = max(stratified_effects, key=lambda x: x.iv_estimate)
        min_effect = min(stratified_effects, key=lambda x: x.iv_estimate)
        interpretation = (
            f"Effect of EGFR on {outcome_name} varies by region. "
            f"Strongest effect at {max_effect.stratum_name} "
            f"(δ={max_effect.iv_estimate:.4f}, 95% CI: [{max_effect.iv_ci_lower:.4f}, {max_effect.iv_ci_upper:.4f}]), "
            f"weakest at {min_effect.stratum_name} "
            f"(δ={min_effect.iv_estimate:.4f}). "
        )
        if heterogeneity_test.get("significant"):
            interpretation += "Heterogeneity is statistically significant (p<0.05)."
        else:
            interpretation += "Heterogeneity is not statistically significant."
    else:
        interpretation = "Insufficient data for regional stratification."

    return HeterogeneityAnalysis(
        analysis_type="spatial_region",
        outcome_name=outcome_name,
        overall_effect=overall_effect,
        stratified_effects=stratified_effects,
        heterogeneity_test=heterogeneity_test,
        interpretation=interpretation,
    )


def stratified_iv_by_hypoxia(
    cells_data: list[dict],
    outcome_name: str = "migration_rate",
    hypoxia_threshold: float = 10.0,
) -> HeterogeneityAnalysis:
    """Estimate causal effects stratified by hypoxia status.

    Implements analysis of effect modification by oxygen status:
    - Do EGFR effects differ in hypoxic vs normoxic cells?
    - Is the "Go-or-Grow" switch captured in the data?

    Args:
        cells_data: Cell data from simulation
        outcome_name: Outcome variable
        hypoxia_threshold: O2 threshold for hypoxia (mmHg)

    Returns:
        HeterogeneityAnalysis with hypoxia-stratified effects
    """
    tumor_cells = [c for c in cells_data if c.get("cell_type") == 6]

    if len(tumor_cells) < 50:
        return HeterogeneityAnalysis(
            analysis_type="hypoxia_status",
            outcome_name=outcome_name,
            overall_effect=0.0,
            stratified_effects=[],
            heterogeneity_test={},
            interpretation="Insufficient tumor cells",
        )

    # Stratify by hypoxia
    hypoxic_cells = [c for c in tumor_cells if c.get("O2_local", 100) < hypoxia_threshold]
    normoxic_cells = [c for c in tumor_cells if c.get("O2_local", 100) >= hypoxia_threshold]

    # Overall effect
    overall_result = _estimate_iv_for_stratum(tumor_cells, outcome_name)
    overall_effect = overall_result.iv_estimate if overall_result else 0.0

    stratified_effects = []

    # Normoxic effect
    if len(normoxic_cells) >= 20:
        result = _estimate_iv_for_stratum(normoxic_cells, outcome_name)
        if result:
            stratified_effects.append(StratifiedEffect(
                stratum_name="normoxic",
                n_observations=result.n_observations,
                iv_estimate=result.iv_estimate,
                iv_se=result.iv_se,
                iv_ci_lower=result.iv_ci_lower,
                iv_ci_upper=result.iv_ci_upper,
                ols_estimate=result.ols_estimate,
                ols_se=result.ols_se,
                first_stage_f=result.first_stage_f,
            ))

    # Hypoxic effect
    if len(hypoxic_cells) >= 20:
        result = _estimate_iv_for_stratum(hypoxic_cells, outcome_name)
        if result:
            stratified_effects.append(StratifiedEffect(
                stratum_name="hypoxic",
                n_observations=result.n_observations,
                iv_estimate=result.iv_estimate,
                iv_se=result.iv_se,
                iv_ci_lower=result.iv_ci_lower,
                iv_ci_upper=result.iv_ci_upper,
                ols_estimate=result.ols_estimate,
                ols_se=result.ols_se,
                first_stage_f=result.first_stage_f,
            ))

    # Heterogeneity test
    heterogeneity_test = {}
    if len(stratified_effects) == 2:
        diff = stratified_effects[0].iv_estimate - stratified_effects[1].iv_estimate
        se_diff = np.sqrt(stratified_effects[0].iv_se ** 2 + stratified_effects[1].iv_se ** 2)
        z_stat = abs(diff) / se_diff if se_diff > 0 else 0
        p_value = 2 * (1 - stats.norm.cdf(z_stat))
        heterogeneity_test = {
            "difference": diff,
            "se_difference": se_diff,
            "z_statistic": z_stat,
            "p_value": p_value,
            "significant": p_value < 0.05,
        }

    interpretation = f"Analysis of {outcome_name} by hypoxia status: "
    if len(stratified_effects) == 2:
        interpretation += (
            f"Normoxic effect: {stratified_effects[0].iv_estimate:.4f}, "
            f"Hypoxic effect: {stratified_effects[1].iv_estimate:.4f}. "
        )
        if heterogeneity_test.get("significant"):
            interpretation += "Effect modification by hypoxia is significant."
        else:
            interpretation += "No significant effect modification."
    else:
        interpretation += "Insufficient cells in one or both strata."

    return HeterogeneityAnalysis(
        analysis_type="hypoxia_status",
        outcome_name=outcome_name,
        overall_effect=overall_effect,
        stratified_effects=stratified_effects,
        heterogeneity_test=heterogeneity_test,
        interpretation=interpretation,
    )


def stratified_iv_by_ecDNA_level(
    cells_data: list[dict],
    outcome_name: str = "migration_rate",
    percentiles: list[float] | None = None,
) -> HeterogeneityAnalysis:
    """Estimate effects stratified by ecDNA copy number.

    Tests whether effect magnitudes differ by ecDNA level
    (non-linear dose-response).

    Args:
        cells_data: Cell data
        outcome_name: Outcome variable
        percentiles: Percentile cutoffs for stratification

    Returns:
        HeterogeneityAnalysis
    """
    if percentiles is None:
        percentiles = [33, 66]

    tumor_cells = [c for c in cells_data if c.get("cell_type") == 6]

    if len(tumor_cells) < 60:
        return HeterogeneityAnalysis(
            analysis_type="ecDNA_level",
            outcome_name=outcome_name,
            overall_effect=0.0,
            stratified_effects=[],
            heterogeneity_test={},
            interpretation="Insufficient cells",
        )

    ecDNA_values = [c["ecDNA_count"] for c in tumor_cells]
    cutoffs = [np.percentile(ecDNA_values, p) for p in percentiles]

    # Stratify
    strata = {
        "low": [c for c in tumor_cells if c["ecDNA_count"] < cutoffs[0]],
        "medium": [c for c in tumor_cells if cutoffs[0] <= c["ecDNA_count"] < cutoffs[1]],
        "high": [c for c in tumor_cells if c["ecDNA_count"] >= cutoffs[1]],
    }

    overall_result = _estimate_iv_for_stratum(tumor_cells, outcome_name)
    overall_effect = overall_result.iv_estimate if overall_result else 0.0

    stratified_effects = []
    for stratum_name, stratum_cells in strata.items():
        if len(stratum_cells) >= 20:
            result = _estimate_iv_for_stratum(stratum_cells, outcome_name)
            if result:
                stratified_effects.append(StratifiedEffect(
                    stratum_name=stratum_name,
                    n_observations=result.n_observations,
                    iv_estimate=result.iv_estimate,
                    iv_se=result.iv_se,
                    iv_ci_lower=result.iv_ci_lower,
                    iv_ci_upper=result.iv_ci_upper,
                    ols_estimate=result.ols_estimate,
                    ols_se=result.ols_se,
                    first_stage_f=result.first_stage_f,
                ))

    interpretation = f"Effect of EGFR on {outcome_name} by ecDNA level: "
    for s in stratified_effects:
        interpretation += f"{s.stratum_name}: {s.iv_estimate:.4f} (n={s.n_observations}), "

    return HeterogeneityAnalysis(
        analysis_type="ecDNA_level",
        outcome_name=outcome_name,
        overall_effect=overall_effect,
        stratified_effects=stratified_effects,
        heterogeneity_test={},
        interpretation=interpretation,
    )


def temporal_effect_evolution(
    cells_by_timestep: dict[int, list[dict]],
    outcome_name: str = "migration_rate",
    min_cells_per_step: int = 30,
) -> dict[str, Any]:
    """Analyze how causal effects evolve over simulation time.

    Tracks effect estimates across timesteps to detect:
    - Effect stability over time
    - Selection-driven changes
    - Microenvironmental evolution effects

    Args:
        cells_by_timestep: Dictionary mapping timestep to cell list
        outcome_name: Outcome variable
        min_cells_per_step: Minimum cells required per timestep

    Returns:
        Dictionary with temporal effect trajectory
    """
    timesteps = sorted(cells_by_timestep.keys())
    results = {
        "timesteps": [],
        "iv_estimates": [],
        "iv_se": [],
        "ols_estimates": [],
        "n_cells": [],
    }

    for t in timesteps:
        cells = cells_by_timestep[t]
        tumor_cells = [c for c in cells if c.get("cell_type") == 6]

        if len(tumor_cells) < min_cells_per_step:
            continue

        result = _estimate_iv_for_stratum(tumor_cells, outcome_name)
        if result:
            results["timesteps"].append(t)
            results["iv_estimates"].append(result.iv_estimate)
            results["iv_se"].append(result.iv_se)
            results["ols_estimates"].append(result.ols_estimate)
            results["n_cells"].append(result.n_observations)

    if len(results["timesteps"]) < 2:
        results["trend_analysis"] = {"error": "Insufficient timesteps"}
        return results

    # Trend analysis
    iv_est = np.array(results["iv_estimates"])
    times = np.array(results["timesteps"])

    # Linear regression for trend
    slope, intercept, r_value, p_value, std_err = stats.linregress(times, iv_est)

    results["trend_analysis"] = {
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_value ** 2,
        "p_value": p_value,
        "significant_trend": p_value < 0.05,
        "interpretation": (
            f"Effect {'increases' if slope > 0 else 'decreases'} over time "
            f"(slope={slope:.6f}/hour, p={p_value:.3f})"
        ),
    }

    return results


def radial_effect_profile(
    cells_data: list[dict],
    outcome_name: str = "migration_rate",
    n_bins: int = 5,
) -> dict[str, Any]:
    """Compute effect as function of distance from tumor center.

    Creates a radial profile showing how the causal effect varies
    with distance, useful for identifying spatial gradients.

    Args:
        cells_data: Cell data
        outcome_name: Outcome variable
        n_bins: Number of radial bins

    Returns:
        Dictionary with radial effect profile
    """
    tumor_cells = [c for c in cells_data if c.get("cell_type") == 6]

    if len(tumor_cells) < 50:
        return {"error": "Insufficient cells"}

    center = _compute_tumor_center(cells_data)

    # Compute distances
    for cell in tumor_cells:
        cell["_distance"] = np.sqrt(
            (cell["x"] - center[0]) ** 2 + (cell["y"] - center[1]) ** 2
        )

    distances = [c["_distance"] for c in tumor_cells]
    bin_edges = np.linspace(0, np.max(distances), n_bins + 1)

    results = {
        "bin_centers": [],
        "bin_ranges": [],
        "iv_estimates": [],
        "iv_se": [],
        "n_cells": [],
    }

    for i in range(n_bins):
        bin_cells = [
            c for c in tumor_cells
            if bin_edges[i] <= c["_distance"] < bin_edges[i + 1]
        ]

        if len(bin_cells) < 20:
            continue

        result = _estimate_iv_for_stratum(bin_cells, outcome_name)
        if result:
            results["bin_centers"].append((bin_edges[i] + bin_edges[i + 1]) / 2)
            results["bin_ranges"].append((bin_edges[i], bin_edges[i + 1]))
            results["iv_estimates"].append(result.iv_estimate)
            results["iv_se"].append(result.iv_se)
            results["n_cells"].append(result.n_observations)

    # Clean up temporary field
    for cell in tumor_cells:
        del cell["_distance"]

    return results


def comprehensive_heterogeneity_analysis(
    cells_data: list[dict],
    outcomes: list[str] | None = None,
) -> dict[str, Any]:
    """Run all heterogeneity analyses for comprehensive report.

    Args:
        cells_data: Cell data from final timestep
        outcomes: List of outcome variables

    Returns:
        Comprehensive results dictionary
    """
    if outcomes is None:
        outcomes = ["migration_rate", "VEGF_secretion"]

    results = {
        "by_region": {},
        "by_hypoxia": {},
        "by_ecDNA_level": {},
        "radial_profiles": {},
    }

    for outcome in outcomes:
        if outcome not in cells_data[0]:
            continue

        # Regional analysis
        regional = stratified_iv_by_region(cells_data, outcome)
        results["by_region"][outcome] = regional.to_dict()

        # Hypoxia analysis
        hypoxia = stratified_iv_by_hypoxia(cells_data, outcome)
        results["by_hypoxia"][outcome] = hypoxia.to_dict()

        # ecDNA level analysis
        ecdna = stratified_iv_by_ecDNA_level(cells_data, outcome)
        results["by_ecDNA_level"][outcome] = ecdna.to_dict()

        # Radial profile
        radial = radial_effect_profile(cells_data, outcome)
        results["radial_profiles"][outcome] = radial

    return results
