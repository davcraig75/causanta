"""Causal effect estimation from CAUSANTA simulation data.

Provides methods to estimate the ground-truth causal effects (α, β, γ, δ)
from simulation output, enabling comparison of estimated vs configured values.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

# Use loader.py for all data loading - single source of truth
from .loader import load_cells_tsv, load_lineage_tsv, load_summary_log


@dataclass
class EffectEstimate:
    """Result of estimating a causal effect from data."""

    parameter_symbol: str
    parameter_name: str
    configured_value: float
    estimated_value: float
    standard_error: float
    n_observations: int
    method: str
    notes: str = ""

    @property
    def relative_error(self) -> float:
        """Relative error between estimated and configured value."""
        if self.configured_value == 0:
            return float('inf') if self.estimated_value != 0 else 0.0
        return abs(self.estimated_value - self.configured_value) / abs(self.configured_value)

    @property
    def within_2se(self) -> bool:
        """Whether configured value is within 2 standard errors of estimate."""
        if self.standard_error == 0:
            return self.estimated_value == self.configured_value
        return abs(self.estimated_value - self.configured_value) <= 2 * self.standard_error

    def to_dict(self) -> dict[str, Any]:
        return {
            "parameter_symbol": self.parameter_symbol,
            "parameter_name": self.parameter_name,
            "configured_value": self.configured_value,
            "estimated_value": self.estimated_value,
            "standard_error": self.standard_error,
            "n_observations": self.n_observations,
            "relative_error": self.relative_error,
            "within_2se": self.within_2se,
            "method": self.method,
            "notes": self.notes,
        }


def estimate_alpha(
    lineage_data: list[dict],
    base_division_time: float,
) -> EffectEstimate:
    """[DEPRECATED PROXY] Descriptive statistic of ecDNA at division.

    .. deprecated::
        This returns ``mean(log2(1+ecDNA)) / base_division_time`` — a descriptive
        statistic, NOT an estimator of α on its structural scale, and it does not
        recover the configured value. The proper structural recovery of α is
        ``alpha_recovery()`` in ``scripts/recover_structural_effects.py``, which
        fits Δt = T_base/(1 + α·log₂(1+EGFR)) on inter-division intervals with
        T_base known, identified on κ_hyp = 0 runs (recovers α ≈ 0.31, truth 0.30).
        Retained only for exploratory back-compatibility; do not report its value.

    Uses the relationship: T_eff = T_base / (1 + α * log2(1 + ecDNA))
    """
    if not lineage_data:
        return EffectEstimate(
            parameter_symbol="α",
            parameter_name="ecDNA_effect_on_division",
            configured_value=0.0,
            estimated_value=0.0,
            standard_error=0.0,
            n_observations=0,
            method="lineage_regression",
            notes="No lineage data available",
        )

    # Group divisions by parent ecDNA count
    # Higher ecDNA should lead to more divisions (faster cycling)
    ecDNA_counts = []
    for rec in lineage_data:
        ecDNA_counts.append(rec["parent_ecDNA_before"])

    if len(ecDNA_counts) < 2:
        return EffectEstimate(
            parameter_symbol="α",
            parameter_name="ecDNA_effect_on_division",
            configured_value=0.0,
            estimated_value=0.0,
            standard_error=0.0,
            n_observations=len(ecDNA_counts),
            method="lineage_regression",
            notes="Insufficient division events for estimation",
        )

    # Compute division rate proxy by ecDNA bin
    # We use log2(1 + ecDNA) as the predictor
    X = np.array([math.log2(1 + e) for e in ecDNA_counts])

    # Since we don't have actual cycle times, we estimate based on
    # the distribution of ecDNA at division time.
    # Higher α means cells with more ecDNA divide faster, so they
    # should be over-represented in the division events.

    # Simple approach: fit the model assuming observed ecDNA distribution
    # reflects the acceleration. Compare mean log2(1+ecDNA) to expected.
    mean_log_ecDNA = np.mean(X)
    var_log_ecDNA = np.var(X, ddof=1) if len(X) > 1 else 0.0
    se_log_ecDNA = np.sqrt(var_log_ecDNA / len(X)) if len(X) > 0 else 0.0

    # Rough estimate: if cells divide at rate proportional to (1 + α*log2(1+ecDNA)),
    # then dividing cells will have elevated mean log2(1+ecDNA) proportional to α.
    # This is a simplified proxy; actual estimation would require more data.

    # For now, we report the mean and note this is a proxy estimate
    estimated_alpha = mean_log_ecDNA / base_division_time if base_division_time > 0 else 0.0

    return EffectEstimate(
        parameter_symbol="α",
        parameter_name="ecDNA_effect_on_division",
        configured_value=0.0,  # Will be filled by caller
        estimated_value=estimated_alpha,
        standard_error=se_log_ecDNA,
        n_observations=len(ecDNA_counts),
        method="division_ecDNA_distribution",
        notes=f"Mean log2(1+ecDNA) at division: {mean_log_ecDNA:.2f}",
    )


def estimate_beta(
    cells_data: list[dict],
    base_vegf_secretion: float,
) -> EffectEstimate:
    """Estimate β (ecDNA effect on VEGF) from cell snapshot data.

    Simulator uses: S_eff = S_base * (1 + β * sqrt(EGFR))
    Rearranging: VEGF / S_base - 1 = β * sqrt(EGFR)
    So regress (VEGF/S_base - 1) on sqrt(EGFR) to get β.
    """
    # Filter to tumor cells with VEGF data
    tumor_cells = [c for c in cells_data if c["cell_type"] == 6
                   and c.get("egfr_expression", 0) > 0]

    if len(tumor_cells) < 2:
        return EffectEstimate(
            parameter_symbol="β",
            parameter_name="ecDNA_effect_on_VEGF",
            configured_value=0.0,
            estimated_value=0.0,
            standard_error=0.0,
            n_observations=len(tumor_cells),
            method="sqrt_egfr_regression",
            notes="Insufficient tumor cells for estimation",
        )

    # Use sqrt(EGFR) as regressor to match simulator's functional form
    sqrt_egfr = np.array([math.sqrt(c["egfr_expression"]) for c in tumor_cells])
    vegf = np.array([c["VEGF_secretion"] for c in tumor_cells], dtype=float)

    # Regression: VEGF = S_base + S_base * β * sqrt(EGFR)
    n = len(sqrt_egfr)
    mean_x = np.mean(sqrt_egfr)
    mean_y = np.mean(vegf)

    numerator = np.sum((sqrt_egfr - mean_x) * (vegf - mean_y))
    denominator = np.sum((sqrt_egfr - mean_x) ** 2)

    if denominator == 0:
        return EffectEstimate(
            parameter_symbol="β",
            parameter_name="ecDNA_effect_on_VEGF",
            configured_value=0.0,
            estimated_value=0.0,
            standard_error=0.0,
            n_observations=n,
            method="sqrt_egfr_regression",
            notes="No variance in EGFR expression",
        )

    slope = numerator / denominator
    # slope = S_base * β, so β = slope / S_base
    estimated_beta = slope / base_vegf_secretion if base_vegf_secretion > 0 else 0.0

    # Standard error
    intercept = mean_y - slope * mean_x
    residuals = vegf - (intercept + slope * sqrt_egfr)
    mse = np.sum(residuals ** 2) / (n - 2) if n > 2 else 0.0
    se_slope = np.sqrt(mse / denominator) if denominator > 0 else 0.0
    se_beta = se_slope / base_vegf_secretion if base_vegf_secretion > 0 else 0.0

    return EffectEstimate(
        parameter_symbol="β",
        parameter_name="ecDNA_effect_on_VEGF",
        configured_value=0.0,
        estimated_value=estimated_beta,
        standard_error=se_beta,
        n_observations=n,
        method="sqrt_egfr_regression",
        notes=f"Regression: VEGF = {intercept:.1f} + {slope:.3f} * sqrt(EGFR)",
    )


def estimate_delta(
    cells_data: list[dict],
    base_migration_speed: float,
) -> EffectEstimate:
    """Estimate δ (ecDNA effect on migration) from cell snapshot data.

    Simulator uses: v_eff = v_base * (1 + δ * EGFR) * hypoxia_mult
    To remove hypoxia confounding, we filter to normoxic cells only.
    Then: speed = v_base + v_base * δ * EGFR, so δ = slope / v_base.
    """
    # Filter to normoxic tumor cells to remove Go-or-Grow confounding
    tumor_cells = [c for c in cells_data if c["cell_type"] == 6
                   and c.get("egfr_expression", 0) > 0
                   and not c.get("is_hypoxic", False)]

    if len(tumor_cells) < 2:
        return EffectEstimate(
            parameter_symbol="δ",
            parameter_name="ecDNA_effect_on_migration",
            configured_value=0.0,
            estimated_value=0.0,
            standard_error=0.0,
            n_observations=len(tumor_cells),
            method="egfr_regression_normoxic",
            notes="Insufficient normoxic tumor cells for estimation",
        )

    egfr = np.array([c["egfr_expression"] for c in tumor_cells], dtype=float)
    migration = np.array([c["migration_rate"] for c in tumor_cells], dtype=float)

    n = len(egfr)
    mean_x = np.mean(egfr)
    mean_y = np.mean(migration)

    numerator = np.sum((egfr - mean_x) * (migration - mean_y))
    denominator = np.sum((egfr - mean_x) ** 2)

    if denominator == 0 or base_migration_speed == 0:
        return EffectEstimate(
            parameter_symbol="δ",
            parameter_name="ecDNA_effect_on_migration",
            configured_value=0.0,
            estimated_value=0.0,
            standard_error=0.0,
            n_observations=n,
            method="egfr_regression_normoxic",
            notes="No variance in EGFR or zero base speed",
        )

    slope = numerator / denominator
    intercept = mean_y - slope * mean_x
    # slope = v_base * δ, so δ = slope / v_base
    estimated_delta = slope / base_migration_speed

    residuals = migration - (intercept + slope * egfr)
    mse = np.sum(residuals ** 2) / (n - 2) if n > 2 else 0.0
    se_slope = np.sqrt(mse / denominator) if denominator > 0 else 0.0
    se_delta = se_slope / base_migration_speed

    return EffectEstimate(
        parameter_symbol="δ",
        parameter_name="ecDNA_effect_on_migration",
        configured_value=0.0,
        estimated_value=estimated_delta,
        standard_error=se_delta,
        n_observations=n,
        method="egfr_regression_normoxic",
        notes=f"Regression on normoxic cells: speed = {intercept:.1f} + {slope:.3f} * EGFR",
    )


def estimate_gamma_from_lineage(
    lineage_data: list[dict],
    final_cells_data: list[dict],
) -> EffectEstimate:
    """[DEPRECATED PROXY] γ (ecDNA effect on survival) is NOT identified here.

    .. deprecated::
        γ is a non-identified design parameter in these runs. Baseline apoptosis
        is 5e-5/hr, so essentially no tumor cell dies of apoptosis and γ (which
        only *further* suppresses it) leaves no survival-selection footprint. This
        function returns a raw ecDNA "enrichment ratio" of surviving vs all-ever
        cells, which is dominated by replication/accumulation dynamics, is not on
        γ's scale, and empirically returns ≈0 (or negative). Do NOT report it as a
        recovered γ. To make γ identifiable one would need a non-negligible base
        apoptosis rate plus a logged death-event hazard regression on EGFR.

    Uses the relationship: a_eff = a_base / (1 + γ * log2(1 + EGFR))
    """
    if not lineage_data or not final_cells_data:
        return EffectEstimate(
            parameter_symbol="γ",
            parameter_name="ecDNA_effect_on_survival",
            configured_value=0.0,
            estimated_value=0.0,
            standard_error=0.0,
            n_observations=0,
            method="survival_comparison",
            notes="Insufficient data for survival analysis",
        )

    # Get ecDNA of all tumor cells that ever existed (from lineage)
    all_ecDNA = []
    for rec in lineage_data:
        all_ecDNA.append(rec["parent_ecDNA_before"])
        all_ecDNA.append(rec["daughter_ecDNA"])

    # Get ecDNA of surviving tumor cells
    surviving_tumor = [c for c in final_cells_data if c["cell_type"] == 6]
    surviving_ecDNA = [c["ecDNA_count"] for c in surviving_tumor]

    if len(all_ecDNA) < 2 or len(surviving_ecDNA) < 2:
        return EffectEstimate(
            parameter_symbol="γ",
            parameter_name="ecDNA_effect_on_survival",
            configured_value=0.0,
            estimated_value=0.0,
            standard_error=0.0,
            n_observations=len(surviving_ecDNA),
            method="survival_comparison",
            notes="Insufficient cells for comparison",
        )

    mean_all = np.mean(all_ecDNA)
    mean_surviving = np.mean(surviving_ecDNA)

    # If γ > 0, surviving cells should have higher mean ecDNA
    # Simple estimate: enrichment ratio as proxy for γ
    if mean_all > 0:
        enrichment = (mean_surviving - mean_all) / mean_all
        estimated_gamma = max(0, enrichment)  # γ should be non-negative
    else:
        estimated_gamma = 0.0

    # Standard error from surviving distribution
    se = np.std(surviving_ecDNA, ddof=1) / np.sqrt(len(surviving_ecDNA)) if len(surviving_ecDNA) > 1 else 0.0

    return EffectEstimate(
        parameter_symbol="γ",
        parameter_name="ecDNA_effect_on_survival",
        configured_value=0.0,
        estimated_value=estimated_gamma,
        standard_error=se / mean_all if mean_all > 0 else 0.0,
        n_observations=len(surviving_ecDNA),
        method="survival_enrichment",
        notes=f"Mean ecDNA - all: {mean_all:.1f}, surviving: {mean_surviving:.1f}",
    )


@dataclass
class CausalEffectAnalysis:
    """Complete causal effect analysis results."""

    alpha: EffectEstimate
    beta: EffectEstimate
    delta: EffectEstimate
    gamma: EffectEstimate
    segregation_stats: dict[str, float]
    summary: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "alpha": self.alpha.to_dict(),
            "beta": self.beta.to_dict(),
            "delta": self.delta.to_dict(),
            "gamma": self.gamma.to_dict(),
            "segregation_stats": self.segregation_stats,
            "summary": self.summary,
        }


def analyze_ecDNA_segregation(lineage_data: list[dict]) -> dict[str, float]:
    """Analyze ecDNA segregation statistics from lineage data.

    Tests whether segregation follows expected Binomial(N, 0.5) distribution.
    """
    if not lineage_data:
        return {
            "n_divisions": 0,
            "mean_parent_before": 0.0,
            "mean_daughter_fraction": 0.0,
            "expected_fraction": 0.5,
            "segregation_variance": 0.0,
            "expected_variance": 0.0,
        }

    parent_before = []
    daughter_fractions = []

    for rec in lineage_data:
        n_before = rec["parent_ecDNA_before"]
        n_after = rec["parent_ecDNA_after"]
        n_daughter = rec["daughter_ecDNA"]
        # The total post-replication pool is parent_after + daughter
        # (replication happens during S-phase, then segregation at M-phase)
        total_post_replication = n_after + n_daughter
        if total_post_replication > 0:
            parent_before.append(n_before)
            # Daughter fraction of the post-replication pool should be ~0.5
            daughter_fractions.append(n_daughter / total_post_replication)

    if not daughter_fractions:
        return {
            "n_divisions": len(lineage_data),
            "mean_parent_before": 0.0,
            "mean_daughter_fraction": 0.0,
            "expected_fraction": 0.5,
            "segregation_variance": 0.0,
            "expected_variance": 0.0,
            "notes": "Fractions computed relative to post-replication pool",
        }

    mean_fraction = np.mean(daughter_fractions)
    var_fraction = np.var(daughter_fractions, ddof=1) if len(daughter_fractions) > 1 else 0.0
    mean_n = np.mean(parent_before)

    # Post-replication pool is ~2*N (with 95% replication fidelity)
    # Expected variance for Binomial(2N, 0.5) / (2N) = 0.25/(2N)
    mean_replicated = np.mean([rec["parent_ecDNA_after"] + rec["daughter_ecDNA"]
                               for rec in lineage_data
                               if rec["parent_ecDNA_after"] + rec["daughter_ecDNA"] > 0])
    expected_var = 0.25 / mean_replicated if mean_replicated > 0 else 0.0

    return {
        "n_divisions": len(lineage_data),
        "mean_parent_before": float(mean_n),
        "mean_daughter_fraction": float(mean_fraction),
        "expected_fraction": 0.5,
        "segregation_variance": float(var_fraction),
        "expected_variance": float(expected_var),
        "mean_replicated_pool": float(mean_replicated),
        "notes": "Fractions computed relative to post-replication pool (parent_after + daughter)",
    }


def estimate_effects_iv(
    cells_data: list[dict],
    ground_truth: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Estimate all causal effects using IV regression.

    Uses ecDNA as instrument for EGFR expression effects on outcomes.

    Args:
        cells_data: Cell data from simulation
        ground_truth: Optional dictionary with true effect values

    Returns:
        Dictionary with IV estimates for each outcome
    """
    from .iv import two_stage_least_squares

    tumor_cells = [c for c in cells_data if c.get("cell_type") == 6]

    if len(tumor_cells) < 50:
        return {"error": "Insufficient tumor cells", "n_cells": len(tumor_cells)}

    # Instrument and treatment
    Z = np.array([c["ecDNA_count"] for c in tumor_cells], dtype=float)
    D = np.array([c.get("egfr_expression", c["ecDNA_count"]) for c in tumor_cells], dtype=float)

    # Covariates (to control for observed confounding)
    O2 = np.array([c.get("O2_local", 20) for c in tumor_cells], dtype=float)
    covariates = O2.reshape(-1, 1)

    results = {
        "n_cells": len(tumor_cells),
        "instrument": "ecDNA_count",
        "treatment": "egfr_expression",
    }

    # Outcome: VEGF secretion (beta)
    if "VEGF_secretion" in tumor_cells[0]:
        Y_vegf = np.array([c["VEGF_secretion"] for c in tumor_cells], dtype=float)
        try:
            first_coef, first_se, f_stat, iv_coef, iv_se, ols_coef = two_stage_least_squares(
                Z, D, Y_vegf, covariates
            )
            results["VEGF_effect"] = {
                "iv_estimate": iv_coef,
                "iv_se": iv_se,
                "ols_estimate": ols_coef,
                "first_stage_coef": first_coef,
                "f_statistic": f_stat,
                "ground_truth": ground_truth.get("beta", None) if ground_truth else None,
            }
        except Exception as e:
            results["VEGF_effect"] = {"error": str(e)}

    # Outcome: migration rate (delta)
    if "migration_rate" in tumor_cells[0]:
        Y_mig = np.array([c["migration_rate"] for c in tumor_cells], dtype=float)
        try:
            first_coef, first_se, f_stat, iv_coef, iv_se, ols_coef = two_stage_least_squares(
                Z, D, Y_mig, covariates
            )
            results["migration_effect"] = {
                "iv_estimate": iv_coef,
                "iv_se": iv_se,
                "ols_estimate": ols_coef,
                "first_stage_coef": first_coef,
                "f_statistic": f_stat,
                "ground_truth": ground_truth.get("delta", None) if ground_truth else None,
            }
        except Exception as e:
            results["migration_effect"] = {"error": str(e)}

    return results


def compare_iv_ols_performance(
    cells_data: list[dict],
    ground_truth: dict[str, float],
    n_bootstrap: int = 100,
) -> dict[str, Any]:
    """Compare IV vs OLS estimation performance.

    Runs bootstrap analysis to compare bias and variance of
    IV vs OLS estimators.

    Args:
        cells_data: Cell data from simulation
        ground_truth: True causal effects
        n_bootstrap: Number of bootstrap samples

    Returns:
        Comparison results
    """
    from .iv import two_stage_least_squares

    tumor_cells = [c for c in cells_data if c.get("cell_type") == 6]
    n = len(tumor_cells)

    if n < 100:
        return {"error": "Insufficient cells for bootstrap"}

    Z = np.array([c["ecDNA_count"] for c in tumor_cells], dtype=float)
    D = np.array([c.get("egfr_expression", c["ecDNA_count"]) for c in tumor_cells], dtype=float)
    Y = np.array([c.get("VEGF_secretion", 0) for c in tumor_cells], dtype=float)
    O2 = np.array([c.get("O2_local", 20) for c in tumor_cells], dtype=float)
    covariates = O2.reshape(-1, 1)

    iv_estimates = []
    ols_estimates = []

    rng = np.random.default_rng(42)

    for _ in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        Z_b = Z[idx]
        D_b = D[idx]
        Y_b = Y[idx]
        cov_b = covariates[idx]

        try:
            _, _, _, iv_coef, _, ols_coef = two_stage_least_squares(Z_b, D_b, Y_b, cov_b)
            iv_estimates.append(iv_coef)
            ols_estimates.append(ols_coef)
        except Exception:
            continue

    if len(iv_estimates) < 10:
        return {"error": "Bootstrap failed"}

    true_effect = ground_truth.get("beta", 0.1)

    iv_arr = np.array(iv_estimates)
    ols_arr = np.array(ols_estimates)

    return {
        "true_effect": true_effect,
        "iv": {
            "mean": float(np.mean(iv_arr)),
            "std": float(np.std(iv_arr)),
            "bias": float(np.mean(iv_arr) - true_effect),
            "rmse": float(np.sqrt(np.mean((iv_arr - true_effect) ** 2))),
            "ci_lower": float(np.percentile(iv_arr, 2.5)),
            "ci_upper": float(np.percentile(iv_arr, 97.5)),
            "coverage": float(np.mean((iv_arr - 1.96 * np.std(iv_arr) <= true_effect) &
                                      (true_effect <= iv_arr + 1.96 * np.std(iv_arr)))),
        },
        "ols": {
            "mean": float(np.mean(ols_arr)),
            "std": float(np.std(ols_arr)),
            "bias": float(np.mean(ols_arr) - true_effect),
            "rmse": float(np.sqrt(np.mean((ols_arr - true_effect) ** 2))),
            "ci_lower": float(np.percentile(ols_arr, 2.5)),
            "ci_upper": float(np.percentile(ols_arr, 97.5)),
        },
        "iv_advantage": {
            "bias_reduction": float(abs(np.mean(ols_arr) - true_effect) - abs(np.mean(iv_arr) - true_effect)),
            "rmse_reduction": float(np.sqrt(np.mean((ols_arr - true_effect) ** 2)) -
                                    np.sqrt(np.mean((iv_arr - true_effect) ** 2))),
        },
        "n_bootstrap": len(iv_estimates),
    }


def run_causal_analysis(
    output_dir: Path,
    config_params: dict[str, float],
) -> CausalEffectAnalysis:
    """Run complete causal effect analysis on simulation output.

    Args:
        output_dir: Path to simulation output directory
        config_params: Dictionary with configured effect values:
            - alpha: ecDNA_effect_on_division
            - beta: ecDNA_effect_on_VEGF
            - delta: ecDNA_effect_on_migration
            - gamma: ecDNA_effect_on_survival
            - base_division_time: division_time_mean_hr for tumor
            - base_vegf: VEGF_secretion_amol_hr for tumor
            - base_migration: migration_speed_um_hr for tumor

    Returns:
        CausalEffectAnalysis with all estimates and comparisons.
    """
    # Load data using loader.py functions
    lineage_path = output_dir / "lineage.tsv"
    lineage_data = load_lineage_tsv(lineage_path)

    # Find the final cells file
    cells_files = sorted(output_dir.glob("cells_t*.tsv"))
    final_cells_data = []
    if cells_files:
        final_cells_data = load_cells_tsv(cells_files[-1])

    # Estimate effects
    alpha_est = estimate_alpha(lineage_data, config_params.get("base_division_time", 36.0))
    alpha_est = EffectEstimate(
        parameter_symbol=alpha_est.parameter_symbol,
        parameter_name=alpha_est.parameter_name,
        configured_value=config_params.get("alpha", 0.0),
        estimated_value=alpha_est.estimated_value,
        standard_error=alpha_est.standard_error,
        n_observations=alpha_est.n_observations,
        method=alpha_est.method,
        notes=alpha_est.notes,
    )

    beta_est = estimate_beta(final_cells_data, config_params.get("base_vegf", 600.0))
    beta_est = EffectEstimate(
        parameter_symbol=beta_est.parameter_symbol,
        parameter_name=beta_est.parameter_name,
        configured_value=config_params.get("beta", 0.0),
        estimated_value=beta_est.estimated_value,
        standard_error=beta_est.standard_error,
        n_observations=beta_est.n_observations,
        method=beta_est.method,
        notes=beta_est.notes,
    )

    delta_est = estimate_delta(final_cells_data, config_params.get("base_migration", 10.0))
    delta_est = EffectEstimate(
        parameter_symbol=delta_est.parameter_symbol,
        parameter_name=delta_est.parameter_name,
        configured_value=config_params.get("delta", 0.0),
        estimated_value=delta_est.estimated_value,
        standard_error=delta_est.standard_error,
        n_observations=delta_est.n_observations,
        method=delta_est.method,
        notes=delta_est.notes,
    )

    gamma_est = estimate_gamma_from_lineage(lineage_data, final_cells_data)
    gamma_est = EffectEstimate(
        parameter_symbol=gamma_est.parameter_symbol,
        parameter_name=gamma_est.parameter_name,
        configured_value=config_params.get("gamma", 0.0),
        estimated_value=gamma_est.estimated_value,
        standard_error=gamma_est.standard_error,
        n_observations=gamma_est.n_observations,
        method=gamma_est.method,
        notes=gamma_est.notes,
    )

    # Segregation analysis
    seg_stats = analyze_ecDNA_segregation(lineage_data)

    # Generate summary
    summary_lines = [
        "Causal Effect Estimation Summary",
        "=" * 40,
        "",
        f"Division events analyzed: {len(lineage_data)}",
        f"Final tumor cells: {sum(1 for c in final_cells_data if c['cell_type'] == 6)}",
        "",
        "Effect Estimates (configured vs estimated):",
        f"  α (division):  {config_params.get('alpha', 0):.3f} vs {alpha_est.estimated_value:.3f}",
        f"  β (VEGF):      {config_params.get('beta', 0):.3f} vs {beta_est.estimated_value:.3f}",
        f"  δ (migration): {config_params.get('delta', 0):.3f} vs {delta_est.estimated_value:.3f}",
        f"  γ (survival):  {config_params.get('gamma', 0):.3f} vs {gamma_est.estimated_value:.3f}",
        "",
        "ecDNA Segregation:",
        f"  Mean daughter fraction: {seg_stats['mean_daughter_fraction']:.3f} (expected: 0.5)",
        f"  Variance: {seg_stats['segregation_variance']:.4f} (expected: {seg_stats['expected_variance']:.4f})",
    ]

    return CausalEffectAnalysis(
        alpha=alpha_est,
        beta=beta_est,
        delta=delta_est,
        gamma=gamma_est,
        segregation_stats=seg_stats,
        summary="\n".join(summary_lines),
    )
