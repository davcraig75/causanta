"""Sensitivity analysis for causal inference.

Implements methods to assess how robust causal estimates are to
potential unmeasured confounding and model specification.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import stats

from .loader import SimulationData, cells_to_array


@dataclass
class SensitivityResult:
    """Result of sensitivity analysis."""

    analysis_type: str
    original_estimate: float
    original_se: float
    robustness_value: float  # RV or similar metric
    breakdown_point: float   # At what level of confounding does effect disappear?
    sensitivity_plot_data: dict[str, list[float]]
    interpretation: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "analysis_type": self.analysis_type,
            "original_estimate": self.original_estimate,
            "original_se": self.original_se,
            "robustness_value": self.robustness_value,
            "breakdown_point": self.breakdown_point,
            "plot_data": self.sensitivity_plot_data,
            "interpretation": self.interpretation,
        }


def rosenbaum_bounds(
    data: SimulationData,
    treatment_threshold: float = 10.0,
    outcome_name: str = "VEGF_secretion",
    gamma_range: list[float] | None = None,
) -> SensitivityResult:
    """Rosenbaum sensitivity analysis for matching estimators.

    Assesses how sensitive the treatment effect is to hidden bias
    (unmeasured confounding).

    Args:
        data: Loaded simulation data
        treatment_threshold: ecDNA threshold for treatment
        outcome_name: Outcome variable
        gamma_range: Range of gamma values to test (default: 1 to 3)

    Returns:
        SensitivityResult with bounds at each gamma level.
    """
    if gamma_range is None:
        # Extended range to find actual breakdown point
        gamma_range = [1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0]

    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 20:
        return SensitivityResult(
            analysis_type="rosenbaum_bounds",
            original_estimate=0.0,
            original_se=0.0,
            robustness_value=0.0,
            breakdown_point=0.0,
            sensitivity_plot_data={},
            interpretation="Insufficient data",
        )

    # Get treatment and outcome
    ecDNA = np.array([c["ecDNA_count"] for c in tumor_cells], dtype=float)
    treatment = (ecDNA > treatment_threshold).astype(float)
    outcome = np.array([c[outcome_name] for c in tumor_cells], dtype=float)

    # Simple mean difference
    treated_idx = treatment == 1
    control_idx = treatment == 0

    if np.sum(treated_idx) < 5 or np.sum(control_idx) < 5:
        return SensitivityResult(
            analysis_type="rosenbaum_bounds",
            original_estimate=0.0,
            original_se=0.0,
            robustness_value=0.0,
            breakdown_point=0.0,
            sensitivity_plot_data={},
            interpretation="Insufficient treated or control units",
        )

    effect = np.mean(outcome[treated_idx]) - np.mean(outcome[control_idx])
    se = np.sqrt(
        np.var(outcome[treated_idx], ddof=1) / np.sum(treated_idx) +
        np.var(outcome[control_idx], ddof=1) / np.sum(control_idx)
    )

    # Wilcoxon signed-rank test statistic
    # For simplicity, use paired differences if available
    n_pairs = min(np.sum(treated_idx), np.sum(control_idx))
    treated_sample = outcome[treated_idx][:n_pairs]
    control_sample = outcome[control_idx][:n_pairs]
    diffs = treated_sample - control_sample

    # Compute bounds at each gamma
    upper_bounds = []
    lower_bounds = []
    p_values_upper = []
    p_values_lower = []

    for gamma in gamma_range:
        # Rosenbaum bounds on the Wilcoxon test
        # Upper bound: assume hidden bias makes treatment more likely for high outcomes
        # Lower bound: opposite direction

        # Simplified approximation using normal distribution
        # True Rosenbaum bounds require more complex computation
        ranks = stats.rankdata(np.abs(diffs))
        W_plus = np.sum(ranks[diffs > 0])
        n = len(diffs)

        # Expected value and variance under null
        E_W = n * (n + 1) / 4
        Var_W = n * (n + 1) * (2 * n + 1) / 24

        # Adjust for gamma
        # Upper bound: inflate variance
        Var_upper = Var_W * gamma
        Var_lower = Var_W / gamma

        z_upper = (W_plus - E_W) / np.sqrt(Var_upper)
        z_lower = (W_plus - E_W) / np.sqrt(Var_lower)

        upper_bounds.append(float(z_upper))
        lower_bounds.append(float(z_lower))

        p_values_upper.append(float(1 - stats.norm.cdf(z_upper)))
        p_values_lower.append(float(1 - stats.norm.cdf(z_lower)))

    # Find breakdown point (gamma where p > 0.05)
    breakdown = gamma_range[-1]  # Default to max
    for i, p in enumerate(p_values_upper):
        if p > 0.05:
            breakdown = gamma_range[i]
            break

    # Robustness value
    # How much confounding would need to exist to explain away the effect
    rv = breakdown if breakdown > 1.0 else 0.0

    interpretation = (
        f"Original effect: {effect:.3f} (SE: {se:.3f}). "
        f"Effect remains significant (p<0.05) until gamma={breakdown:.2f}. "
        f"This means an unmeasured confounder would need to change the odds of "
        f"treatment by a factor of {breakdown:.1f} to explain away the effect."
    )

    return SensitivityResult(
        analysis_type="rosenbaum_bounds",
        original_estimate=float(effect),
        original_se=float(se),
        robustness_value=float(rv),
        breakdown_point=float(breakdown),
        sensitivity_plot_data={
            "gamma": gamma_range,
            "upper_bound": upper_bounds,
            "lower_bound": lower_bounds,
            "p_value_upper": p_values_upper,
            "p_value_lower": p_values_lower,
        },
        interpretation=interpretation,
    )


def e_value_analysis(
    estimate: float,
    se: float,
    null_value: float = 0.0,
    log_scale: bool = True,
) -> SensitivityResult:
    """Compute E-value for unmeasured confounding (VanderWeele & Ding 2017).

    Delegates to the canonical implementation in causanta.analyze.iv.compute_e_value,
    which implements the formula E = RR + sqrt(RR(RR - 1)) per VanderWeele &
    Ding, Ann Intern Med 167(4):268-274 (2017), with the standard recipes for
    RR < 1 (flip to 1/RR) and for the CI E-value (bound closest to the null).

    Args:
        estimate: Point estimate. By default treated as log-RR (and
            exponentiated); set log_scale=False to pass an estimate already
            on the RR scale.
        se: Standard error of the estimate (same scale as `estimate`).
        null_value: Null hypothesis value on the supplied scale (default 0
            for log-RR, i.e. RR = 1). Currently used only for documentation;
            the formula uses the canonical null RR = 1 internally.
        log_scale: If True (default), treat `estimate` ± 1.96·SE as log-RR
            and exponentiate before applying the formula.

    Returns:
        SensitivityResult with E-value for the point estimate, the E-value
        for the CI bound closest to the null, and an interpretation string.
    """
    from .iv import compute_e_value

    ci_lo = estimate - 1.96 * se
    ci_hi = estimate + 1.96 * se
    ev = compute_e_value(estimate, ci_lo, ci_hi, log_scale=log_scale)
    e_val = ev["e_value_point"]
    e_val_ci = ev.get("e_value_ci", 1.0)

    interpretation = (
        f"E-value = {e_val:.2f}. An unmeasured confounder would need to be "
        f"associated with both treatment and outcome by a risk-ratio factor "
        f"of {e_val:.2f} (above-and-beyond measured covariates) to fully "
        f"explain the observed effect. E-value for the 95% CI bound closest "
        f"to the null = {e_val_ci:.2f}."
    )

    return SensitivityResult(
        analysis_type="e_value",
        original_estimate=float(estimate),
        original_se=float(se),
        robustness_value=float(e_val),
        breakdown_point=float(e_val_ci),
        sensitivity_plot_data={
            "e_value_point": [float(e_val)],
            "e_value_ci": [float(e_val_ci)],
        },
        interpretation=interpretation,
    )


def omitted_variable_bias(
    data: SimulationData,
    treatment_name: str = "ecDNA_count",
    outcome_name: str = "VEGF_secretion",
    covariates: list[str] | None = None,
    r2_values: list[float] | None = None,
) -> SensitivityResult:
    """Assess potential bias from omitted variables.

    Based on Cinelli & Hazlett (2020) framework for omitted variable bias.

    Args:
        data: Loaded simulation data
        treatment_name: Treatment variable
        outcome_name: Outcome variable
        covariates: Controlled covariates
        r2_values: Hypothetical R² values for unmeasured confounder

    Returns:
        SensitivityResult with bias analysis.
    """
    if covariates is None:
        covariates = ["O2_local", "glucose_local"]

    if r2_values is None:
        r2_values = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3]

    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 20:
        return SensitivityResult(
            analysis_type="omitted_variable_bias",
            original_estimate=0.0,
            original_se=0.0,
            robustness_value=0.0,
            breakdown_point=0.0,
            sensitivity_plot_data={},
            interpretation="Insufficient data",
        )

    # Get data
    all_vars = [treatment_name, outcome_name] + covariates
    X = cells_to_array(tumor_cells, all_vars)
    n = len(tumor_cells)

    treatment = X[:, 0]
    outcome = X[:, 1]
    controls = X[:, 2:]

    # Controlled regression
    X_full = np.column_stack([np.ones(n), treatment, controls])
    beta_full = np.linalg.lstsq(X_full, outcome, rcond=None)[0]
    y_hat = X_full @ beta_full
    resid = outcome - y_hat

    original_effect = beta_full[1]

    # R² of full model
    ss_res = np.sum(resid ** 2)
    ss_tot = np.sum((outcome - np.mean(outcome)) ** 2)
    r2_full = 1 - ss_res / ss_tot

    # Compute standard error
    mse = ss_res / (n - X_full.shape[1])
    var_beta = mse * np.linalg.pinv(X_full.T @ X_full)
    se = np.sqrt(var_beta[1, 1])

    # Residual variation in treatment after controlling
    X_treat_ctrl = np.column_stack([np.ones(n), controls])
    beta_t = np.linalg.lstsq(X_treat_ctrl, treatment, rcond=None)[0]
    resid_t = treatment - X_treat_ctrl @ beta_t
    var_t_resid = np.var(resid_t)

    # Compute bias for each hypothetical confounder strength
    biases = []
    adjusted_effects = []

    for r2_yd in r2_values:  # R² of confounder with outcome
        for r2_td in r2_values:  # R² of confounder with treatment
            # Bias approximation (simplified)
            # bias ≈ sqrt(r2_yd) * sqrt(r2_td) * sd(Y) / sd(T_resid)
            bias = np.sqrt(r2_yd) * np.sqrt(r2_td) * np.std(outcome) / np.sqrt(var_t_resid)
            biases.append(bias)
            adjusted_effects.append(original_effect - bias)

    # Find R² combination that would nullify the effect
    breakdown_r2 = 0.0
    for r2 in r2_values:
        bias = r2 * np.std(outcome) / np.sqrt(var_t_resid)
        if bias >= abs(original_effect):
            breakdown_r2 = r2
            break
    else:
        breakdown_r2 = r2_values[-1]

    # Robustness value: what proportion of residual variance would confounder need?
    rv = (original_effect * np.sqrt(var_t_resid) / np.std(outcome)) ** 2

    interpretation = (
        f"Original effect: {original_effect:.4f} (SE: {se:.4f}). "
        f"To explain away the effect, an omitted variable would need to explain "
        f"at least {breakdown_r2*100:.1f}% of variation in both treatment and outcome. "
        f"Robustness value (RV) = {rv:.3f}."
    )

    return SensitivityResult(
        analysis_type="omitted_variable_bias",
        original_estimate=float(original_effect),
        original_se=float(se),
        robustness_value=float(rv),
        breakdown_point=float(breakdown_r2),
        sensitivity_plot_data={
            "r2_values": r2_values,
            "sample_biases": biases[:len(r2_values)],
        },
        interpretation=interpretation,
    )


def placebo_test(
    data: SimulationData,
    treatment_name: str = "ecDNA_count",
    outcome_name: str = "VEGF_secretion",
    placebo_outcomes: list[str] | None = None,
) -> dict[str, Any]:
    """Placebo/falsification tests using outcomes that shouldn't be affected.

    If treatment affects "placebo" outcomes, this suggests confounding.

    Args:
        data: Loaded simulation data
        treatment_name: Treatment variable
        outcome_name: True outcome
        placebo_outcomes: Variables that shouldn't be affected by treatment

    Returns:
        Dictionary with placebo test results.
    """
    if placebo_outcomes is None:
        # Use position as placebo (ecDNA shouldn't directly cause position)
        placebo_outcomes = ["x", "y"]

    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 20:
        return {"error": "insufficient_data"}

    results = {"treatment": treatment_name, "tests": []}

    # True effect
    treatment = np.array([c[treatment_name] for c in tumor_cells], dtype=float)
    outcome = np.array([c[outcome_name] for c in tumor_cells], dtype=float)

    n = len(treatment)
    X_true = np.column_stack([np.ones(n), treatment])
    beta_true = np.linalg.lstsq(X_true, outcome, rcond=None)[0]

    results["true_effect"] = {
        "outcome": outcome_name,
        "coefficient": float(beta_true[1]),
    }

    # Placebo tests
    for placebo in placebo_outcomes:
        if placebo not in tumor_cells[0]:
            continue

        placebo_outcome = np.array([c[placebo] for c in tumor_cells], dtype=float)
        beta_placebo = np.linalg.lstsq(X_true, placebo_outcome, rcond=None)[0]

        resid = placebo_outcome - X_true @ beta_placebo
        mse = np.sum(resid ** 2) / (n - 2)
        se = np.sqrt(mse * np.linalg.pinv(X_true.T @ X_true)[1, 1])

        t_stat = beta_placebo[1] / se if se > 0 else 0
        p_val = 2 * (1 - stats.t.cdf(abs(t_stat), n - 2))

        results["tests"].append({
            "placebo_outcome": placebo,
            "coefficient": float(beta_placebo[1]),
            "se": float(se),
            "p_value": float(p_val),
            "significant": p_val < 0.05,
        })

    # Interpretation
    n_significant = sum(1 for t in results["tests"] if t["significant"])
    if n_significant == 0:
        results["interpretation"] = "No placebo effects detected. Results appear robust."
    else:
        results["interpretation"] = (
            f"{n_significant} of {len(results['tests'])} placebo tests significant. "
            "This may indicate confounding or model misspecification."
        )

    return results
