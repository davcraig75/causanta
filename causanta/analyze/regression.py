"""Regression adjustment for causal inference.

Implements regression-based methods for estimating causal effects
while controlling for observed confounders.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import stats

from .loader import SimulationData, cells_to_array


@dataclass
class RegressionResult:
    """Result of regression-based causal analysis."""

    outcome_name: str
    treatment_name: str
    treatment_coef: float
    treatment_se: float
    treatment_pvalue: float
    treatment_ci_lower: float
    treatment_ci_upper: float
    r_squared: float
    adj_r_squared: float
    n_observations: int
    covariate_coefs: dict[str, float]
    covariate_ses: dict[str, float]
    model_type: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "outcome": self.outcome_name,
            "treatment": self.treatment_name,
            "treatment_coef": self.treatment_coef,
            "treatment_se": self.treatment_se,
            "treatment_pvalue": self.treatment_pvalue,
            "treatment_ci": [self.treatment_ci_lower, self.treatment_ci_upper],
            "r_squared": self.r_squared,
            "adj_r_squared": self.adj_r_squared,
            "n_observations": self.n_observations,
            "covariates": self.covariate_coefs,
            "model_type": self.model_type,
        }


def ols_with_inference(
    X: np.ndarray,
    y: np.ndarray,
    feature_names: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    """OLS regression with standard errors, p-values, and R-squared.

    Args:
        X: Design matrix with intercept (n x p)
        y: Response vector (n,)
        feature_names: Names for each feature (excluding intercept)

    Returns:
        (coefficients, standard_errors, p_values, r_squared, adj_r_squared)
    """
    n, p = X.shape

    # OLS estimate
    XtX = X.T @ X
    XtX_inv = np.linalg.pinv(XtX)
    beta = XtX_inv @ X.T @ y

    # Residuals
    y_hat = X @ beta
    residuals = y - y_hat

    # MSE and standard errors
    dof = n - p
    mse = np.sum(residuals ** 2) / dof if dof > 0 else 0

    var_beta = mse * XtX_inv
    se = np.sqrt(np.diag(var_beta))

    # T-statistics and p-values
    t_stats = beta / (se + 1e-10)
    p_values = 2 * (1 - stats.t.cdf(np.abs(t_stats), dof))

    # R-squared
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    # Adjusted R-squared
    adj_r_sq = 1 - (1 - r_sq) * (n - 1) / (n - p) if n > p else r_sq

    return beta, se, p_values, r_sq, adj_r_sq


def linear_regression_adjustment(
    data: SimulationData,
    outcome_name: str = "VEGF_secretion",
    treatment_name: str = "ecDNA_count",
    covariates: list[str] | None = None,
) -> RegressionResult:
    """Estimate treatment effect using linear regression adjustment.

    Y = b0 + b1*Treatment + b2*X1 + ... + error

    The treatment coefficient b1 is interpreted as the causal effect
    if all confounders are controlled.

    Args:
        data: Loaded simulation data
        outcome_name: Outcome variable name
        treatment_name: Treatment variable name
        covariates: List of covariate names to control for

    Returns:
        RegressionResult with effect estimates.
    """
    if covariates is None:
        covariates = ["O2_local", "glucose_local"]

    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 10:
        return RegressionResult(
            outcome_name=outcome_name,
            treatment_name=treatment_name,
            treatment_coef=0.0,
            treatment_se=0.0,
            treatment_pvalue=1.0,
            treatment_ci_lower=0.0,
            treatment_ci_upper=0.0,
            r_squared=0.0,
            adj_r_squared=0.0,
            n_observations=0,
            covariate_coefs={},
            covariate_ses={},
            model_type="insufficient_data",
        )

    # Build arrays
    all_vars = [treatment_name] + covariates
    X_raw = cells_to_array(tumor_cells, all_vars)
    y = np.array([c[outcome_name] for c in tumor_cells], dtype=float)

    n = len(y)

    # Add intercept
    X = np.column_stack([np.ones(n), X_raw])
    feature_names = ["intercept"] + all_vars

    # Fit model
    beta, se, pvalues, r_sq, adj_r_sq = ols_with_inference(X, y, feature_names)

    # 95% CI for treatment effect
    t_crit = stats.t.ppf(0.975, n - X.shape[1])
    ci_lower = beta[1] - t_crit * se[1]
    ci_upper = beta[1] + t_crit * se[1]

    # Build covariate results
    cov_coefs = {}
    cov_ses = {}
    for i, cov in enumerate(covariates):
        cov_coefs[cov] = float(beta[2 + i])
        cov_ses[cov] = float(se[2 + i])

    return RegressionResult(
        outcome_name=outcome_name,
        treatment_name=treatment_name,
        treatment_coef=float(beta[1]),
        treatment_se=float(se[1]),
        treatment_pvalue=float(pvalues[1]),
        treatment_ci_lower=float(ci_lower),
        treatment_ci_upper=float(ci_upper),
        r_squared=float(r_sq),
        adj_r_squared=float(adj_r_sq),
        n_observations=n,
        covariate_coefs=cov_coefs,
        covariate_ses=cov_ses,
        model_type="linear_regression",
    )


def polynomial_regression(
    data: SimulationData,
    outcome_name: str = "VEGF_secretion",
    treatment_name: str = "ecDNA_count",
    covariates: list[str] | None = None,
    degree: int = 2,
) -> RegressionResult:
    """Estimate treatment effect with polynomial terms.

    Allows for non-linear relationships between treatment and outcome.

    Args:
        data: Loaded simulation data
        outcome_name: Outcome variable name
        treatment_name: Treatment variable name
        covariates: Covariate names
        degree: Polynomial degree for treatment

    Returns:
        RegressionResult with effect estimates.
    """
    if covariates is None:
        covariates = ["O2_local", "glucose_local"]

    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 10:
        return RegressionResult(
            outcome_name=outcome_name,
            treatment_name=treatment_name,
            treatment_coef=0.0,
            treatment_se=0.0,
            treatment_pvalue=1.0,
            treatment_ci_lower=0.0,
            treatment_ci_upper=0.0,
            r_squared=0.0,
            adj_r_squared=0.0,
            n_observations=0,
            covariate_coefs={},
            covariate_ses={},
            model_type="insufficient_data",
        )

    # Build arrays
    treatment = np.array([c[treatment_name] for c in tumor_cells], dtype=float)
    cov_data = cells_to_array(tumor_cells, covariates)
    y = np.array([c[outcome_name] for c in tumor_cells], dtype=float)

    n = len(y)

    # Build design matrix with polynomial terms
    X_list = [np.ones(n)]
    feature_names = ["intercept"]

    for d in range(1, degree + 1):
        X_list.append(treatment ** d)
        feature_names.append(f"{treatment_name}^{d}" if d > 1 else treatment_name)

    X_list.append(cov_data)
    feature_names.extend(covariates)

    X = np.column_stack(X_list)

    # Fit model
    beta, se, pvalues, r_sq, adj_r_sq = ols_with_inference(X, y, feature_names)

    # Treatment effect is sum of polynomial coefficients (approximate marginal effect at mean)
    # For simplicity, report linear coefficient
    t_crit = stats.t.ppf(0.975, n - X.shape[1])
    ci_lower = beta[1] - t_crit * se[1]
    ci_upper = beta[1] + t_crit * se[1]

    cov_coefs = {}
    cov_ses = {}
    cov_start_idx = 1 + degree
    for i, cov in enumerate(covariates):
        cov_coefs[cov] = float(beta[cov_start_idx + i])
        cov_ses[cov] = float(se[cov_start_idx + i])

    return RegressionResult(
        outcome_name=outcome_name,
        treatment_name=treatment_name,
        treatment_coef=float(beta[1]),
        treatment_se=float(se[1]),
        treatment_pvalue=float(pvalues[1]),
        treatment_ci_lower=float(ci_lower),
        treatment_ci_upper=float(ci_upper),
        r_squared=float(r_sq),
        adj_r_squared=float(adj_r_sq),
        n_observations=n,
        covariate_coefs=cov_coefs,
        covariate_ses=cov_ses,
        model_type=f"polynomial_degree_{degree}",
    )


def interaction_regression(
    data: SimulationData,
    outcome_name: str = "VEGF_secretion",
    treatment_name: str = "ecDNA_count",
    modifier_name: str = "O2_local",
    covariates: list[str] | None = None,
) -> dict[str, Any]:
    """Estimate heterogeneous treatment effects with interaction terms.

    Y = b0 + b1*T + b2*M + b3*T*M + ...

    Args:
        data: Loaded simulation data
        outcome_name: Outcome variable
        treatment_name: Treatment variable
        modifier_name: Effect modifier variable
        covariates: Additional covariates

    Returns:
        Dictionary with main and interaction effects.
    """
    if covariates is None:
        covariates = ["glucose_local"]

    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 10:
        return {"error": "insufficient_data"}

    treatment = np.array([c[treatment_name] for c in tumor_cells], dtype=float)
    modifier = np.array([c[modifier_name] for c in tumor_cells], dtype=float)
    y = np.array([c[outcome_name] for c in tumor_cells], dtype=float)

    # Covariates excluding modifier
    covariates = [c for c in covariates if c != modifier_name]
    cov_data = cells_to_array(tumor_cells, covariates) if covariates else np.zeros((len(y), 0))

    n = len(y)

    # Build design matrix
    interaction = treatment * modifier
    X = np.column_stack([
        np.ones(n),
        treatment,
        modifier,
        interaction,
        cov_data,
    ])

    feature_names = ["intercept", treatment_name, modifier_name,
                     f"{treatment_name}*{modifier_name}"] + covariates

    # Fit
    beta, se, pvalues, r_sq, adj_r_sq = ols_with_inference(X, y, feature_names)

    return {
        "outcome": outcome_name,
        "treatment": treatment_name,
        "modifier": modifier_name,
        "main_effect": float(beta[1]),
        "main_effect_se": float(se[1]),
        "main_effect_pvalue": float(pvalues[1]),
        "modifier_effect": float(beta[2]),
        "modifier_effect_se": float(se[2]),
        "interaction_effect": float(beta[3]),
        "interaction_effect_se": float(se[3]),
        "interaction_pvalue": float(pvalues[3]),
        "r_squared": float(r_sq),
        "n_observations": n,
        "interpretation": (
            f"At mean {modifier_name}, effect of {treatment_name} on {outcome_name} is "
            f"{beta[1]:.4f} + {beta[3]:.4f}*{modifier_name}"
        ),
    }


def doubly_robust_estimation(
    data: SimulationData,
    treatment_threshold: float = 10.0,
    outcome_name: str = "VEGF_secretion",
    covariates: list[str] | None = None,
) -> dict[str, Any]:
    """Doubly robust estimator combining regression and propensity weighting.

    Provides consistent estimates if either the outcome model or propensity
    model is correctly specified.

    Args:
        data: Loaded simulation data
        treatment_threshold: ecDNA threshold for binary treatment
        outcome_name: Outcome variable
        covariates: Covariate names

    Returns:
        Dictionary with doubly robust estimates.
    """
    if covariates is None:
        covariates = ["O2_local", "glucose_local"]

    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 20:
        return {"error": "insufficient_data"}

    # Build arrays
    ecDNA = np.array([c["ecDNA_count"] for c in tumor_cells], dtype=float)
    treatment = (ecDNA > treatment_threshold).astype(float)
    outcome = np.array([c[outcome_name] for c in tumor_cells], dtype=float)
    X = cells_to_array(tumor_cells, covariates)

    n = len(outcome)

    # Import propensity function
    from .matching import logistic_propensity

    # Estimate propensity scores
    propensity = logistic_propensity(X, treatment)

    # Outcome models (separate for treated and control)
    X_int = np.column_stack([np.ones(n), X])

    # Model for treated
    treated_idx = treatment == 1
    n_treated = np.sum(treated_idx)
    if n_treated > 2:
        beta_t = np.linalg.lstsq(X_int[treated_idx], outcome[treated_idx], rcond=None)[0]
        mu_1 = X_int @ beta_t  # Predicted outcome if treated
    elif n_treated > 0:
        mu_1 = np.mean(outcome[treated_idx]) * np.ones(n)
    else:
        mu_1 = np.zeros(n)

    # Model for control
    control_idx = treatment == 0
    n_control = np.sum(control_idx)
    if n_control > 2:
        beta_c = np.linalg.lstsq(X_int[control_idx], outcome[control_idx], rcond=None)[0]
        mu_0 = X_int @ beta_c  # Predicted outcome if control
    elif n_control > 0:
        mu_0 = np.mean(outcome[control_idx]) * np.ones(n)
    else:
        mu_0 = np.zeros(n)

    # Doubly robust estimator (clip propensities to avoid division issues)
    propensity_clipped = np.clip(propensity, 0.01, 0.99)
    dr_1 = mu_1 + treatment * (outcome - mu_1) / propensity_clipped
    dr_0 = mu_0 + (1 - treatment) * (outcome - mu_0) / (1 - propensity_clipped)

    ate_dr = np.mean(dr_1 - dr_0)

    # Bootstrap SE
    n_boot = 200
    boot_ates = []

    for _ in range(n_boot):
        idx = np.random.choice(n, size=n, replace=True)
        boot_ates.append(np.mean(dr_1[idx] - dr_0[idx]))

    ate_se = np.std(boot_ates)

    return {
        "method": "doubly_robust",
        "treatment": f"ecDNA>{treatment_threshold}",
        "outcome": outcome_name,
        "ate": float(ate_dr),
        "ate_se": float(ate_se),
        "ate_ci_lower": float(ate_dr - 1.96 * ate_se),
        "ate_ci_upper": float(ate_dr + 1.96 * ate_se),
        "n_observations": n,
        "n_treated": int(np.sum(treatment)),
        "n_control": int(np.sum(1 - treatment)),
    }
