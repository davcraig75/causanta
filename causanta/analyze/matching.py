"""Propensity score matching for causal inference.

Implements propensity score estimation and matching methods to
estimate causal effects while controlling for observed confounders.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import stats

from .loader import SimulationData, cells_to_array


@dataclass
class MatchingResult:
    """Result of propensity score matching analysis."""

    treatment_name: str
    outcome_name: str
    n_treated: int
    n_control: int
    n_matched: int
    ate: float  # Average Treatment Effect
    att: float  # Average Treatment Effect on Treated
    ate_se: float
    att_se: float
    balance_before: dict[str, float]  # Standardized mean differences before
    balance_after: dict[str, float]  # Standardized mean differences after
    caliper: float
    method: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "treatment": self.treatment_name,
            "outcome": self.outcome_name,
            "n_treated": self.n_treated,
            "n_control": self.n_control,
            "n_matched": self.n_matched,
            "ate": self.ate,
            "att": self.att,
            "ate_se": self.ate_se,
            "att_se": self.att_se,
            "balance_before": self.balance_before,
            "balance_after": self.balance_after,
            "caliper": self.caliper,
            "method": self.method,
        }


def logistic_propensity(
    X: np.ndarray,
    treatment: np.ndarray,
    max_iter: int = 100,
) -> np.ndarray:
    """Estimate propensity scores using logistic regression.

    Uses gradient descent for simplicity (no external dependencies).

    Args:
        X: Covariate matrix (n x p)
        treatment: Binary treatment indicator (n,)
        max_iter: Maximum iterations for optimization

    Returns:
        Propensity scores (n,)
    """
    n, p = X.shape

    # Add intercept
    X_int = np.column_stack([np.ones(n), X])

    # Initialize coefficients
    beta = np.zeros(p + 1)

    # Learning rate
    lr = 0.1

    for _ in range(max_iter):
        # Predictions
        z = X_int @ beta
        # Clip to avoid overflow
        z = np.clip(z, -20, 20)
        p_hat = 1 / (1 + np.exp(-z))

        # Gradient
        gradient = X_int.T @ (p_hat - treatment) / n

        # Update
        beta -= lr * gradient

        # Check convergence
        if np.max(np.abs(gradient)) < 1e-6:
            break

    # Final propensity scores
    z = X_int @ beta
    z = np.clip(z, -20, 20)
    propensity = 1 / (1 + np.exp(-z))

    # Clip to avoid exact 0 or 1
    return np.clip(propensity, 0.01, 0.99)


def nearest_neighbor_matching(
    propensity: np.ndarray,
    treatment: np.ndarray,
    caliper: float = 0.1,
    replace: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Match treated units to control units based on propensity score.

    Args:
        propensity: Propensity scores (n,)
        treatment: Binary treatment indicator (n,)
        caliper: Maximum allowed distance for matching
        replace: Whether to match with replacement

    Returns:
        (treated_indices, matched_control_indices)
    """
    treated_idx = np.where(treatment == 1)[0]
    control_idx = np.where(treatment == 0)[0]

    if len(control_idx) == 0:
        return np.array([]), np.array([])

    matched_treated = []
    matched_control = []
    used_controls = set()

    for t_idx in treated_idx:
        t_prop = propensity[t_idx]

        # Find closest control
        best_dist = float('inf')
        best_c_idx = None

        for c_idx in control_idx:
            if not replace and c_idx in used_controls:
                continue

            dist = abs(propensity[c_idx] - t_prop)
            if dist < best_dist and dist <= caliper:
                best_dist = dist
                best_c_idx = c_idx

        if best_c_idx is not None:
            matched_treated.append(t_idx)
            matched_control.append(best_c_idx)
            used_controls.add(best_c_idx)

    return np.array(matched_treated), np.array(matched_control)


def standardized_mean_diff(
    X: np.ndarray,
    treatment: np.ndarray,
    covariate_names: list[str],
) -> dict[str, float]:
    """Compute standardized mean differences for covariates.

    SMD = (mean_treated - mean_control) / pooled_sd

    Args:
        X: Covariate matrix (n x p)
        treatment: Binary treatment indicator
        covariate_names: Names for each covariate

    Returns:
        Dictionary mapping covariate name to SMD.
    """
    treated = treatment == 1
    control = treatment == 0
    n_treated = np.sum(treated)
    n_control = np.sum(control)

    result = {}
    for j, name in enumerate(covariate_names):
        # Guard against empty groups
        if n_treated == 0 or n_control == 0:
            result[name] = 0.0
            continue

        mean_t = np.mean(X[treated, j])
        mean_c = np.mean(X[control, j])

        var_t = np.var(X[treated, j], ddof=1) if n_treated > 1 else 0.0
        var_c = np.var(X[control, j], ddof=1) if n_control > 1 else 0.0

        pooled_sd = np.sqrt((var_t + var_c) / 2)

        if pooled_sd > 0:
            smd = (mean_t - mean_c) / pooled_sd
        else:
            smd = 0.0

        result[name] = smd

    return result


def propensity_score_matching(
    data: SimulationData,
    treatment_threshold: float = 10.0,
    outcome_name: str = "VEGF_secretion",
    covariates: list[str] | None = None,
    caliper: float = 0.1,
) -> MatchingResult:
    """Perform propensity score matching analysis.

    Treats high-ecDNA cells (above threshold) as "treated" and estimates
    the causal effect on the outcome.

    Args:
        data: Loaded simulation data
        treatment_threshold: ecDNA count threshold for treatment
        outcome_name: Outcome variable
        covariates: Covariate names (default: O2, glucose)
        caliper: Matching caliper

    Returns:
        MatchingResult with effect estimates.
    """
    if covariates is None:
        covariates = ["O2_local", "glucose_local"]

    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 20:
        return MatchingResult(
            treatment_name=f"ecDNA>{treatment_threshold}",
            outcome_name=outcome_name,
            n_treated=0,
            n_control=0,
            n_matched=0,
            ate=0.0,
            att=0.0,
            ate_se=0.0,
            att_se=0.0,
            balance_before={},
            balance_after={},
            caliper=caliper,
            method="insufficient_data",
        )

    # Build arrays
    ecDNA = np.array([c["ecDNA_count"] for c in tumor_cells], dtype=float)
    treatment = (ecDNA > treatment_threshold).astype(float)
    outcome = np.array([c[outcome_name] for c in tumor_cells], dtype=float)
    X = cells_to_array(tumor_cells, covariates)

    n_treated = int(np.sum(treatment))
    n_control = len(treatment) - n_treated

    if n_treated < 5 or n_control < 5:
        return MatchingResult(
            treatment_name=f"ecDNA>{treatment_threshold}",
            outcome_name=outcome_name,
            n_treated=n_treated,
            n_control=n_control,
            n_matched=0,
            ate=0.0,
            att=0.0,
            ate_se=0.0,
            att_se=0.0,
            balance_before={},
            balance_after={},
            caliper=caliper,
            method="insufficient_treated_or_control",
        )

    # Compute balance before matching
    balance_before = standardized_mean_diff(X, treatment, covariates)

    # Estimate propensity scores
    propensity = logistic_propensity(X, treatment)

    # Match
    treated_idx, control_idx = nearest_neighbor_matching(
        propensity, treatment, caliper=caliper
    )

    n_matched = len(treated_idx)

    if n_matched < 5:
        return MatchingResult(
            treatment_name=f"ecDNA>{treatment_threshold}",
            outcome_name=outcome_name,
            n_treated=n_treated,
            n_control=n_control,
            n_matched=n_matched,
            ate=0.0,
            att=0.0,
            ate_se=0.0,
            att_se=0.0,
            balance_before=balance_before,
            balance_after={},
            caliper=caliper,
            method="insufficient_matches",
        )

    # Compute balance after matching
    matched_treatment = np.zeros(n_matched * 2)
    matched_treatment[:n_matched] = 1
    matched_X = np.vstack([X[treated_idx], X[control_idx]])
    balance_after = standardized_mean_diff(matched_X, matched_treatment, covariates)

    # Estimate effects
    Y_treated = outcome[treated_idx]
    Y_control = outcome[control_idx]

    # ATT: Average Treatment Effect on Treated
    att = np.mean(Y_treated - Y_control)
    att_se = np.std(Y_treated - Y_control, ddof=1) / np.sqrt(n_matched)

    # ATE: Using all matched pairs
    ate = att  # With 1-1 matching, ATE = ATT
    ate_se = att_se

    return MatchingResult(
        treatment_name=f"ecDNA>{treatment_threshold}",
        outcome_name=outcome_name,
        n_treated=n_treated,
        n_control=n_control,
        n_matched=n_matched,
        ate=ate,
        att=att,
        ate_se=ate_se,
        att_se=att_se,
        balance_before=balance_before,
        balance_after=balance_after,
        caliper=caliper,
        method="nearest_neighbor",
    )


def inverse_propensity_weighting(
    data: SimulationData,
    treatment_threshold: float = 10.0,
    outcome_name: str = "VEGF_secretion",
    covariates: list[str] | None = None,
) -> dict[str, Any]:
    """Estimate effects using inverse propensity weighting (IPW).

    Args:
        data: Loaded simulation data
        treatment_threshold: ecDNA threshold for treatment
        outcome_name: Outcome variable
        covariates: Covariate names

    Returns:
        Dictionary with IPW estimates.
    """
    if covariates is None:
        covariates = ["O2_local", "glucose_local"]

    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 20:
        return {"error": "insufficient_data"}

    ecDNA = np.array([c["ecDNA_count"] for c in tumor_cells], dtype=float)
    treatment = (ecDNA > treatment_threshold).astype(float)
    outcome = np.array([c[outcome_name] for c in tumor_cells], dtype=float)
    X = cells_to_array(tumor_cells, covariates)

    # Estimate propensity scores
    propensity = logistic_propensity(X, treatment)

    # IPW weights
    # For treated: 1/p(X)
    # For control: 1/(1-p(X))
    # Clip propensities to avoid division by zero
    propensity_clipped = np.clip(propensity, 0.01, 0.99)
    weights_treated = treatment / propensity_clipped
    weights_control = (1 - treatment) / (1 - propensity_clipped)

    # Normalize weights (with guards against empty groups)
    sum_wt = np.sum(weights_treated)
    sum_wc = np.sum(weights_control)
    if sum_wt > 0:
        weights_treated = weights_treated / sum_wt * np.sum(treatment)
    if sum_wc > 0:
        weights_control = weights_control / sum_wc * np.sum(1 - treatment)

    # IPW estimate of ATE
    n_treated = np.sum(treatment)
    n_control = np.sum(1 - treatment)
    if n_treated > 0 and n_control > 0:
        ate_ipw = np.sum(weights_treated * outcome) / n_treated - \
                  np.sum(weights_control * outcome) / n_control
    else:
        ate_ipw = 0.0

    # Bootstrap SE
    n_boot = 200
    boot_ates = []
    n = len(outcome)

    for _ in range(n_boot):
        idx = np.random.choice(n, size=n, replace=True)
        t_b = treatment[idx]
        y_b = outcome[idx]
        p_b = np.clip(propensity_clipped[idx], 0.01, 0.99)

        w_t = t_b / p_b
        w_c = (1 - t_b) / (1 - p_b)

        if np.sum(t_b) > 0 and np.sum(1 - t_b) > 0:
            ate_b = np.sum(w_t * y_b) / np.sum(t_b) - np.sum(w_c * y_b) / np.sum(1 - t_b)
            boot_ates.append(ate_b)

    ate_se = np.std(boot_ates) if boot_ates else 0.0

    # Compute mean propensities (with guards for empty groups)
    mean_prop_treated = float(np.mean(propensity[treatment == 1])) if n_treated > 0 else 0.0
    mean_prop_control = float(np.mean(propensity[treatment == 0])) if n_control > 0 else 0.0

    return {
        "method": "IPW",
        "treatment": f"ecDNA>{treatment_threshold}",
        "outcome": outcome_name,
        "ate": float(ate_ipw),
        "ate_se": float(ate_se),
        "n_treated": int(n_treated),
        "n_control": int(n_control),
        "mean_propensity_treated": mean_prop_treated,
        "mean_propensity_control": mean_prop_control,
    }
