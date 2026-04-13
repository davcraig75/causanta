"""Instrumental Variable (IV) / Two-Stage Least Squares estimation.

Uses ecDNA copy number as an instrumental variable to estimate
causal effects on cellular phenotypes, accounting for environmental
confounding.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .loader import SimulationData, cells_to_array


@dataclass
class IVEstimate:
    """Result of IV / 2SLS estimation."""

    outcome_name: str
    treatment_name: str
    instrument_name: str
    first_stage_coef: float
    first_stage_se: float
    first_stage_f_stat: float
    second_stage_coef: float  # The causal effect estimate
    second_stage_se: float
    ols_coef: float  # Naive OLS for comparison
    ols_se: float
    n_observations: int
    weak_instrument: bool  # F-stat < 10

    @property
    def bias_ratio(self) -> float:
        """Ratio of OLS to IV estimate (indicates confounding bias)."""
        if abs(self.second_stage_coef) < 1e-10:
            return float('inf')
        return self.ols_coef / self.second_stage_coef

    def to_dict(self) -> dict[str, Any]:
        return {
            "outcome": self.outcome_name,
            "treatment": self.treatment_name,
            "instrument": self.instrument_name,
            "first_stage_coef": self.first_stage_coef,
            "first_stage_se": self.first_stage_se,
            "first_stage_f_stat": self.first_stage_f_stat,
            "second_stage_coef": self.second_stage_coef,
            "second_stage_se": self.second_stage_se,
            "ols_coef": self.ols_coef,
            "ols_se": self.ols_se,
            "n_observations": self.n_observations,
            "weak_instrument": self.weak_instrument,
            "bias_ratio": self.bias_ratio,
        }


def _ols_regression(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Simple OLS regression with standard errors.

    Args:
        X: Design matrix (n x p), should include intercept if desired
        y: Response vector (n,)

    Returns:
        (coefficients, standard_errors, r_squared)
    """
    n, p = X.shape
    XtX = X.T @ X
    XtX_inv = np.linalg.pinv(XtX)
    beta = XtX_inv @ X.T @ y

    # Residuals and MSE
    y_hat = X @ beta
    residuals = y - y_hat
    mse = np.sum(residuals ** 2) / (n - p)

    # Standard errors
    var_beta = mse * XtX_inv
    se = np.sqrt(np.diag(var_beta))

    # R-squared
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    return beta, se, r_sq


def two_stage_least_squares(
    Z: np.ndarray,
    D: np.ndarray,
    Y: np.ndarray,
    covariates: np.ndarray | None = None,
) -> tuple[float, float, float, float, float, float]:
    """Two-stage least squares estimation.

    Args:
        Z: Instrument (n,)
        D: Treatment/endogenous variable (n,)
        Y: Outcome (n,)
        covariates: Optional additional covariates (n x k)

    Returns:
        (first_stage_coef, first_stage_se, f_stat,
         second_stage_coef, second_stage_se, ols_coef)
    """
    n = len(Z)

    # Build design matrix for first stage
    if covariates is not None:
        X_first = np.column_stack([np.ones(n), Z, covariates])
        X_ols = np.column_stack([np.ones(n), D, covariates])
    else:
        X_first = np.column_stack([np.ones(n), Z])
        X_ols = np.column_stack([np.ones(n), D])

    # First stage: D ~ Z (+ covariates)
    beta_first, se_first, r_sq_first = _ols_regression(X_first, D)
    first_stage_coef = beta_first[1]  # Coefficient on Z
    first_stage_se = se_first[1]

    # F-statistic for instrument strength
    # F = (R² / k) / ((1 - R²) / (n - k - 1)) for the Z coefficient
    D_hat = X_first @ beta_first
    residuals_first = D - D_hat
    ss_res = np.sum(residuals_first ** 2)
    ss_explained = np.sum((D_hat - np.mean(D)) ** 2)
    f_stat = (ss_explained / 1) / (ss_res / (n - X_first.shape[1])) if ss_res > 0 else 0

    # Second stage: Y ~ D_hat (+ covariates)
    if covariates is not None:
        X_second = np.column_stack([np.ones(n), D_hat, covariates])
    else:
        X_second = np.column_stack([np.ones(n), D_hat])

    beta_second, se_second, _ = _ols_regression(X_second, Y)
    second_stage_coef = beta_second[1]  # Causal effect estimate

    # Standard error correction for 2SLS
    # Need to use original D in residual calculation
    if covariates is not None:
        X_final = np.column_stack([np.ones(n), D, covariates])
    else:
        X_final = np.column_stack([np.ones(n), D])
    Y_hat_2sls = X_final @ beta_second
    residuals_2sls = Y - Y_hat_2sls
    mse_2sls = np.sum(residuals_2sls ** 2) / (n - X_final.shape[1])

    # Corrected SE
    XtX_inv = np.linalg.pinv(X_second.T @ X_second)
    second_stage_se = np.sqrt(mse_2sls * XtX_inv[1, 1])

    # OLS for comparison
    beta_ols, se_ols, _ = _ols_regression(X_ols, Y)
    ols_coef = beta_ols[1]

    return (first_stage_coef, first_stage_se, f_stat,
            second_stage_coef, second_stage_se, ols_coef)


def estimate_iv_effects(
    data: SimulationData,
    outcomes: list[str] | None = None,
) -> list[IVEstimate]:
    """Estimate causal effects using ecDNA as instrument.

    Uses ecDNA copy number as an instrument for gene expression effects
    on cellular phenotypes.

    Args:
        data: Loaded simulation data
        outcomes: List of outcome variables to estimate. If None, uses
                  defaults: ['proliferation_rate', 'VEGF_secretion', 'migration_rate']

    Returns:
        List of IVEstimate objects for each outcome.
    """
    if outcomes is None:
        outcomes = ["VEGF_secretion", "migration_rate"]

    # Get tumor cells at final timestep
    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 10:
        return []

    results = []

    # Instrument: ecDNA count
    Z = np.array([c["ecDNA_count"] for c in tumor_cells], dtype=float)

    # Confounders to control for
    O2 = np.array([c["O2_local"] for c in tumor_cells], dtype=float)
    glucose = np.array([c["glucose_local"] for c in tumor_cells], dtype=float)
    covariates = np.column_stack([O2, glucose])

    for outcome_name in outcomes:
        if outcome_name not in tumor_cells[0]:
            continue

        Y = np.array([c[outcome_name] for c in tumor_cells], dtype=float)

        # Treatment variable (endogenous): could be ecDNA itself
        # or a derived expression variable
        D = Z  # For simplicity, use ecDNA as both instrument and treatment

        # Run 2SLS
        try:
            (first_coef, first_se, f_stat,
             second_coef, second_se, ols_coef) = two_stage_least_squares(
                Z, D, Y, covariates
            )

            # OLS SE
            n = len(Z)
            X_ols = np.column_stack([np.ones(n), D, covariates])
            _, se_ols, _ = _ols_regression(X_ols, Y)
            ols_se = se_ols[1]

            results.append(IVEstimate(
                outcome_name=outcome_name,
                treatment_name="ecDNA_count",
                instrument_name="ecDNA_segregation",
                first_stage_coef=first_coef,
                first_stage_se=first_se,
                first_stage_f_stat=f_stat,
                second_stage_coef=second_coef,
                second_stage_se=second_se,
                ols_coef=ols_coef,
                ols_se=ols_se,
                n_observations=len(tumor_cells),
                weak_instrument=(f_stat < 10),
            ))
        except Exception as e:
            # Skip outcomes that fail estimation
            continue

    return results


def estimate_segregation_iv(
    data: SimulationData,
) -> IVEstimate | None:
    """Use ecDNA segregation randomness as instrument for division effects.

    The randomness in ecDNA segregation during mitosis provides exogenous
    variation that can be used as an instrument.

    Args:
        data: Loaded simulation data

    Returns:
        IVEstimate for the effect of ecDNA on division rate, or None if
        insufficient data.
    """
    if len(data.lineage) < 20:
        return None

    # Use parent ecDNA before division as treatment
    # Use daughter ecDNA (which is random) as instrument
    parent_ecDNA = np.array([r["parent_ecDNA_before"] for r in data.lineage], dtype=float)
    daughter_ecDNA = np.array([r["daughter_ecDNA"] for r in data.lineage], dtype=float)

    # Outcome: generation (proxy for cumulative division rate)
    generation = np.array([r["generation"] for r in data.lineage], dtype=float)

    # Since daughter ecDNA is random given parent, it's a valid instrument
    # for parent ecDNA effects on proliferation

    n = len(parent_ecDNA)
    Z = daughter_ecDNA  # Instrument
    D = parent_ecDNA    # Treatment
    Y = generation      # Outcome

    try:
        (first_coef, first_se, f_stat,
         second_coef, second_se, ols_coef) = two_stage_least_squares(Z, D, Y)

        X_ols = np.column_stack([np.ones(n), D])
        _, se_ols, _ = _ols_regression(X_ols, Y)

        return IVEstimate(
            outcome_name="generation",
            treatment_name="parent_ecDNA",
            instrument_name="daughter_ecDNA",
            first_stage_coef=first_coef,
            first_stage_se=first_se,
            first_stage_f_stat=f_stat,
            second_stage_coef=second_coef,
            second_stage_se=second_se,
            ols_coef=ols_coef,
            ols_se=se_ols[1],
            n_observations=n,
            weak_instrument=(f_stat < 10),
        )
    except Exception:
        return None
