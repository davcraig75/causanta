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

        # Treatment variable (endogenous): EGFR expression (instrumented by ecDNA)
        # This is the correct IV setup: Z=ecDNA (instrument), D=EGFR (treatment)
        if "egfr_expression" in tumor_cells[0] and tumor_cells[0]["egfr_expression"] > 0:
            D = np.array([c["egfr_expression"] for c in tumor_cells], dtype=float)
            treatment_name = "egfr_expression"
        else:
            # Fallback: reduced-form analysis (ecDNA as both Z and D)
            # This collapses 2SLS to OLS — results should be interpreted with caution
            D = Z
            treatment_name = "ecDNA_count"

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
                treatment_name=treatment_name,
                instrument_name="ecDNA_count",
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


# =============================================================================
# Additional IV Diagnostics
# =============================================================================

@dataclass
class IVDiagnostics:
    """Comprehensive IV diagnostics."""

    # Instrument strength
    f_statistic: float
    first_stage_r2: float
    partial_r2: float  # R² of Z after partialing out covariates

    # Endogeneity test
    wu_hausman_stat: float
    wu_hausman_pvalue: float
    endogeneity_detected: bool

    # Weak instrument robust inference
    anderson_rubin_stat: float
    anderson_rubin_pvalue: float
    anderson_rubin_ci: tuple[float, float]

    # Stock-Yogo critical values for weak IV
    stock_yogo_5pct: float  # Critical value for 5% maximal IV bias
    stock_yogo_10pct: float  # Critical value for 10% maximal IV bias

    def to_dict(self) -> dict[str, Any]:
        return {
            "f_statistic": self.f_statistic,
            "first_stage_r2": self.first_stage_r2,
            "partial_r2": self.partial_r2,
            "wu_hausman_stat": self.wu_hausman_stat,
            "wu_hausman_pvalue": self.wu_hausman_pvalue,
            "endogeneity_detected": self.endogeneity_detected,
            "anderson_rubin_stat": self.anderson_rubin_stat,
            "anderson_rubin_pvalue": self.anderson_rubin_pvalue,
            "anderson_rubin_ci": self.anderson_rubin_ci,
            "stock_yogo_5pct": self.stock_yogo_5pct,
            "stock_yogo_10pct": self.stock_yogo_10pct,
        }


def compute_iv_diagnostics(
    Z: np.ndarray,
    D: np.ndarray,
    Y: np.ndarray,
    covariates: np.ndarray | None = None,
    alpha: float = 0.05,
) -> IVDiagnostics:
    """Compute comprehensive IV diagnostics.

    Args:
        Z: Instrument (n,)
        D: Treatment (n,)
        Y: Outcome (n,)
        covariates: Optional covariates (n x k)
        alpha: Significance level for tests

    Returns:
        IVDiagnostics object
    """
    from scipy import stats

    n = len(Z)

    # Build design matrices
    if covariates is not None:
        X_full = np.column_stack([np.ones(n), Z, covariates])
        X_partial = np.column_stack([np.ones(n), covariates])
    else:
        X_full = np.column_stack([np.ones(n), Z])
        X_partial = np.ones((n, 1))

    # First-stage regression
    beta_first, se_first, r2_first = _ols_regression(X_full, D)
    D_hat = X_full @ beta_first

    # Partial R² (contribution of Z after covariates)
    if covariates is not None:
        beta_partial, _, r2_partial_only = _ols_regression(X_partial, D)
        partial_r2 = r2_first - r2_partial_only
    else:
        partial_r2 = r2_first

    # F-statistic for instrument
    residuals_first = D - D_hat
    ss_res = np.sum(residuals_first ** 2)
    ss_explained = np.sum((D_hat - np.mean(D)) ** 2)
    k = 1  # Number of instruments
    f_stat = (ss_explained / k) / (ss_res / (n - X_full.shape[1])) if ss_res > 0 else 0

    # Wu-Hausman endogeneity test
    # Compare OLS and 2SLS estimates
    first_stage_resid = D - D_hat

    # Augmented regression: Y ~ D + v_hat
    if covariates is not None:
        X_augmented = np.column_stack([np.ones(n), D, covariates, first_stage_resid])
    else:
        X_augmented = np.column_stack([np.ones(n), D, first_stage_resid])

    beta_aug, se_aug, _ = _ols_regression(X_augmented, Y)
    # Test coefficient on first_stage_resid
    resid_coef_idx = -1
    t_stat_hausman = beta_aug[resid_coef_idx] / se_aug[resid_coef_idx] if se_aug[resid_coef_idx] > 0 else 0
    wu_hausman_pvalue = 2 * (1 - stats.t.cdf(abs(t_stat_hausman), n - X_augmented.shape[1]))

    # Anderson-Rubin test (robust to weak instruments)
    # Test H0: beta = 0 using reduced form
    if covariates is not None:
        X_rf = np.column_stack([np.ones(n), Z, covariates])
    else:
        X_rf = np.column_stack([np.ones(n), Z])

    beta_rf, se_rf, _ = _ols_regression(X_rf, Y)
    ar_stat = (beta_rf[1] / se_rf[1]) ** 2 if se_rf[1] > 0 else 0
    ar_pvalue = 1 - stats.chi2.cdf(ar_stat, 1)

    # Anderson-Rubin confidence interval
    # AR CI inverts the AR test
    # For simplicity, use normal approximation
    z_crit = stats.norm.ppf(1 - alpha / 2)
    if abs(beta_first[1]) > 1e-10:
        ar_point = beta_rf[1] / beta_first[1]
        # Delta method SE
        ar_se = abs(ar_point) * np.sqrt(
            (se_rf[1] / beta_rf[1]) ** 2 + (se_first[1] / beta_first[1]) ** 2
            if abs(beta_rf[1]) > 1e-10 else 1
        ) if abs(beta_rf[1]) > 1e-10 else float('inf')
        ar_ci = (ar_point - z_crit * ar_se, ar_point + z_crit * ar_se)
    else:
        ar_ci = (float('-inf'), float('inf'))

    # Stock-Yogo critical values (approximate for single instrument)
    # For one instrument, one endogenous regressor
    stock_yogo_5pct = 16.38   # 5% maximal IV size
    stock_yogo_10pct = 8.96   # 10% maximal IV size

    return IVDiagnostics(
        f_statistic=f_stat,
        first_stage_r2=r2_first,
        partial_r2=partial_r2,
        wu_hausman_stat=t_stat_hausman ** 2,
        wu_hausman_pvalue=wu_hausman_pvalue,
        endogeneity_detected=wu_hausman_pvalue < alpha,
        anderson_rubin_stat=ar_stat,
        anderson_rubin_pvalue=ar_pvalue,
        anderson_rubin_ci=ar_ci,
        stock_yogo_5pct=stock_yogo_5pct,
        stock_yogo_10pct=stock_yogo_10pct,
    )


def compute_e_value(
    point_estimate: float,
    ci_lower: float | None = None,
    ci_upper: float | None = None,
) -> dict[str, float]:
    """Compute E-value for sensitivity analysis.

    The E-value is the minimum strength of association that an unmeasured
    confounder would need with both treatment and outcome to explain away
    the observed effect.

    Args:
        point_estimate: Point estimate of causal effect (risk ratio scale)
        ci_lower: Lower confidence interval bound
        ci_upper: Upper confidence interval bound

    Returns:
        Dictionary with E-values for point estimate and CI bound
    """
    # Convert to risk ratio if needed (assume on RR scale)
    rr = max(point_estimate, 1 / point_estimate) if point_estimate > 0 else 1.0

    # E-value formula: E = RR + sqrt(RR * (RR - 1))
    e_point = rr + np.sqrt(rr * (rr - 1)) if rr > 1 else 1.0

    result = {"e_value_point": e_point}

    # E-value for CI bound closest to null
    if ci_lower is not None and ci_upper is not None:
        # Find bound closest to 1 (null)
        if ci_lower > 1:
            rr_ci = ci_lower
        elif ci_upper < 1:
            rr_ci = 1 / ci_upper
        else:
            # CI includes null
            rr_ci = 1.0

        e_ci = rr_ci + np.sqrt(rr_ci * (rr_ci - 1)) if rr_ci > 1 else 1.0
        result["e_value_ci"] = e_ci

    return result


def sensitivity_analysis_rosenbaum(
    Z: np.ndarray,
    D: np.ndarray,
    Y: np.ndarray,
    gamma_values: list[float] | None = None,
) -> dict[str, Any]:
    """Rosenbaum bounds sensitivity analysis.

    Computes how sensitive the IV results are to violations of
    the exclusion restriction (unobserved confounding).

    Args:
        Z: Instrument
        D: Treatment
        Y: Outcome
        gamma_values: Values of Gamma to evaluate

    Returns:
        Dictionary with sensitivity analysis results
    """
    from scipy import stats

    if gamma_values is None:
        gamma_values = [1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0]

    n = len(Z)

    # Basic IV estimate
    X = np.column_stack([np.ones(n), Z])
    beta_first, _, _ = _ols_regression(X, D)
    D_hat = X @ beta_first

    X_second = np.column_stack([np.ones(n), D_hat])
    beta_iv, se_iv, _ = _ols_regression(X_second, Y)

    iv_estimate = beta_iv[1]
    iv_se = se_iv[1]

    # For each gamma, compute p-value bounds
    p_upper = []
    p_lower = []

    for gamma in gamma_values:
        # Adjust test statistic for potential confounding
        # This is an approximation of Rosenbaum bounds
        t_stat = iv_estimate / iv_se

        # Upper bound: assume confounding works against our estimate
        t_adjusted_upper = t_stat / gamma
        p_up = 2 * (1 - stats.norm.cdf(abs(t_adjusted_upper)))

        # Lower bound: assume confounding works in favor
        t_adjusted_lower = t_stat * gamma
        p_lo = 2 * (1 - stats.norm.cdf(abs(t_adjusted_lower)))

        p_upper.append(p_up)
        p_lower.append(p_lo)

    # Find critical gamma (where upper p-value crosses 0.05)
    critical_gamma = None
    for i, (g, p) in enumerate(zip(gamma_values, p_upper)):
        if p > 0.05:
            if i > 0:
                # Interpolate
                g_prev = gamma_values[i - 1]
                p_prev = p_upper[i - 1]
                critical_gamma = g_prev + (g - g_prev) * (0.05 - p_prev) / (p - p_prev)
            else:
                critical_gamma = g
            break

    return {
        "gamma_values": gamma_values,
        "p_upper": p_upper,
        "p_lower": p_lower,
        "critical_gamma": critical_gamma,
        "iv_estimate": iv_estimate,
        "iv_se": iv_se,
        "interpretation": (
            f"Results remain significant at alpha=0.05 for unobserved confounding "
            f"up to Gamma={critical_gamma:.2f}" if critical_gamma else
            "Results are not robust to unobserved confounding"
        ),
    }
