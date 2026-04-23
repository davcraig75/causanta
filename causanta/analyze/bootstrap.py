"""Bootstrap confidence intervals for causal inference.

Implements BCa (bias-corrected and accelerated) bootstrap for:
- IV/2SLS estimates
- Causal effect estimates
- Segregation statistics
- Effect comparisons (OLS vs IV)

References:
- Efron & Tibshirani (1993) "An Introduction to the Bootstrap"
- DiCiccio & Efron (1996) "Bootstrap Confidence Intervals"
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np
from scipy import stats


@dataclass
class BootstrapResult:
    """Result of bootstrap analysis."""

    point_estimate: float
    se_bootstrap: float
    ci_lower: float
    ci_upper: float
    ci_level: float
    n_bootstrap: int
    method: str
    bootstrap_distribution: np.ndarray
    bias: float
    acceleration: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "point_estimate": self.point_estimate,
            "se_bootstrap": self.se_bootstrap,
            "ci_lower": self.ci_lower,
            "ci_upper": self.ci_upper,
            "ci_level": self.ci_level,
            "n_bootstrap": self.n_bootstrap,
            "method": self.method,
            "bias": self.bias,
            "acceleration": self.acceleration,
        }


def _jackknife_influence(
    data: np.ndarray,
    statistic_fn: Callable[[np.ndarray], float],
) -> np.ndarray:
    """Compute jackknife influence values for acceleration calculation."""
    n = len(data)
    theta_hat = statistic_fn(data)

    # Leave-one-out estimates
    theta_i = np.zeros(n)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        theta_i[i] = statistic_fn(data[mask])

    # Influence values
    theta_bar = np.mean(theta_i)
    influence = (n - 1) * (theta_bar - theta_i)

    return influence


def bootstrap_bca(
    data: np.ndarray,
    statistic_fn: Callable[[np.ndarray], float],
    n_bootstrap: int = 1000,
    ci_level: float = 0.95,
    rng: np.random.Generator | None = None,
) -> BootstrapResult:
    """BCa (bias-corrected and accelerated) bootstrap confidence interval.

    The BCa method corrects for both bias and skewness in the bootstrap
    distribution, providing more accurate coverage than percentile methods.

    Args:
        data: Input data array (n,) or (n, p) for multivariate
        statistic_fn: Function that computes the statistic from data
        n_bootstrap: Number of bootstrap samples
        ci_level: Confidence level (default 0.95)
        rng: Random number generator

    Returns:
        BootstrapResult with CI and diagnostics
    """
    if rng is None:
        rng = np.random.default_rng()

    n = len(data)

    # Original estimate
    theta_hat = statistic_fn(data)

    # Bootstrap distribution
    theta_boot = np.zeros(n_bootstrap)
    for b in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        theta_boot[b] = statistic_fn(data[idx])

    # Bias correction factor (z0)
    prop_below = np.mean(theta_boot < theta_hat)
    # Handle edge cases
    prop_below = np.clip(prop_below, 0.001, 0.999)
    z0 = stats.norm.ppf(prop_below)

    # Acceleration factor (a) from jackknife
    try:
        influence = _jackknife_influence(data, statistic_fn)
        a = np.sum(influence ** 3) / (6 * (np.sum(influence ** 2) ** 1.5))
        a = np.clip(a, -0.5, 0.5)  # Bound for stability
    except Exception:
        a = 0.0  # Fall back to BC (no acceleration)

    # BCa percentiles
    alpha = 1 - ci_level
    z_alpha_lower = stats.norm.ppf(alpha / 2)
    z_alpha_upper = stats.norm.ppf(1 - alpha / 2)

    # Adjusted percentiles
    def adjusted_percentile(z_alpha: float) -> float:
        numerator = z0 + z_alpha
        denominator = 1 - a * numerator
        if abs(denominator) < 1e-10:
            return stats.norm.cdf(z_alpha)
        return stats.norm.cdf(z0 + numerator / denominator)

    p_lower = adjusted_percentile(z_alpha_lower)
    p_upper = adjusted_percentile(z_alpha_upper)

    # Bound percentiles
    p_lower = np.clip(p_lower, 0.001, 0.999)
    p_upper = np.clip(p_upper, 0.001, 0.999)

    ci_lower = np.percentile(theta_boot, 100 * p_lower)
    ci_upper = np.percentile(theta_boot, 100 * p_upper)

    return BootstrapResult(
        point_estimate=theta_hat,
        se_bootstrap=np.std(theta_boot, ddof=1),
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        ci_level=ci_level,
        n_bootstrap=n_bootstrap,
        method="BCa",
        bootstrap_distribution=theta_boot,
        bias=z0,
        acceleration=a,
    )


def bootstrap_percentile(
    data: np.ndarray,
    statistic_fn: Callable[[np.ndarray], float],
    n_bootstrap: int = 1000,
    ci_level: float = 0.95,
    rng: np.random.Generator | None = None,
) -> BootstrapResult:
    """Simple percentile bootstrap (faster, less accurate than BCa)."""
    if rng is None:
        rng = np.random.default_rng()

    n = len(data)
    theta_hat = statistic_fn(data)

    # Bootstrap distribution
    theta_boot = np.zeros(n_bootstrap)
    for b in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        theta_boot[b] = statistic_fn(data[idx])

    alpha = 1 - ci_level
    ci_lower = np.percentile(theta_boot, 100 * alpha / 2)
    ci_upper = np.percentile(theta_boot, 100 * (1 - alpha / 2))

    return BootstrapResult(
        point_estimate=theta_hat,
        se_bootstrap=np.std(theta_boot, ddof=1),
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        ci_level=ci_level,
        n_bootstrap=n_bootstrap,
        method="percentile",
        bootstrap_distribution=theta_boot,
        bias=0.0,
        acceleration=0.0,
    )


def bootstrap_iv_estimate(
    Z: np.ndarray,
    D: np.ndarray,
    Y: np.ndarray,
    covariates: np.ndarray | None = None,
    n_bootstrap: int = 1000,
    ci_level: float = 0.95,
    rng: np.random.Generator | None = None,
) -> BootstrapResult:
    """Bootstrap confidence interval for 2SLS IV estimate.

    Args:
        Z: Instrument (ecDNA count)
        D: Treatment (EGFR expression)
        Y: Outcome
        covariates: Additional control variables
        n_bootstrap: Number of bootstrap samples
        ci_level: Confidence level
        rng: Random number generator

    Returns:
        BootstrapResult for the IV coefficient
    """
    if rng is None:
        rng = np.random.default_rng()

    n = len(Z)

    def iv_statistic(idx: np.ndarray) -> float:
        """Compute 2SLS coefficient for given indices."""
        Z_b = Z[idx]
        D_b = D[idx]
        Y_b = Y[idx]

        if covariates is not None:
            cov_b = covariates[idx]
            X_first = np.column_stack([np.ones(len(idx)), Z_b, cov_b])
        else:
            X_first = np.column_stack([np.ones(len(idx)), Z_b])

        # First stage
        try:
            beta_first = np.linalg.lstsq(X_first, D_b, rcond=None)[0]
            D_hat = X_first @ beta_first
        except Exception:
            return np.nan

        # Second stage
        if covariates is not None:
            X_second = np.column_stack([np.ones(len(idx)), D_hat, cov_b])
        else:
            X_second = np.column_stack([np.ones(len(idx)), D_hat])

        try:
            beta_second = np.linalg.lstsq(X_second, Y_b, rcond=None)[0]
            return beta_second[1]  # Coefficient on D_hat
        except Exception:
            return np.nan

    # Original estimate
    full_idx = np.arange(n)
    theta_hat = iv_statistic(full_idx)

    # Bootstrap
    theta_boot = np.zeros(n_bootstrap)
    for b in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        theta_boot[b] = iv_statistic(idx)

    # Remove NaN values
    theta_boot = theta_boot[~np.isnan(theta_boot)]

    if len(theta_boot) < n_bootstrap * 0.5:
        # Too many failures
        return BootstrapResult(
            point_estimate=theta_hat,
            se_bootstrap=np.nan,
            ci_lower=np.nan,
            ci_upper=np.nan,
            ci_level=ci_level,
            n_bootstrap=len(theta_boot),
            method="BCa (failed)",
            bootstrap_distribution=theta_boot,
            bias=np.nan,
            acceleration=np.nan,
        )

    # BCa calculations
    prop_below = np.mean(theta_boot < theta_hat)
    prop_below = np.clip(prop_below, 0.001, 0.999)
    z0 = stats.norm.ppf(prop_below)

    # Simplified acceleration (skip jackknife for IV which is expensive)
    a = 0.0

    alpha = 1 - ci_level
    z_alpha_lower = stats.norm.ppf(alpha / 2)
    z_alpha_upper = stats.norm.ppf(1 - alpha / 2)

    p_lower = stats.norm.cdf(2 * z0 + z_alpha_lower)
    p_upper = stats.norm.cdf(2 * z0 + z_alpha_upper)

    p_lower = np.clip(p_lower, 0.001, 0.999)
    p_upper = np.clip(p_upper, 0.001, 0.999)

    ci_lower = np.percentile(theta_boot, 100 * p_lower)
    ci_upper = np.percentile(theta_boot, 100 * p_upper)

    return BootstrapResult(
        point_estimate=theta_hat,
        se_bootstrap=np.std(theta_boot, ddof=1),
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        ci_level=ci_level,
        n_bootstrap=len(theta_boot),
        method="BCa",
        bootstrap_distribution=theta_boot,
        bias=z0,
        acceleration=a,
    )


def bootstrap_segregation_test(
    daughter_fractions: np.ndarray,
    null_value: float = 0.5,
    n_bootstrap: int = 1000,
    ci_level: float = 0.95,
    rng: np.random.Generator | None = None,
) -> dict[str, Any]:
    """Bootstrap test for ecDNA segregation following Binomial(N, 0.5).

    Tests whether the mean daughter fraction equals 0.5 as expected
    under random segregation.

    Args:
        daughter_fractions: Array of daughter ecDNA fractions
        null_value: Expected mean under null (0.5)
        n_bootstrap: Number of bootstrap samples
        ci_level: Confidence level
        rng: Random number generator

    Returns:
        Dictionary with test results and CI
    """
    if rng is None:
        rng = np.random.default_rng()

    n = len(daughter_fractions)
    observed_mean = np.mean(daughter_fractions)

    # Bootstrap CI for mean
    def mean_fn(data: np.ndarray) -> float:
        return np.mean(data)

    result = bootstrap_bca(
        daughter_fractions, mean_fn, n_bootstrap, ci_level, rng
    )

    # P-value for two-sided test of H0: mean = 0.5
    # Use bootstrap distribution centered at null
    centered_boot = result.bootstrap_distribution - observed_mean + null_value
    p_value = np.mean(np.abs(centered_boot - null_value) >= np.abs(observed_mean - null_value))

    # Also compute variance test
    observed_var = np.var(daughter_fractions, ddof=1)

    return {
        "mean": observed_mean,
        "mean_ci_lower": result.ci_lower,
        "mean_ci_upper": result.ci_upper,
        "expected_mean": null_value,
        "p_value_mean": p_value,
        "variance": observed_var,
        "n_observations": n,
        "ci_level": ci_level,
        "contains_null": result.ci_lower <= null_value <= result.ci_upper,
        "interpretation": (
            "Segregation consistent with Binomial(N, 0.5)"
            if result.ci_lower <= null_value <= result.ci_upper
            else "Segregation deviates from expected"
        ),
    }


def bootstrap_effect_comparison(
    Z: np.ndarray,
    D: np.ndarray,
    Y: np.ndarray,
    covariates: np.ndarray | None = None,
    n_bootstrap: int = 1000,
    ci_level: float = 0.95,
    rng: np.random.Generator | None = None,
) -> dict[str, Any]:
    """Bootstrap comparison of OLS vs IV estimates.

    Computes bootstrap CIs for both estimators and their difference,
    allowing formal testing of whether confounding bias is significant.

    Returns:
        Dictionary with OLS, IV, and difference results
    """
    if rng is None:
        rng = np.random.default_rng()

    n = len(Z)

    def compute_both(idx: np.ndarray) -> tuple[float, float]:
        """Compute both OLS and IV for given indices."""
        Z_b = Z[idx]
        D_b = D[idx]
        Y_b = Y[idx]

        if covariates is not None:
            cov_b = covariates[idx]
            X_ols = np.column_stack([np.ones(len(idx)), D_b, cov_b])
            X_first = np.column_stack([np.ones(len(idx)), Z_b, cov_b])
        else:
            X_ols = np.column_stack([np.ones(len(idx)), D_b])
            X_first = np.column_stack([np.ones(len(idx)), Z_b])

        # OLS
        try:
            beta_ols = np.linalg.lstsq(X_ols, Y_b, rcond=None)[0]
            ols_coef = beta_ols[1]
        except Exception:
            ols_coef = np.nan

        # IV (2SLS)
        try:
            beta_first = np.linalg.lstsq(X_first, D_b, rcond=None)[0]
            D_hat = X_first @ beta_first

            if covariates is not None:
                X_second = np.column_stack([np.ones(len(idx)), D_hat, cov_b])
            else:
                X_second = np.column_stack([np.ones(len(idx)), D_hat])

            beta_second = np.linalg.lstsq(X_second, Y_b, rcond=None)[0]
            iv_coef = beta_second[1]
        except Exception:
            iv_coef = np.nan

        return ols_coef, iv_coef

    # Original estimates
    full_idx = np.arange(n)
    ols_hat, iv_hat = compute_both(full_idx)

    # Bootstrap
    ols_boot = np.zeros(n_bootstrap)
    iv_boot = np.zeros(n_bootstrap)

    for b in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        ols_boot[b], iv_boot[b] = compute_both(idx)

    # Remove NaN values
    valid = ~(np.isnan(ols_boot) | np.isnan(iv_boot))
    ols_boot = ols_boot[valid]
    iv_boot = iv_boot[valid]

    if len(ols_boot) < n_bootstrap * 0.5:
        return {"error": "Too many bootstrap failures"}

    # Difference distribution (measures confounding bias)
    diff_boot = ols_boot - iv_boot
    diff_hat = ols_hat - iv_hat

    alpha = 1 - ci_level

    return {
        "ols_estimate": ols_hat,
        "ols_se": np.std(ols_boot, ddof=1),
        "ols_ci_lower": np.percentile(ols_boot, 100 * alpha / 2),
        "ols_ci_upper": np.percentile(ols_boot, 100 * (1 - alpha / 2)),
        "iv_estimate": iv_hat,
        "iv_se": np.std(iv_boot, ddof=1),
        "iv_ci_lower": np.percentile(iv_boot, 100 * alpha / 2),
        "iv_ci_upper": np.percentile(iv_boot, 100 * (1 - alpha / 2)),
        "difference": diff_hat,
        "difference_se": np.std(diff_boot, ddof=1),
        "difference_ci_lower": np.percentile(diff_boot, 100 * alpha / 2),
        "difference_ci_upper": np.percentile(diff_boot, 100 * (1 - alpha / 2)),
        "bias_significant": (
            np.percentile(diff_boot, 100 * alpha / 2) > 0 or
            np.percentile(diff_boot, 100 * (1 - alpha / 2)) < 0
        ),
        "n_bootstrap": len(ols_boot),
        "ci_level": ci_level,
    }


def bootstrap_causal_effects(
    cells_data: list[dict],
    outcomes: list[str] | None = None,
    n_bootstrap: int = 1000,
    ci_level: float = 0.95,
    rng: np.random.Generator | None = None,
) -> dict[str, BootstrapResult]:
    """Bootstrap CIs for all causal effect estimates.

    Args:
        cells_data: List of cell dictionaries from final timestep
        outcomes: List of outcome variable names
        n_bootstrap: Number of bootstrap samples
        ci_level: Confidence level
        rng: Random number generator

    Returns:
        Dictionary mapping outcome names to BootstrapResult
    """
    if rng is None:
        rng = np.random.default_rng()

    if outcomes is None:
        outcomes = ["VEGF_secretion", "migration_rate"]

    # Filter to tumor cells
    tumor_cells = [c for c in cells_data if c.get("cell_type") == 6]
    if len(tumor_cells) < 30:
        return {"error": "Insufficient tumor cells"}

    # Extract arrays
    Z = np.array([c["ecDNA_count"] for c in tumor_cells], dtype=float)
    D = np.array([c.get("egfr_expression", c["ecDNA_count"]) for c in tumor_cells], dtype=float)
    O2 = np.array([c["O2_local"] for c in tumor_cells], dtype=float)
    glucose = np.array([c["glucose_local"] for c in tumor_cells], dtype=float)
    covariates = np.column_stack([O2, glucose])

    results = {}

    for outcome_name in outcomes:
        if outcome_name not in tumor_cells[0]:
            continue

        Y = np.array([c[outcome_name] for c in tumor_cells], dtype=float)

        result = bootstrap_iv_estimate(
            Z, D, Y, covariates, n_bootstrap, ci_level, rng
        )
        results[outcome_name] = result

    return results
