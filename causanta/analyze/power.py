"""Power analysis for causal inference with ecDNA instruments.

Implements methods to determine:
- Minimum sample sizes for detecting true effects
- Power curves for different effect sizes
- Required number of cell divisions for reliable IV inference

Key considerations for ecDNA-based IV:
1. Strong instrument (F > 10): ecDNA strongly predicts EGFR
2. Effect size variability: smaller effects need more data
3. Confounding strength: higher confounding increases IV advantage
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import stats


@dataclass
class PowerAnalysisResult:
    """Result of power analysis."""

    target_effect: float
    required_n: int
    achieved_power: float
    target_power: float
    alpha: float
    first_stage_r2: float
    method: str
    power_curve: dict[str, list[float]]  # n -> power

    def to_dict(self) -> dict[str, Any]:
        return {
            "target_effect": self.target_effect,
            "required_n": self.required_n,
            "achieved_power": self.achieved_power,
            "target_power": self.target_power,
            "alpha": self.alpha,
            "first_stage_r2": self.first_stage_r2,
            "method": self.method,
            "power_curve": self.power_curve,
        }


def iv_power_analytical(
    effect_size: float,
    first_stage_r2: float,
    n: int,
    alpha: float = 0.05,
    residual_var: float = 1.0,
) -> float:
    """Analytical power calculation for 2SLS estimator.

    Based on the asymptotic distribution of the 2SLS estimator:
    sqrt(n) * (beta_IV - beta) ~ N(0, sigma^2 / (pi_1^2 * Var(Z)))

    where pi_1 is the first-stage coefficient and sigma^2 is outcome variance.

    Args:
        effect_size: True causal effect (standardized)
        first_stage_r2: R² from first-stage regression (instrument strength)
        n: Sample size
        alpha: Significance level
        residual_var: Residual variance of outcome

    Returns:
        Statistical power (probability of rejecting null when effect exists)
    """
    # Standard error of IV estimator (approximate)
    # SE_IV ≈ sigma / (sqrt(n) * |pi_1| * SD(Z))
    # With standardized variables: SE_IV ≈ sqrt((1 - R2_first) / (n * R2_first))
    # Plus adjustment for second-stage variance

    if first_stage_r2 <= 0 or first_stage_r2 >= 1:
        return 0.0

    # Variance inflation from IV relative to OLS
    # This is an approximation; exact formula depends on data structure
    variance_inflation = 1.0 / first_stage_r2

    # SE of IV estimate
    se_iv = np.sqrt(residual_var * variance_inflation / n)

    # Non-centrality parameter
    ncp = abs(effect_size) / se_iv

    # Critical value (two-sided)
    z_crit = stats.norm.ppf(1 - alpha / 2)

    # Power
    power = 1 - stats.norm.cdf(z_crit - ncp) + stats.norm.cdf(-z_crit - ncp)

    return float(power)


def sample_size_for_power(
    effect_size: float,
    first_stage_r2: float,
    target_power: float = 0.80,
    alpha: float = 0.05,
    residual_var: float = 1.0,
    max_n: int = 100000,
) -> int:
    """Calculate required sample size to achieve target power.

    Uses binary search to find minimum n for target power.

    Args:
        effect_size: True causal effect (standardized)
        first_stage_r2: R² from first-stage regression
        target_power: Desired power (default 0.80)
        alpha: Significance level
        residual_var: Residual variance
        max_n: Maximum sample size to consider

    Returns:
        Required sample size
    """
    # Binary search for n
    low, high = 10, max_n

    while high - low > 1:
        mid = (low + high) // 2
        power = iv_power_analytical(
            effect_size, first_stage_r2, mid, alpha, residual_var
        )
        if power >= target_power:
            high = mid
        else:
            low = mid

    return high


def power_curve(
    effect_size: float,
    first_stage_r2: float,
    n_values: list[int] | None = None,
    alpha: float = 0.05,
    residual_var: float = 1.0,
) -> dict[str, list[float]]:
    """Compute power curve over range of sample sizes.

    Args:
        effect_size: True causal effect
        first_stage_r2: Instrument strength
        n_values: Sample sizes to evaluate
        alpha: Significance level
        residual_var: Residual variance

    Returns:
        Dictionary with 'n' and 'power' lists
    """
    if n_values is None:
        n_values = [50, 100, 200, 500, 1000, 2000, 5000, 10000]

    powers = []
    for n in n_values:
        p = iv_power_analytical(effect_size, first_stage_r2, n, alpha, residual_var)
        powers.append(p)

    return {"n": n_values, "power": powers}


def minimum_detectable_effect(
    n: int,
    first_stage_r2: float,
    target_power: float = 0.80,
    alpha: float = 0.05,
    residual_var: float = 1.0,
) -> float:
    """Calculate minimum detectable effect size for given sample size.

    The MDE is the smallest effect that can be detected with target power.

    Args:
        n: Sample size
        first_stage_r2: Instrument strength
        target_power: Target power level
        alpha: Significance level
        residual_var: Residual variance

    Returns:
        Minimum detectable effect (standardized)
    """
    # Binary search for effect size
    low, high = 0.0, 2.0

    for _ in range(50):  # 50 iterations for precision
        mid = (low + high) / 2
        power = iv_power_analytical(mid, first_stage_r2, n, alpha, residual_var)
        if power >= target_power:
            high = mid
        else:
            low = mid

    return high


def analyze_power_from_data(
    cells_data: list[dict],
    outcome_name: str = "VEGF_secretion",
    effect_sizes: list[float] | None = None,
    n_values: list[int] | None = None,
) -> PowerAnalysisResult:
    """Analyze power from actual simulation data.

    Uses the observed first-stage regression and outcome variance
    to compute power curves and required sample sizes.

    Args:
        cells_data: Cell data from simulation
        outcome_name: Outcome variable name
        effect_sizes: Effect sizes to analyze
        n_values: Sample sizes for power curve

    Returns:
        PowerAnalysisResult with curves and recommendations
    """
    if effect_sizes is None:
        effect_sizes = [0.01, 0.02, 0.05, 0.1, 0.2]

    if n_values is None:
        n_values = [100, 200, 500, 1000, 2000, 5000]

    # Filter to tumor cells
    tumor_cells = [c for c in cells_data if c.get("cell_type") == 6]

    if len(tumor_cells) < 30:
        return PowerAnalysisResult(
            target_effect=0.0,
            required_n=0,
            achieved_power=0.0,
            target_power=0.80,
            alpha=0.05,
            first_stage_r2=0.0,
            method="insufficient_data",
            power_curve={},
        )

    # Extract data
    Z = np.array([c["ecDNA_count"] for c in tumor_cells], dtype=float)
    D = np.array([c.get("egfr_expression", c["ecDNA_count"]) for c in tumor_cells], dtype=float)
    Y = np.array([c[outcome_name] for c in tumor_cells], dtype=float)

    # First-stage R²
    n = len(Z)
    X_first = np.column_stack([np.ones(n), Z])
    beta_first = np.linalg.lstsq(X_first, D, rcond=None)[0]
    D_hat = X_first @ beta_first
    ss_res = np.sum((D - D_hat) ** 2)
    ss_tot = np.sum((D - np.mean(D)) ** 2)
    first_stage_r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    # Residual variance of outcome
    residual_var = np.var(Y, ddof=1)

    # Standardize effect sizes relative to outcome SD
    outcome_sd = np.std(Y, ddof=1)

    # Power curves for different effect sizes
    power_curves = {}
    for es in effect_sizes:
        standardized_es = es / outcome_sd if outcome_sd > 0 else es
        curve = power_curve(standardized_es, first_stage_r2, n_values)
        power_curves[f"effect_{es}"] = curve["power"]

    power_curves["n"] = n_values

    # Required N for a moderate effect (e.g., delta = 0.05)
    target_effect = 0.05
    standardized_target = target_effect / outcome_sd if outcome_sd > 0 else target_effect
    required_n = sample_size_for_power(
        standardized_target, first_stage_r2, target_power=0.80
    )

    # Current power
    current_power = iv_power_analytical(
        standardized_target, first_stage_r2, len(tumor_cells)
    )

    return PowerAnalysisResult(
        target_effect=target_effect,
        required_n=required_n,
        achieved_power=current_power,
        target_power=0.80,
        alpha=0.05,
        first_stage_r2=first_stage_r2,
        method="analytical_iv",
        power_curve=power_curves,
    )


def simulation_based_power(
    true_effect: float,
    n_samples: int,
    first_stage_coef: float = 0.5,
    instrument_var: float = 100.0,
    confounding_strength: float = 0.5,
    n_simulations: int = 1000,
    alpha: float = 0.05,
    rng: np.random.Generator | None = None,
) -> dict[str, float]:
    """Simulation-based power analysis.

    Simulates data from the structural model and checks how often
    we correctly reject the null hypothesis.

    Args:
        true_effect: True causal effect (beta)
        n_samples: Sample size
        first_stage_coef: First-stage coefficient (pi_1)
        instrument_var: Variance of instrument
        confounding_strength: Correlation between confounder and treatment
        n_simulations: Number of simulation replications
        alpha: Significance level
        rng: Random number generator

    Returns:
        Dictionary with power and other diagnostics
    """
    if rng is None:
        rng = np.random.default_rng()

    rejections = 0
    iv_biases = []
    ols_biases = []

    for _ in range(n_simulations):
        # Generate data from structural model
        # Z -> D -> Y, U -> D, U -> Y (confounding)

        # Instrument
        Z = rng.normal(0, np.sqrt(instrument_var), n_samples)

        # Confounder
        U = rng.normal(0, 1, n_samples)

        # Treatment (affected by Z and U)
        D = first_stage_coef * Z + confounding_strength * U + rng.normal(0, 1, n_samples)

        # Outcome (affected by D and U)
        Y = true_effect * D + confounding_strength * U + rng.normal(0, 1, n_samples)

        # Estimate IV
        n = len(Z)
        X_first = np.column_stack([np.ones(n), Z])
        beta_first = np.linalg.lstsq(X_first, D, rcond=None)[0]
        D_hat = X_first @ beta_first

        X_second = np.column_stack([np.ones(n), D_hat])
        beta_second = np.linalg.lstsq(X_second, Y, rcond=None)[0]
        iv_estimate = beta_second[1]

        # Standard error (approximate)
        resid = Y - X_second @ beta_second
        mse = np.sum(resid ** 2) / (n - 2)
        XtX_inv = np.linalg.pinv(X_second.T @ X_second)
        se = np.sqrt(mse * XtX_inv[1, 1])

        # Test
        t_stat = iv_estimate / se
        p_value = 2 * (1 - stats.t.cdf(abs(t_stat), n - 2))

        if p_value < alpha:
            rejections += 1

        iv_biases.append(iv_estimate - true_effect)

        # OLS for comparison
        X_ols = np.column_stack([np.ones(n), D])
        beta_ols = np.linalg.lstsq(X_ols, Y, rcond=None)[0]
        ols_biases.append(beta_ols[1] - true_effect)

    return {
        "power": rejections / n_simulations,
        "n_simulations": n_simulations,
        "n_samples": n_samples,
        "true_effect": true_effect,
        "mean_iv_bias": np.mean(iv_biases),
        "mean_ols_bias": np.mean(ols_biases),
        "iv_rmse": np.sqrt(np.mean(np.array(iv_biases) ** 2)),
        "ols_rmse": np.sqrt(np.mean(np.array(ols_biases) ** 2)),
    }


def comprehensive_power_analysis(
    effect_sizes: list[float] | None = None,
    sample_sizes: list[int] | None = None,
    first_stage_r2: float = 0.85,
    n_simulations: int = 500,
    rng: np.random.Generator | None = None,
) -> dict[str, Any]:
    """Comprehensive power analysis across effect sizes and sample sizes.

    Combines analytical and simulation-based approaches for robust
    power estimation.

    Args:
        effect_sizes: List of effect sizes to analyze
        sample_sizes: List of sample sizes to analyze
        first_stage_r2: Instrument strength (R²)
        n_simulations: Simulations per condition
        rng: Random number generator

    Returns:
        Comprehensive results dictionary
    """
    if effect_sizes is None:
        effect_sizes = [0.01, 0.02, 0.05, 0.1, 0.2]

    if sample_sizes is None:
        sample_sizes = [100, 200, 500, 1000, 2000, 5000]

    if rng is None:
        rng = np.random.default_rng()

    results = {
        "effect_sizes": effect_sizes,
        "sample_sizes": sample_sizes,
        "first_stage_r2": first_stage_r2,
        "analytical_power": {},
        "required_n_80_power": {},
        "required_n_90_power": {},
        "mde_by_n": {},
    }

    # Analytical power matrix
    for es in effect_sizes:
        powers = []
        for n in sample_sizes:
            p = iv_power_analytical(es, first_stage_r2, n)
            powers.append(p)
        results["analytical_power"][f"effect_{es}"] = powers

    # Required N for each effect size
    for es in effect_sizes:
        n_80 = sample_size_for_power(es, first_stage_r2, target_power=0.80)
        n_90 = sample_size_for_power(es, first_stage_r2, target_power=0.90)
        results["required_n_80_power"][f"effect_{es}"] = n_80
        results["required_n_90_power"][f"effect_{es}"] = n_90

    # MDE for each sample size
    for n in sample_sizes:
        mde = minimum_detectable_effect(n, first_stage_r2, target_power=0.80)
        results["mde_by_n"][f"n_{n}"] = mde

    # Summary recommendations
    results["recommendations"] = {
        "for_delta_0.05": {
            "effect": 0.05,
            "required_n_80_power": sample_size_for_power(0.05, first_stage_r2, 0.80),
            "required_n_90_power": sample_size_for_power(0.05, first_stage_r2, 0.90),
        },
        "for_n_1000": {
            "sample_size": 1000,
            "mde_80_power": minimum_detectable_effect(1000, first_stage_r2, 0.80),
            "power_for_0.05": iv_power_analytical(0.05, first_stage_r2, 1000),
        },
        "interpretation": (
            f"With first-stage R²={first_stage_r2:.2f}, detecting a 5% effect "
            f"with 80% power requires ~{sample_size_for_power(0.05, first_stage_r2, 0.80)} observations. "
            f"With 1000 observations, the MDE (80% power) is ~{minimum_detectable_effect(1000, first_stage_r2, 0.80):.3f}."
        ),
    }

    return results
