"""Structural Equation Modeling (SEM) for causal inference.

Implements path analysis and structural equation models for
estimating causal effects in systems with multiple mediators.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy import stats

from .loader import SimulationData, cells_to_array


@dataclass
class PathCoefficient:
    """A single path coefficient in the SEM."""

    from_var: str
    to_var: str
    coefficient: float
    standard_error: float
    z_value: float
    p_value: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "from": self.from_var,
            "to": self.to_var,
            "coefficient": self.coefficient,
            "se": self.standard_error,
            "z": self.z_value,
            "p": self.p_value,
        }


@dataclass
class SEMResult:
    """Result of structural equation model estimation."""

    model_name: str
    paths: list[PathCoefficient]
    direct_effects: dict[str, float]
    indirect_effects: dict[str, float]
    total_effects: dict[str, float]
    fit_indices: dict[str, float]
    n_observations: int
    residual_variances: dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "model": self.model_name,
            "paths": [p.to_dict() for p in self.paths],
            "direct_effects": self.direct_effects,
            "indirect_effects": self.indirect_effects,
            "total_effects": self.total_effects,
            "fit_indices": self.fit_indices,
            "n_observations": self.n_observations,
        }


def estimate_path_model(
    X: np.ndarray,
    outcome_idx: int,
    treatment_idx: int,
    mediator_indices: list[int],
    variable_names: list[str],
) -> tuple[list[PathCoefficient], dict[str, float]]:
    """Estimate a simple path model via sequential regression.

    Treatment -> Mediators -> Outcome

    Args:
        X: Data matrix (n x p)
        outcome_idx: Index of outcome variable
        treatment_idx: Index of treatment variable
        mediator_indices: Indices of mediator variables
        variable_names: Names for each variable

    Returns:
        (path_coefficients, residual_variances)
    """
    n, p = X.shape
    paths = []
    residual_vars = {}

    treatment_name = variable_names[treatment_idx]
    outcome_name = variable_names[outcome_idx]

    # Path from treatment to each mediator
    for m_idx in mediator_indices:
        m_name = variable_names[m_idx]

        # Regression: mediator ~ treatment
        X_t = np.column_stack([np.ones(n), X[:, treatment_idx]])
        y_m = X[:, m_idx]

        beta = np.linalg.lstsq(X_t, y_m, rcond=None)[0]
        y_hat = X_t @ beta
        resid = y_m - y_hat

        mse = np.sum(resid ** 2) / (n - 2)
        var_beta = mse * np.linalg.pinv(X_t.T @ X_t)
        se = np.sqrt(var_beta[1, 1])

        z = beta[1] / se if se > 0 else 0
        pval = 2 * (1 - stats.norm.cdf(abs(z)))

        paths.append(PathCoefficient(
            from_var=treatment_name,
            to_var=m_name,
            coefficient=float(beta[1]),
            standard_error=float(se),
            z_value=float(z),
            p_value=float(pval),
        ))

        residual_vars[m_name] = float(np.var(resid, ddof=1))

    # Path from treatment and mediators to outcome
    X_full = np.column_stack([
        np.ones(n),
        X[:, treatment_idx],
        X[:, mediator_indices],
    ])
    y_out = X[:, outcome_idx]

    beta_full = np.linalg.lstsq(X_full, y_out, rcond=None)[0]
    y_hat = X_full @ beta_full
    resid = y_out - y_hat

    mse = np.sum(resid ** 2) / (n - X_full.shape[1])
    var_beta = mse * np.linalg.pinv(X_full.T @ X_full)

    # Direct effect: treatment -> outcome (controlling for mediators)
    se_direct = np.sqrt(var_beta[1, 1])
    z_direct = beta_full[1] / se_direct if se_direct > 0 else 0
    pval_direct = 2 * (1 - stats.norm.cdf(abs(z_direct)))

    paths.append(PathCoefficient(
        from_var=treatment_name,
        to_var=outcome_name,
        coefficient=float(beta_full[1]),
        standard_error=float(se_direct),
        z_value=float(z_direct),
        p_value=float(pval_direct),
    ))

    # Mediator -> outcome paths
    for i, m_idx in enumerate(mediator_indices):
        m_name = variable_names[m_idx]
        idx = 2 + i  # Position in beta_full

        se_m = np.sqrt(var_beta[idx, idx])
        z_m = beta_full[idx] / se_m if se_m > 0 else 0
        pval_m = 2 * (1 - stats.norm.cdf(abs(z_m)))

        paths.append(PathCoefficient(
            from_var=m_name,
            to_var=outcome_name,
            coefficient=float(beta_full[idx]),
            standard_error=float(se_m),
            z_value=float(z_m),
            p_value=float(pval_m),
        ))

    residual_vars[outcome_name] = float(np.var(resid, ddof=1))

    return paths, residual_vars


def compute_effects(
    paths: list[PathCoefficient],
    treatment_name: str,
    outcome_name: str,
    mediator_names: list[str],
) -> tuple[dict[str, float], dict[str, float], dict[str, float]]:
    """Compute direct, indirect, and total effects from path coefficients.

    Args:
        paths: List of path coefficients
        treatment_name: Name of treatment variable
        outcome_name: Name of outcome variable
        mediator_names: Names of mediator variables

    Returns:
        (direct_effects, indirect_effects, total_effects)
    """
    # Build coefficient lookup
    coef_lookup = {(p.from_var, p.to_var): p.coefficient for p in paths}

    # Direct effect: treatment -> outcome
    direct = coef_lookup.get((treatment_name, outcome_name), 0.0)

    # Indirect effects via each mediator
    indirect_by_mediator = {}
    total_indirect = 0.0

    for m in mediator_names:
        # treatment -> mediator
        a = coef_lookup.get((treatment_name, m), 0.0)
        # mediator -> outcome
        b = coef_lookup.get((m, outcome_name), 0.0)
        indirect = a * b
        indirect_by_mediator[f"via_{m}"] = indirect
        total_indirect += indirect

    direct_effects = {f"{treatment_name}_to_{outcome_name}": direct}
    indirect_effects = {"total_indirect": total_indirect, **indirect_by_mediator}
    total_effects = {f"{treatment_name}_total": direct + total_indirect}

    return direct_effects, indirect_effects, total_effects


def mediation_analysis(
    data: SimulationData,
    treatment: str = "ecDNA_count",
    mediator: str = "VEGF_secretion",
    outcome: str = "migration_rate",
    covariates: list[str] | None = None,
) -> SEMResult:
    """Perform mediation analysis for a single mediator.

    Treatment -> Mediator -> Outcome

    Args:
        data: Loaded simulation data
        treatment: Treatment variable name
        mediator: Mediator variable name
        outcome: Outcome variable name
        covariates: Additional covariates to control for

    Returns:
        SEMResult with decomposed effects.
    """
    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 20:
        return SEMResult(
            model_name="mediation",
            paths=[],
            direct_effects={},
            indirect_effects={},
            total_effects={},
            fit_indices={"error": "insufficient_data"},
            n_observations=0,
        )

    # Build variable list
    variables = [treatment, mediator, outcome]
    if covariates:
        variables.extend(covariates)

    X = cells_to_array(tumor_cells, variables)
    n = len(tumor_cells)

    # Estimate paths
    paths, residual_vars = estimate_path_model(
        X,
        outcome_idx=2,  # outcome
        treatment_idx=0,  # treatment
        mediator_indices=[1],  # mediator
        variable_names=variables,
    )

    # Compute effects
    direct, indirect, total = compute_effects(
        paths, treatment, outcome, [mediator]
    )

    # Simple fit indices
    X_model = np.column_stack([np.ones(n), X[:, 0], X[:, 1]])  # treatment + mediator
    y = X[:, 2]
    beta = np.linalg.lstsq(X_model, y, rcond=None)[0]
    y_hat = X_model @ beta
    ss_res = np.sum((y - y_hat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    fit_indices = {
        "r_squared": float(r_sq),
        "n_paths": len(paths),
    }

    return SEMResult(
        model_name="mediation",
        paths=paths,
        direct_effects=direct,
        indirect_effects=indirect,
        total_effects=total,
        fit_indices=fit_indices,
        n_observations=n,
        residual_variances=residual_vars,
    )


def full_sem(
    data: SimulationData,
    model_spec: dict[str, list[str]] | None = None,
) -> SEMResult:
    """Estimate full SEM based on specification.

    Args:
        data: Loaded simulation data
        model_spec: Dictionary mapping each endogenous variable to its
                    predictors. If None, uses default CAUSANTA model:
                    {
                        "VEGF_secretion": ["ecDNA_count", "O2_local"],
                        "migration_rate": ["ecDNA_count", "VEGF_secretion"],
                    }

    Returns:
        SEMResult with all estimated paths.
    """
    if model_spec is None:
        model_spec = {
            "VEGF_secretion": ["ecDNA_count", "O2_local"],
            "migration_rate": ["ecDNA_count", "VEGF_secretion", "O2_local"],
        }

    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 20:
        return SEMResult(
            model_name="full_sem",
            paths=[],
            direct_effects={},
            indirect_effects={},
            total_effects={},
            fit_indices={"error": "insufficient_data"},
            n_observations=0,
        )

    # Collect all variables
    all_vars = set()
    for outcome, predictors in model_spec.items():
        all_vars.add(outcome)
        all_vars.update(predictors)

    all_vars = sorted(all_vars)
    var_idx = {v: i for i, v in enumerate(all_vars)}

    X = cells_to_array(tumor_cells, all_vars)
    n = len(tumor_cells)

    paths = []
    residual_vars = {}

    # Estimate each equation
    for outcome, predictors in model_spec.items():
        y = X[:, var_idx[outcome]]
        X_pred = np.column_stack([np.ones(n)] + [X[:, var_idx[p]] for p in predictors])

        beta = np.linalg.lstsq(X_pred, y, rcond=None)[0]
        y_hat = X_pred @ beta
        resid = y - y_hat

        mse = np.sum(resid ** 2) / (n - len(predictors) - 1)
        var_beta = mse * np.linalg.pinv(X_pred.T @ X_pred)

        for i, pred in enumerate(predictors):
            coef = beta[i + 1]
            se = np.sqrt(var_beta[i + 1, i + 1])
            z = coef / se if se > 0 else 0
            pval = 2 * (1 - stats.norm.cdf(abs(z)))

            paths.append(PathCoefficient(
                from_var=pred,
                to_var=outcome,
                coefficient=float(coef),
                standard_error=float(se),
                z_value=float(z),
                p_value=float(pval),
            ))

        residual_vars[outcome] = float(np.var(resid, ddof=1))

    # Compute effects for ecDNA -> migration_rate
    direct, indirect, total = compute_effects(
        paths,
        "ecDNA_count",
        "migration_rate",
        ["VEGF_secretion"],
    )

    fit_indices = {
        "n_paths": len(paths),
        "n_equations": len(model_spec),
    }

    return SEMResult(
        model_name="full_sem",
        paths=paths,
        direct_effects=direct,
        indirect_effects=indirect,
        total_effects=total,
        fit_indices=fit_indices,
        n_observations=n,
        residual_variances=residual_vars,
    )


def sobel_test(
    a: float,
    b: float,
    se_a: float,
    se_b: float,
) -> tuple[float, float]:
    """Sobel test for mediation significance.

    Tests whether the indirect effect a*b is significantly different from 0.

    Args:
        a: Path coefficient from treatment to mediator
        b: Path coefficient from mediator to outcome
        se_a: Standard error of a
        se_b: Standard error of b

    Returns:
        (z_statistic, p_value)
    """
    # Sobel standard error
    se_indirect = np.sqrt(a**2 * se_b**2 + b**2 * se_a**2)

    # Z-statistic
    z = (a * b) / se_indirect if se_indirect > 0 else 0

    # Two-tailed p-value
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))

    return float(z), float(p_value)
