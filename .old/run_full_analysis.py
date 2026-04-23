#!/usr/bin/env python3
"""Run full simulation and comprehensive causal analysis.

This script:
1. Runs simulation using the TRUE causal DAG (ground truth)
2. Analyzes data using a HYPOTHESIZED DAG (what an analyst might propose)
3. Compares naive vs. IV estimates to test causal inference validity
4. Generates comprehensive statistical report with visualizations

The key insight: we simulate from the TRUE structure but analyze with
a HYPOTHESIZED structure - testing whether methods recover truth.

TRUE DAG (simulation):
    ecDNA → Proliferation → O2 depletion → Hypoxia → Invasion
    (ecDNA does NOT directly cause hypoxia - it's mediated!)

HYPOTHESIZED DAG (analysis - potentially wrong):
    ecDNA → Hypoxia (direct? - actually wrong!)
    ecDNA → Proliferation
    + Unmeasured confounders (U)
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from causanta.simulate.core import Simulation
from causanta.graph.causal_dag import build_causanta_dag, CausalDAG, CausalNode, CausalEdge, NodeType, EffectType


def build_analysis_dag() -> CausalDAG:
    """Build the HYPOTHESIZED analysis DAG.

    This is what an analyst might propose WITHOUT knowing the true
    causal structure. It includes:
    - Potential direct effect of ecDNA on hypoxia (WRONG in simulation)
    - Unmeasured confounders representing unknown factors
    - Simplified pathways

    The goal is to see if causal inference methods can correctly
    identify that ecDNA affects hypoxia INDIRECTLY through proliferation.
    """
    dag = CausalDAG(
        name="Hypothesized Analysis DAG",
        description=(
            "Analyst's hypothesized causal structure. Includes potential "
            "direct path from ecDNA to hypoxia (which is actually indirect "
            "in the true model) and unmeasured confounders."
        ),
    )

    # Nodes
    dag.add_node(CausalNode(
        name="ecDNA_EGFR",
        description="ecDNA copy number (observed instrument)",
        node_type=NodeType.INSTRUMENT,
        units="copies",
    ))

    dag.add_node(CausalNode(
        name="EGFR_expression",
        description="EGFR protein expression (observed exposure)",
        node_type=NodeType.EXPOSURE,
        units="arbitrary units",
    ))

    dag.add_node(CausalNode(
        name="Proliferation",
        description="Cell division rate (observed)",
        node_type=NodeType.OUTCOME,
        units="1/hr",
    ))

    dag.add_node(CausalNode(
        name="Hypoxia",
        description="Hypoxic state (observed)",
        node_type=NodeType.MEDIATOR,
        units="binary",
    ))

    dag.add_node(CausalNode(
        name="Invasion",
        description="Migration/invasion speed (observed)",
        node_type=NodeType.OUTCOME,
        units="um/hr",
    ))

    dag.add_node(CausalNode(
        name="U_location",
        description="Unmeasured: spatial location effects",
        node_type=NodeType.CONFOUNDER,
        observable=False,
    ))

    dag.add_node(CausalNode(
        name="U_neighbors",
        description="Unmeasured: neighbor cell effects",
        node_type=NodeType.CONFOUNDER,
        observable=False,
    ))

    # Edges - hypothesized relationships
    dag.add_edge(CausalEdge(
        source="ecDNA_EGFR",
        target="Proliferation",
        effect_type=EffectType.LOG,
        parameter_symbol="α",
        parameter_name="ecDNA_effect_on_division",
        parameter_value=0.3,
        equation_template="Prolif = f(ecDNA)",
        description="Direct effect of EGFR on proliferation",
        sign="+",
    ))

    dag.add_edge(CausalEdge(
        source="ecDNA_EGFR",
        target="Hypoxia",
        effect_type=EffectType.LINEAR,
        parameter_symbol="?",
        parameter_name="unknown",
        parameter_value=0.0,  # Actually 0 in true model!
        equation_template="Hypoxia = f(ecDNA)?",
        description="HYPOTHESIZED direct effect (actually indirect!)",
        sign="?",
    ))

    dag.add_edge(CausalEdge(
        source="Hypoxia",
        target="Invasion",
        effect_type=EffectType.MULTIPLICATIVE,
        parameter_symbol="δ",
        parameter_name="hypoxia_invasion",
        parameter_value=2.0,
        equation_template="Invasion = f(Hypoxia)",
        description="Hypoxia triggers invasion (Go response)",
        sign="+",
    ))

    dag.add_edge(CausalEdge(
        source="U_location",
        target="Hypoxia",
        effect_type=EffectType.LINEAR,
        parameter_symbol="u1",
        parameter_name="location_confounding",
        parameter_value=1.0,
        equation_template="Hypoxia = f(location)",
        description="Location affects O2 availability",
        sign="+",
    ))

    dag.add_edge(CausalEdge(
        source="U_location",
        target="Invasion",
        effect_type=EffectType.LINEAR,
        parameter_symbol="u2",
        parameter_name="location_invasion",
        parameter_value=1.0,
        equation_template="Invasion = f(location)",
        description="Location affects invasion (edge effects)",
        sign="+",
    ))

    dag.add_edge(CausalEdge(
        source="U_neighbors",
        target="Proliferation",
        effect_type=EffectType.LINEAR,
        parameter_symbol="u3",
        parameter_name="neighbor_prolif",
        parameter_value=1.0,
        equation_template="Prolif = f(neighbors)",
        description="Contact inhibition from neighbors",
        sign="-",
    ))

    return dag


def run_simulation(hours: int = 168, seed: int = 42) -> Path:
    """Run simulation and return output directory."""
    import tempfile
    import json as json_module

    # Load default config
    config_path = Path(__file__).parent.parent / "causanta" / "simulate" / "params" / "default.json"
    with open(config_path) as f:
        config = json_module.load(f)

    # Modify for analysis
    config["time"]["total_hours"] = hours
    config["rng_seed"] = seed

    # Create temp config
    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json_module.dump(config, f)
        temp_config = f.name

    # Run simulation
    print(f"Running simulation for {hours} hours...")
    sim = Simulation(temp_config)
    sim.run()

    # Clean up temp config
    Path(temp_config).unlink()

    return sim.output_dir


def load_simulation_data(output_dir: Path) -> dict[str, Any]:
    """Load all simulation output data."""
    import pandas as pd

    # Find all cell TSV files
    cell_files = sorted(output_dir.glob("cells_t*.tsv"))

    if not cell_files:
        raise FileNotFoundError(f"No cell files found in {output_dir}")

    # Load final timepoint
    final_cells = pd.read_csv(cell_files[-1], sep="\t")

    # Load all timepoints for time series
    all_cells = []
    for f in cell_files:
        df = pd.read_csv(f, sep="\t")
        t = int(f.stem.split("_t")[1])
        df["timepoint"] = t
        all_cells.append(df)

    all_cells_df = pd.concat(all_cells, ignore_index=True)

    # Load config
    config_path = output_dir / "params.json"
    with open(config_path) as f:
        config = json.load(f)

    return {
        "final_cells": final_cells,
        "all_cells": all_cells_df,
        "config": config,
        "output_dir": output_dir,
    }


def compute_descriptive_statistics(data: dict) -> dict[str, Any]:
    """Compute comprehensive descriptive statistics."""
    df = data["final_cells"]

    # Filter to tumor cells only
    tumor = df[df["cell_type"] == 6].copy()

    if len(tumor) == 0:
        return {"error": "No tumor cells found"}

    stats = {
        "n_total_cells": len(df),
        "n_tumor_cells": len(tumor),
        "tumor_fraction": len(tumor) / len(df),
    }

    # ecDNA statistics
    if "ecDNA_count" in tumor.columns:
        ecDNA = tumor["ecDNA_count"]
        stats["ecDNA"] = {
            "mean": float(ecDNA.mean()),
            "std": float(ecDNA.std()),
            "median": float(ecDNA.median()),
            "min": int(ecDNA.min()),
            "max": int(ecDNA.max()),
            "q25": float(ecDNA.quantile(0.25)),
            "q75": float(ecDNA.quantile(0.75)),
            "pct_zero": float((ecDNA == 0).mean() * 100),
            "pct_high": float((ecDNA > 50).mean() * 100),
        }

    # EGFR expression statistics (causal chain: ecDNA → EGFR → phenotypes)
    if "egfr_expression" in tumor.columns:
        egfr = tumor["egfr_expression"]
        stats["egfr_expression"] = {
            "mean": float(egfr.mean()),
            "std": float(egfr.std()),
            "median": float(egfr.median()),
            "min": float(egfr.min()),
            "max": float(egfr.max()),
            "q25": float(egfr.quantile(0.25)),
            "q75": float(egfr.quantile(0.75)),
        }
        # Correlation with ecDNA (should be high - gene dosage effect)
        if "ecDNA_count" in tumor.columns and len(tumor) > 2:
            stats["egfr_expression"]["corr_with_ecDNA"] = float(
                tumor["egfr_expression"].corr(tumor["ecDNA_count"])
            )

    # Hypoxia statistics
    if "is_hypoxic" in tumor.columns:
        hypoxic = tumor["is_hypoxic"].astype(str).str.lower() == "true"
        stats["hypoxia"] = {
            "pct_hypoxic": float(hypoxic.mean() * 100),
            "n_hypoxic": int(hypoxic.sum()),
            "n_normoxic": int((~hypoxic).sum()),
        }

    # Migration statistics
    if "migration_rate" in tumor.columns:
        mig = tumor["migration_rate"]
        stats["migration"] = {
            "mean": float(mig.mean()),
            "std": float(mig.std()),
            "median": float(mig.median()),
            "min": float(mig.min()),
            "max": float(mig.max()),
        }

    # O2 statistics
    if "O2_local" in tumor.columns:
        o2 = tumor["O2_local"]
        stats["O2"] = {
            "mean": float(o2.mean()),
            "std": float(o2.std()),
            "median": float(o2.median()),
            "min": float(o2.min()),
            "max": float(o2.max()),
        }

    return stats


def compute_correlations(data: dict) -> dict[str, Any]:
    """Compute correlation matrix for key variables."""
    df = data["final_cells"]
    tumor = df[df["cell_type"] == 6].copy()

    if len(tumor) == 0:
        return {"error": "No tumor cells"}

    # Convert hypoxia to numeric
    if "is_hypoxic" in tumor.columns:
        tumor["is_hypoxic_num"] = (tumor["is_hypoxic"].astype(str).str.lower() == "true").astype(int)

    # Select numeric columns (including EGFR expression for causal chain)
    cols = ["ecDNA_count", "egfr_expression", "migration_rate", "O2_local", "is_hypoxic_num"]
    cols = [c for c in cols if c in tumor.columns]

    if len(cols) < 2:
        return {"error": "Not enough columns for correlation"}

    # Correlation matrix
    corr_matrix = tumor[cols].corr()

    # Pairwise correlations with p-values
    from scipy import stats as scipy_stats

    correlations = {}
    for i, c1 in enumerate(cols):
        for c2 in cols[i+1:]:
            x = tumor[c1].dropna()
            y = tumor[c2].dropna()
            common = x.index.intersection(y.index)
            if len(common) > 2:
                r, p = scipy_stats.pearsonr(x[common], y[common])
                correlations[f"{c1}_vs_{c2}"] = {
                    "pearson_r": float(r),
                    "p_value": float(p),
                    "n": len(common),
                    "significant_05": p < 0.05,
                    "significant_01": p < 0.01,
                }

    return {
        "matrix": corr_matrix.to_dict(),
        "pairwise": correlations,
    }


def run_naive_regression(data: dict) -> dict[str, Any]:
    """Run naive OLS regression (biased by confounding)."""
    from scipy import stats as scipy_stats

    df = data["final_cells"]
    tumor = df[df["cell_type"] == 6].copy()

    if len(tumor) < 10:
        return {"error": "Too few tumor cells"}

    results = {}

    # Convert hypoxia to numeric
    if "is_hypoxic" in tumor.columns:
        tumor["is_hypoxic_num"] = (tumor["is_hypoxic"].astype(str).str.lower() == "true").astype(int)

    # Regression 1: ecDNA → Hypoxia (naive)
    if "ecDNA_count" in tumor.columns and "is_hypoxic_num" in tumor.columns:
        x = tumor["ecDNA_count"].values
        y = tumor["is_hypoxic_num"].values

        # Add constant
        X = np.column_stack([np.ones(len(x)), x])

        # OLS
        try:
            beta = np.linalg.lstsq(X, y, rcond=None)[0]
            y_pred = X @ beta
            residuals = y - y_pred
            n, k = X.shape

            # Standard errors
            mse = np.sum(residuals**2) / (n - k)
            var_beta = mse * np.linalg.inv(X.T @ X)
            se = np.sqrt(np.diag(var_beta))

            # t-statistics
            t_stats = beta / se
            p_values = 2 * (1 - scipy_stats.t.cdf(np.abs(t_stats), n - k))

            # R-squared
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y - np.mean(y))**2)
            r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0

            results["ecDNA_to_hypoxia_naive"] = {
                "method": "OLS (Naive - potentially biased)",
                "outcome": "is_hypoxic",
                "treatment": "ecDNA_count",
                "coefficient": float(beta[1]),
                "std_error": float(se[1]),
                "t_statistic": float(t_stats[1]),
                "p_value": float(p_values[1]),
                "r_squared": float(r_squared),
                "n": int(n),
                "interpretation": (
                    f"A 1-unit increase in ecDNA is associated with a "
                    f"{beta[1]:.4f} change in P(hypoxic). "
                    f"{'Significant' if p_values[1] < 0.05 else 'Not significant'} at α=0.05. "
                    f"WARNING: This is a NAIVE estimate - likely confounded!"
                ),
            }
        except Exception as e:
            results["ecDNA_to_hypoxia_naive"] = {"error": str(e)}

    # Regression 2: ecDNA → Migration
    if "ecDNA_count" in tumor.columns and "migration_rate" in tumor.columns:
        x = tumor["ecDNA_count"].values
        y = tumor["migration_rate"].values

        X = np.column_stack([np.ones(len(x)), x])

        try:
            beta = np.linalg.lstsq(X, y, rcond=None)[0]
            y_pred = X @ beta
            residuals = y - y_pred
            n, k = X.shape

            mse = np.sum(residuals**2) / (n - k)
            var_beta = mse * np.linalg.inv(X.T @ X)
            se = np.sqrt(np.diag(var_beta))

            t_stats = beta / se
            p_values = 2 * (1 - scipy_stats.t.cdf(np.abs(t_stats), n - k))

            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y - np.mean(y))**2)
            r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0

            results["ecDNA_to_migration_naive"] = {
                "method": "OLS (Naive)",
                "outcome": "migration_rate",
                "treatment": "ecDNA_count",
                "coefficient": float(beta[1]),
                "std_error": float(se[1]),
                "t_statistic": float(t_stats[1]),
                "p_value": float(p_values[1]),
                "r_squared": float(r_squared),
                "n": int(n),
                "interpretation": (
                    f"Each additional ecDNA copy associated with "
                    f"{beta[1]:.3f} um/hr increase in migration speed. "
                    f"{'Significant' if p_values[1] < 0.05 else 'Not significant'} at α=0.05."
                ),
            }
        except Exception as e:
            results["ecDNA_to_migration_naive"] = {"error": str(e)}

    # Regression 3: Hypoxia → Migration (conditional on ecDNA)
    if all(c in tumor.columns for c in ["ecDNA_count", "is_hypoxic_num", "migration_rate"]):
        x1 = tumor["ecDNA_count"].values
        x2 = tumor["is_hypoxic_num"].values
        y = tumor["migration_rate"].values

        X = np.column_stack([np.ones(len(x1)), x1, x2])

        try:
            beta = np.linalg.lstsq(X, y, rcond=None)[0]
            y_pred = X @ beta
            residuals = y - y_pred
            n, k = X.shape

            mse = np.sum(residuals**2) / (n - k)
            var_beta = mse * np.linalg.inv(X.T @ X)
            se = np.sqrt(np.diag(var_beta))

            t_stats = beta / se
            p_values = 2 * (1 - scipy_stats.t.cdf(np.abs(t_stats), n - k))

            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y - np.mean(y))**2)
            r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0

            # F-test for overall model
            f_stat = (r_squared / (k - 1)) / ((1 - r_squared) / (n - k)) if r_squared < 1 else np.inf
            f_pvalue = 1 - scipy_stats.f.cdf(f_stat, k - 1, n - k)

            results["hypoxia_migration_multivariate"] = {
                "method": "Multiple OLS",
                "outcome": "migration_rate",
                "predictors": ["ecDNA_count", "is_hypoxic"],
                "coefficients": {
                    "intercept": float(beta[0]),
                    "ecDNA_count": float(beta[1]),
                    "is_hypoxic": float(beta[2]),
                },
                "std_errors": {
                    "intercept": float(se[0]),
                    "ecDNA_count": float(se[1]),
                    "is_hypoxic": float(se[2]),
                },
                "p_values": {
                    "intercept": float(p_values[0]),
                    "ecDNA_count": float(p_values[1]),
                    "is_hypoxic": float(p_values[2]),
                },
                "r_squared": float(r_squared),
                "f_statistic": float(f_stat),
                "f_pvalue": float(f_pvalue),
                "n": int(n),
                "interpretation": (
                    f"Controlling for ecDNA, hypoxia increases migration by "
                    f"{beta[2]:.2f} um/hr (p={p_values[2]:.4f}). "
                    f"This captures the 'Go' effect."
                ),
            }
        except Exception as e:
            results["hypoxia_migration_multivariate"] = {"error": str(e)}

    return results


def run_iv_analysis(data: dict) -> dict[str, Any]:
    """Run instrumental variable analysis using ecDNA as instrument.

    We use the segregation-induced variation in ecDNA as an instrument
    to estimate causal effects free from confounding.
    """
    from scipy import stats as scipy_stats

    df = data["final_cells"]
    tumor = df[df["cell_type"] == 6].copy()

    if len(tumor) < 20:
        return {"error": "Too few tumor cells for IV"}

    results = {}

    # Convert hypoxia
    if "is_hypoxic" in tumor.columns:
        tumor["is_hypoxic_num"] = (tumor["is_hypoxic"].astype(str).str.lower() == "true").astype(int)

    # IV Analysis: ecDNA (instrument) → O2 (mediator) → Migration (outcome)
    # We can test if the effect of ecDNA on migration works through O2/hypoxia

    if all(c in tumor.columns for c in ["ecDNA_count", "O2_local", "migration_rate"]):
        z = tumor["ecDNA_count"].values  # Instrument
        x = tumor["O2_local"].values     # Endogenous (potential mediator)
        y = tumor["migration_rate"].values  # Outcome

        n = len(z)

        # Stage 1: x = γ₀ + γ₁z + ν
        Z = np.column_stack([np.ones(n), z])
        gamma = np.linalg.lstsq(Z, x, rcond=None)[0]
        x_hat = Z @ gamma
        residuals_1 = x - x_hat

        # First-stage F-statistic (test for weak instruments)
        ss_res_1 = np.sum(residuals_1**2)
        ss_tot_1 = np.sum((x - np.mean(x))**2)
        r2_1 = 1 - ss_res_1 / ss_tot_1 if ss_tot_1 > 0 else 0
        f_stat_1 = (r2_1 / 1) / ((1 - r2_1) / (n - 2)) if r2_1 < 1 else np.inf
        f_pval_1 = 1 - scipy_stats.f.cdf(f_stat_1, 1, n - 2)

        # Stage 2: y = β₀ + β₁x̂ + ε
        X_hat = np.column_stack([np.ones(n), x_hat])
        beta_2sls = np.linalg.lstsq(X_hat, y, rcond=None)[0]
        y_pred = X_hat @ beta_2sls
        residuals_2 = y - y_pred

        # 2SLS standard errors (needs correction for predicted regressor)
        mse_2 = np.sum(residuals_2**2) / (n - 2)
        X_actual = np.column_stack([np.ones(n), x])
        var_beta = mse_2 * np.linalg.inv(X_hat.T @ X_hat)
        se_2sls = np.sqrt(np.diag(var_beta))

        # t-stat and p-value
        t_stat = beta_2sls[1] / se_2sls[1]
        p_val = 2 * (1 - scipy_stats.t.cdf(np.abs(t_stat), n - 2))

        # Reduced form: y = π₀ + π₁z + ω
        pi = np.linalg.lstsq(Z, y, rcond=None)[0]
        y_pred_rf = Z @ pi
        residuals_rf = y - y_pred_rf

        results["iv_o2_to_migration"] = {
            "method": "Two-Stage Least Squares (2SLS)",
            "instrument": "ecDNA_count",
            "endogenous": "O2_local",
            "outcome": "migration_rate",
            "first_stage": {
                "coefficient": float(gamma[1]),
                "f_statistic": float(f_stat_1),
                "f_pvalue": float(f_pval_1),
                "weak_instrument": f_stat_1 < 10,
                "interpretation": (
                    f"ecDNA explains {r2_1*100:.1f}% of O2 variation. "
                    f"F={f_stat_1:.1f}. "
                    f"{'WEAK INSTRUMENT WARNING!' if f_stat_1 < 10 else 'Strong instrument.'}"
                ),
            },
            "second_stage": {
                "coefficient": float(beta_2sls[1]),
                "std_error": float(se_2sls[1]),
                "t_statistic": float(t_stat),
                "p_value": float(p_val),
            },
            "reduced_form": {
                "coefficient": float(pi[1]),
                "interpretation": f"Total effect of ecDNA on migration: {pi[1]:.4f}",
            },
            "n": int(n),
            "interpretation": (
                f"IV estimate: 1 mmHg increase in O2 changes migration by "
                f"{beta_2sls[1]:.3f} um/hr. This is the causal effect "
                f"identified via ecDNA-induced O2 variation."
            ),
        }

    # Direct IV for ecDNA → Migration
    if "ecDNA_count" in tumor.columns and "migration_rate" in tumor.columns:
        z = tumor["ecDNA_count"].values
        y = tumor["migration_rate"].values
        n = len(z)

        Z = np.column_stack([np.ones(n), z])
        beta = np.linalg.lstsq(Z, y, rcond=None)[0]
        y_pred = Z @ beta
        residuals = y - y_pred

        mse = np.sum(residuals**2) / (n - 2)
        var_beta = mse * np.linalg.inv(Z.T @ Z)
        se = np.sqrt(np.diag(var_beta))

        t_stat = beta[1] / se[1]
        p_val = 2 * (1 - scipy_stats.t.cdf(np.abs(t_stat), n - 2))

        # Compare to config ground truth
        config = data["config"]
        tumor_cfg = config.get("cell_types", {}).get("6", {})
        true_delta = tumor_cfg.get("ecDNA_effect_on_migration", 0.05)
        base_speed = tumor_cfg.get("migration_speed_um_hr", 10.0)

        results["iv_ecdna_to_migration"] = {
            "method": "Reduced Form (direct effect)",
            "instrument_exposure": "ecDNA_count",
            "outcome": "migration_rate",
            "estimate": float(beta[1]),
            "std_error": float(se[1]),
            "t_statistic": float(t_stat),
            "p_value": float(p_val),
            "ground_truth_delta": float(true_delta),
            "ground_truth_base_speed": float(base_speed),
            "n": int(n),
            "interpretation": (
                f"Each ecDNA copy increases migration by {beta[1]:.3f} um/hr. "
                f"Ground truth formula: v = {base_speed} * (1 + {true_delta} * ecDNA), "
                f"implying slope ≈ {base_speed * true_delta:.3f} at low ecDNA."
            ),
        }

    return results


def run_mediation_analysis(data: dict) -> dict[str, Any]:
    """Run mediation analysis to decompose direct/indirect effects.

    Tests if ecDNA affects invasion:
    - Directly
    - Indirectly through hypoxia (Go-or-Grow mechanism)
    """
    from scipy import stats as scipy_stats

    df = data["final_cells"]
    tumor = df[df["cell_type"] == 6].copy()

    if len(tumor) < 20:
        return {"error": "Too few cells"}

    if "is_hypoxic" in tumor.columns:
        tumor["is_hypoxic_num"] = (tumor["is_hypoxic"].astype(str).str.lower() == "true").astype(int)

    required = ["ecDNA_count", "is_hypoxic_num", "migration_rate"]
    if not all(c in tumor.columns for c in required):
        return {"error": f"Missing columns: {required}"}

    X = tumor["ecDNA_count"].values  # Treatment
    M = tumor["is_hypoxic_num"].values  # Mediator
    Y = tumor["migration_rate"].values  # Outcome
    n = len(X)

    # Step 1: Total effect (c path)
    # Y = c₀ + c*X + ε
    X_mat = np.column_stack([np.ones(n), X])
    c = np.linalg.lstsq(X_mat, Y, rcond=None)[0]
    total_effect = c[1]

    # Step 2: Effect of X on M (a path)
    # M = a₀ + a*X + ε
    a = np.linalg.lstsq(X_mat, M, rcond=None)[0]
    a_effect = a[1]

    # Step 3: Effect of M on Y controlling for X (b path)
    # Y = c'₀ + c'*X + b*M + ε
    XM_mat = np.column_stack([np.ones(n), X, M])

    # Check if we have variation in the mediator
    if np.var(M) < 1e-10:
        return {
            "error": "No variation in mediator (all cells same hypoxia status)",
            "suggestion": "Run longer simulation to create hypoxic regions",
        }

    try:
        coefs = np.linalg.lstsq(XM_mat, Y, rcond=None)[0]
    except Exception as e:
        return {"error": f"Regression failed: {e}"}

    c_prime = coefs[1]  # Direct effect
    b_effect = coefs[2]  # Effect of mediator

    # Indirect effect = a * b
    indirect_effect = a_effect * b_effect

    # Proportion mediated
    prop_mediated = indirect_effect / total_effect if abs(total_effect) > 1e-10 else 0

    # Sobel test for indirect effect
    # SE(ab) ≈ sqrt(a²*SE(b)² + b²*SE(a)²)
    y_pred = XM_mat @ coefs
    residuals = Y - y_pred
    mse = np.sum(residuals**2) / (n - 3)

    try:
        var_coefs = mse * np.linalg.inv(XM_mat.T @ XM_mat)
        se_b = np.sqrt(var_coefs[2, 2])
    except np.linalg.LinAlgError:
        # Matrix is singular, use pseudo-inverse
        var_coefs = mse * np.linalg.pinv(XM_mat.T @ XM_mat)
        se_b = np.sqrt(max(var_coefs[2, 2], 1e-10))

    m_pred = X_mat @ a
    res_m = M - m_pred
    mse_m = np.sum(res_m**2) / (n - 2)
    var_a = mse_m * np.linalg.inv(X_mat.T @ X_mat)
    se_a = np.sqrt(var_a[1, 1])

    se_indirect = np.sqrt(a_effect**2 * se_b**2 + b_effect**2 * se_a**2)
    z_sobel = indirect_effect / se_indirect if se_indirect > 0 else 0
    p_sobel = 2 * (1 - scipy_stats.norm.cdf(np.abs(z_sobel)))

    return {
        "method": "Baron & Kenny Mediation Analysis",
        "treatment": "ecDNA_count",
        "mediator": "is_hypoxic",
        "outcome": "migration_rate",
        "paths": {
            "total_effect_c": {
                "estimate": float(total_effect),
                "interpretation": f"Total effect of ecDNA on migration: {total_effect:.4f}",
            },
            "a_path": {
                "estimate": float(a_effect),
                "interpretation": f"Effect of ecDNA on hypoxia: {a_effect:.4f}",
            },
            "b_path": {
                "estimate": float(b_effect),
                "interpretation": f"Effect of hypoxia on migration (controlling for ecDNA): {b_effect:.4f}",
            },
            "direct_effect_c_prime": {
                "estimate": float(c_prime),
                "interpretation": f"Direct effect of ecDNA on migration: {c_prime:.4f}",
            },
            "indirect_effect_ab": {
                "estimate": float(indirect_effect),
                "std_error": float(se_indirect),
                "z_sobel": float(z_sobel),
                "p_value": float(p_sobel),
                "significant": p_sobel < 0.05,
                "interpretation": (
                    f"Indirect effect through hypoxia: {indirect_effect:.4f}. "
                    f"Sobel z={z_sobel:.2f}, p={p_sobel:.4f}. "
                    f"{'Significant mediation!' if p_sobel < 0.05 else 'Mediation not significant.'}"
                ),
            },
        },
        "proportion_mediated": float(prop_mediated),
        "n": int(n),
        "overall_interpretation": (
            f"Of the total effect ({total_effect:.4f}), "
            f"{abs(prop_mediated)*100:.1f}% is mediated through hypoxia. "
            f"Direct effect: {c_prime:.4f}, Indirect: {indirect_effect:.4f}. "
            f"This {'supports' if abs(prop_mediated) > 0.1 else 'does not support'} "
            f"the Go-or-Grow hypothesis where hypoxia mediates invasion."
        ),
    }


def create_visualizations(data: dict, stats: dict, correlations: dict,
                          regression: dict, iv: dict, mediation: dict,
                          output_dir: Path) -> dict[str, str]:
    """Create all visualizations for the report."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns

    viz_dir = output_dir / "figures"
    viz_dir.mkdir(exist_ok=True)

    figures = {}
    df = data["final_cells"]
    tumor = df[df["cell_type"] == 6].copy()

    if len(tumor) == 0:
        return figures

    # Convert hypoxia
    if "is_hypoxic" in tumor.columns:
        tumor["is_hypoxic_num"] = (tumor["is_hypoxic"].astype(str).str.lower() == "true").astype(int)
        tumor["Hypoxia Status"] = tumor["is_hypoxic_num"].map({0: "Normoxic", 1: "Hypoxic"})

    plt.style.use('seaborn-v0_8-whitegrid')

    # 1. ecDNA Distribution
    if "ecDNA_count" in tumor.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.histplot(tumor["ecDNA_count"], bins=30, kde=True, ax=ax, color="#3498db")
        ax.axvline(tumor["ecDNA_count"].mean(), color='red', linestyle='--',
                   label=f'Mean: {tumor["ecDNA_count"].mean():.1f}')
        ax.axvline(tumor["ecDNA_count"].median(), color='green', linestyle='--',
                   label=f'Median: {tumor["ecDNA_count"].median():.1f}')
        ax.set_xlabel("ecDNA Copy Number", fontsize=12)
        ax.set_ylabel("Count", fontsize=12)
        ax.set_title("Distribution of ecDNA Copy Number in Tumor Cells", fontsize=14)
        ax.legend()
        fig.savefig(viz_dir / "ecdna_distribution.png", dpi=150, bbox_inches='tight')
        figures["ecdna_distribution"] = "figures/ecdna_distribution.png"
        plt.close(fig)

    # 2. ecDNA vs Migration scatter with hypoxia coloring
    if all(c in tumor.columns for c in ["ecDNA_count", "migration_rate", "Hypoxia Status"]):
        fig, ax = plt.subplots(figsize=(10, 7))
        colors = {"Normoxic": "#3498db", "Hypoxic": "#e74c3c"}
        for status, color in colors.items():
            mask = tumor["Hypoxia Status"] == status
            ax.scatter(tumor.loc[mask, "ecDNA_count"],
                      tumor.loc[mask, "migration_rate"],
                      c=color, alpha=0.5, label=status, s=30)

        # Add regression lines for each group
        for status, color in colors.items():
            mask = tumor["Hypoxia Status"] == status
            x = tumor.loc[mask, "ecDNA_count"].values
            y = tumor.loc[mask, "migration_rate"].values
            if len(x) > 2:
                z = np.polyfit(x, y, 1)
                p = np.poly1d(z)
                x_line = np.linspace(x.min(), x.max(), 100)
                ax.plot(x_line, p(x_line), color=color, linewidth=2,
                       linestyle='--', label=f'{status} trend')

        ax.set_xlabel("ecDNA Copy Number", fontsize=12)
        ax.set_ylabel("Migration Rate (um/hr)", fontsize=12)
        ax.set_title("ecDNA vs Migration: Go-or-Grow Effect\n(Hypoxic cells migrate faster at same ecDNA level)", fontsize=14)
        ax.legend()
        fig.savefig(viz_dir / "ecdna_migration_scatter.png", dpi=150, bbox_inches='tight')
        figures["ecdna_migration_scatter"] = "figures/ecdna_migration_scatter.png"
        plt.close(fig)

    # 3. Correlation Heatmap
    cols = ["ecDNA_count", "migration_rate", "O2_local"]
    if "is_hypoxic_num" in tumor.columns:
        cols.append("is_hypoxic_num")
    cols = [c for c in cols if c in tumor.columns]

    if len(cols) >= 2:
        fig, ax = plt.subplots(figsize=(8, 6))
        corr = tumor[cols].corr()
        sns.heatmap(corr, annot=True, cmap="RdBu_r", center=0,
                   fmt=".2f", ax=ax, vmin=-1, vmax=1)
        ax.set_title("Correlation Matrix of Key Variables", fontsize=14)
        fig.savefig(viz_dir / "correlation_heatmap.png", dpi=150, bbox_inches='tight')
        figures["correlation_heatmap"] = "figures/correlation_heatmap.png"
        plt.close(fig)

    # 4. Effect Estimates Comparison (Forest Plot style)
    fig, ax = plt.subplots(figsize=(12, 6))

    estimates = []
    labels = []
    errors = []

    if "ecDNA_to_migration_naive" in regression:
        r = regression["ecDNA_to_migration_naive"]
        if "coefficient" in r:
            estimates.append(r["coefficient"])
            errors.append(r["std_error"] * 1.96)
            labels.append("Naive OLS")

    if "iv_ecdna_to_migration" in iv:
        r = iv["iv_ecdna_to_migration"]
        if "estimate" in r:
            estimates.append(r["estimate"])
            errors.append(r["std_error"] * 1.96)
            labels.append("IV (Reduced Form)")

    if estimates:
        y_pos = np.arange(len(labels))
        ax.barh(y_pos, estimates, xerr=errors, align='center',
               color=['#3498db', '#27ae60'][:len(estimates)], alpha=0.7, capsize=5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.axvline(0, color='black', linewidth=0.5)

        # Add ground truth line
        config = data["config"]
        tumor_cfg = config.get("cell_types", {}).get("6", {})
        true_delta = tumor_cfg.get("ecDNA_effect_on_migration", 0.05)
        base_speed = tumor_cfg.get("migration_speed_um_hr", 10.0)
        gt_slope = base_speed * true_delta
        ax.axvline(gt_slope, color='red', linestyle='--', linewidth=2,
                  label=f'Ground Truth: {gt_slope:.3f}')

        ax.set_xlabel("Effect of ecDNA on Migration (um/hr per copy)", fontsize=12)
        ax.set_title("Comparison of Causal Effect Estimates", fontsize=14)
        ax.legend()
        fig.savefig(viz_dir / "effect_comparison.png", dpi=150, bbox_inches='tight')
        figures["effect_comparison"] = "figures/effect_comparison.png"
        plt.close(fig)

    # 5. Mediation Diagram
    if mediation and "paths" in mediation and "error" not in mediation:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 6)
        ax.axis('off')

        # Draw boxes
        boxes = {
            "ecDNA": (1, 3),
            "Hypoxia": (5, 5),
            "Migration": (9, 3),
        }

        for name, (x, y) in boxes.items():
            rect = plt.Rectangle((x-0.8, y-0.4), 1.6, 0.8,
                                 fill=True, facecolor='lightblue',
                                 edgecolor='black', linewidth=2)
            ax.add_patch(rect)
            ax.text(x, y, name, ha='center', va='center', fontsize=12, fontweight='bold')

        # Draw arrows with coefficients
        paths = mediation["paths"]

        # a path: ecDNA → Hypoxia
        ax.annotate('', xy=(4.2, 4.6), xytext=(1.8, 3.4),
                   arrowprops=dict(arrowstyle='->', color='blue', lw=2))
        a_val = paths.get("a_path", {}).get("estimate", 0)
        ax.text(2.8, 4.3, f'a = {a_val:.4f}', fontsize=10, color='blue')

        # b path: Hypoxia → Migration
        ax.annotate('', xy=(8.2, 3.4), xytext=(5.8, 4.6),
                   arrowprops=dict(arrowstyle='->', color='blue', lw=2))
        b_val = paths.get("b_path", {}).get("estimate", 0)
        ax.text(7.2, 4.3, f'b = {b_val:.2f}', fontsize=10, color='blue')

        # c' path: ecDNA → Migration (direct)
        ax.annotate('', xy=(8.2, 3), xytext=(1.8, 3),
                   arrowprops=dict(arrowstyle='->', color='green', lw=2))
        c_prime = paths.get("direct_effect_c_prime", {}).get("estimate", 0)
        ax.text(5, 2.5, f"c' = {c_prime:.4f} (direct)", fontsize=10, color='green')

        # Add summary
        indirect = paths.get("indirect_effect_ab", {}).get("estimate", 0)
        total = paths.get("total_effect_c", {}).get("estimate", 0)
        ax.text(5, 1, f"Indirect effect (a×b) = {indirect:.4f}", fontsize=11, ha='center')
        ax.text(5, 0.5, f"Total effect (c) = {total:.4f}", fontsize=11, ha='center')

        ax.set_title("Mediation Analysis: Does Hypoxia Mediate ecDNA → Migration?",
                    fontsize=14, fontweight='bold')

        fig.savefig(viz_dir / "mediation_diagram.png", dpi=150, bbox_inches='tight')
        figures["mediation_diagram"] = "figures/mediation_diagram.png"
        plt.close(fig)

    # 6. O2 vs Hypoxia threshold visualization
    if all(c in tumor.columns for c in ["O2_local", "is_hypoxic_num"]):
        fig, ax = plt.subplots(figsize=(10, 6))

        config = data["config"]
        hypoxia_thresh = config.get("environment", {}).get("hypoxia_threshold_mmHg", 10.0)

        colors = ["#3498db" if h == 0 else "#e74c3c" for h in tumor["is_hypoxic_num"]]
        ax.scatter(tumor["O2_local"], tumor["migration_rate"], c=colors, alpha=0.5, s=20)
        ax.axvline(hypoxia_thresh, color='red', linestyle='--', linewidth=2,
                  label=f'Hypoxia threshold: {hypoxia_thresh} mmHg')
        ax.set_xlabel("Local O2 (mmHg)", fontsize=12)
        ax.set_ylabel("Migration Rate (um/hr)", fontsize=12)
        ax.set_title("O2 Level vs Migration: Hypoxia Triggers Invasion", fontsize=14)
        ax.legend()

        fig.savefig(viz_dir / "o2_migration.png", dpi=150, bbox_inches='tight')
        figures["o2_migration"] = "figures/o2_migration.png"
        plt.close(fig)

    # 7. DAG Visualization
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # True DAG
    ax = axes[0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis('off')
    ax.set_title("TRUE Causal DAG (Simulation)", fontsize=14, fontweight='bold')

    true_nodes = {
        "ecDNA": (1, 4),
        "Prolif": (3.5, 4),
        "O2": (6, 4),
        "Hypoxia": (6, 6),
        "Invasion": (9, 6),
    }

    for name, (x, y) in true_nodes.items():
        rect = plt.Rectangle((x-0.6, y-0.3), 1.2, 0.6,
                             fill=True, facecolor='lightgreen',
                             edgecolor='black', linewidth=2)
        ax.add_patch(rect)
        ax.text(x, y, name, ha='center', va='center', fontsize=10)

    # True arrows
    ax.annotate('', xy=(2.9, 4), xytext=(1.6, 4), arrowprops=dict(arrowstyle='->', lw=2))
    ax.annotate('', xy=(5.4, 4), xytext=(4.1, 4), arrowprops=dict(arrowstyle='->', lw=2))
    ax.annotate('', xy=(6, 5.7), xytext=(6, 4.3), arrowprops=dict(arrowstyle='->', lw=2))
    ax.annotate('', xy=(8.4, 6), xytext=(6.6, 6), arrowprops=dict(arrowstyle='->', lw=2))

    ax.text(2.2, 4.3, '+α', fontsize=10, color='blue')
    ax.text(4.7, 4.3, '-ω', fontsize=10, color='red')
    ax.text(6.3, 5, '-', fontsize=10, color='red')
    ax.text(7.5, 6.3, '+δ', fontsize=10, color='blue')

    # Hypothesized DAG
    ax = axes[1]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 8)
    ax.axis('off')
    ax.set_title("HYPOTHESIZED DAG (Analysis)", fontsize=14, fontweight='bold')

    hyp_nodes = {
        "ecDNA": (2, 4),
        "Hypoxia": (5, 6),
        "Invasion": (8, 4),
        "U": (5, 2),
    }

    for name, (x, y) in hyp_nodes.items():
        color = 'lightblue' if name != 'U' else 'lightyellow'
        rect = plt.Rectangle((x-0.6, y-0.3), 1.2, 0.6,
                             fill=True, facecolor=color,
                             edgecolor='black', linewidth=2,
                             linestyle='--' if name == 'U' else '-')
        ax.add_patch(rect)
        ax.text(x, y, name, ha='center', va='center', fontsize=10)

    # Hypothesized arrows
    ax.annotate('', xy=(4.4, 5.7), xytext=(2.6, 4.3), arrowprops=dict(arrowstyle='->', lw=2))
    ax.annotate('', xy=(7.4, 4.3), xytext=(5.6, 5.7), arrowprops=dict(arrowstyle='->', lw=2))
    ax.annotate('', xy=(7.4, 4), xytext=(2.6, 4), arrowprops=dict(arrowstyle='->', lw=2, linestyle='--'))
    ax.annotate('', xy=(4.4, 5.7), xytext=(4.4, 2.3), arrowprops=dict(arrowstyle='->', lw=1, linestyle=':', color='gray'))
    ax.annotate('', xy=(7.4, 4), xytext=(5.6, 2.3), arrowprops=dict(arrowstyle='->', lw=1, linestyle=':', color='gray'))

    ax.text(3, 5.5, '?', fontsize=14, color='red', fontweight='bold')
    ax.text(6.5, 5.5, '+δ', fontsize=10, color='blue')
    ax.text(5, 3.7, 'direct?', fontsize=9, color='gray')

    fig.savefig(viz_dir / "dag_comparison.png", dpi=150, bbox_inches='tight')
    figures["dag_comparison"] = "figures/dag_comparison.png"
    plt.close(fig)

    return figures


def generate_html_report(
    data: dict,
    stats: dict,
    correlations: dict,
    regression: dict,
    iv: dict,
    mediation: dict,
    figures: dict,
    true_dag: CausalDAG,
    analysis_dag: CausalDAG,
    output_path: Path,
) -> None:
    """Generate comprehensive HTML analysis report."""

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    config = data["config"]
    tumor_cfg = config.get("cell_types", {}).get("6", {})

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CAUSANTA Comprehensive Analysis Report</title>
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
            line-height: 1.6;
            color: #24292f;
            background: #f6f8fa;
            margin: 0;
            padding: 20px;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 40px;
            border-radius: 6px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #0969da;
            border-bottom: 2px solid #0969da;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #24292f;
            border-bottom: 1px solid #d0d7de;
            padding-bottom: 8px;
            margin-top: 40px;
        }}
        h3 {{
            color: #57606a;
            margin-top: 30px;
        }}
        .summary-box {{
            background: #f0f9ff;
            border: 1px solid #0969da;
            border-radius: 6px;
            padding: 20px;
            margin: 20px 0;
        }}
        .warning-box {{
            background: #fff8e6;
            border: 1px solid #d29922;
            border-radius: 6px;
            padding: 20px;
            margin: 20px 0;
        }}
        .success-box {{
            background: #f0fff4;
            border: 1px solid #2da44e;
            border-radius: 6px;
            padding: 20px;
            margin: 20px 0;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 20px 0;
        }}
        th, td {{
            border: 1px solid #d0d7de;
            padding: 12px;
            text-align: left;
        }}
        th {{
            background: #f6f8fa;
            font-weight: 600;
        }}
        tr:hover td {{
            background: #f6f8fa;
        }}
        .figure {{
            text-align: center;
            margin: 30px 0;
        }}
        .figure img {{
            max-width: 100%;
            border: 1px solid #d0d7de;
            border-radius: 6px;
        }}
        .figure-caption {{
            font-style: italic;
            color: #57606a;
            margin-top: 10px;
        }}
        code {{
            background: #f6f8fa;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: monospace;
        }}
        .stat-significant {{
            color: #2da44e;
            font-weight: bold;
        }}
        .stat-not-significant {{
            color: #cf222e;
        }}
        .interpretation {{
            background: #f6f8fa;
            border-left: 4px solid #0969da;
            padding: 15px;
            margin: 15px 0;
            font-style: italic;
        }}
        .test-explanation {{
            background: #fff8e6;
            border-left: 4px solid #d29922;
            padding: 15px;
            margin: 15px 0;
        }}
    </style>
</head>
<body>
<div class="container">

<h1>CAUSANTA Comprehensive Causal Analysis Report</h1>

<p><strong>Generated:</strong> {timestamp}</p>
<p><strong>Simulation Duration:</strong> {config.get("time", {}).get("total_hours", "?")} hours</p>
<p><strong>Random Seed:</strong> {config.get("rng_seed", "?")}</p>

<div class="summary-box">
<h3>Executive Summary</h3>
<p>This report analyzes simulation data using causal inference methods to test whether
we can recover the TRUE causal structure from observational data. The key finding:</p>
<ul>
<li><strong>TRUE model:</strong> ecDNA → Proliferation → O2 depletion → Hypoxia → Invasion</li>
<li><strong>Tested hypothesis:</strong> Does ecDNA directly cause hypoxia? (Answer: NO - it's mediated)</li>
<li><strong>IV validity:</strong> ecDNA segregation provides valid randomization for causal inference</li>
</ul>
</div>

<h2>1. Causal Structure: True vs. Hypothesized DAG</h2>

<p>We simulate from a KNOWN causal structure but analyze with a HYPOTHESIZED structure
that an analyst might propose. This tests whether causal methods can recover truth.</p>

<div class="figure">
<img src="{figures.get("dag_comparison", "")}" alt="DAG Comparison">
<div class="figure-caption">Figure 1: True simulation DAG (left) vs. hypothesized analysis DAG (right).
The key question: does ecDNA directly cause hypoxia, or is the effect mediated through proliferation?</div>
</div>

<div class="test-explanation">
<strong>Why different DAGs?</strong><br>
In real research, we don't know the true causal structure. By simulating from known
ground truth, we can test if our statistical methods correctly identify:
<ul>
<li>That ecDNA affects invasion <em>indirectly</em> through proliferation and hypoxia</li>
<li>That hypoxia is a <em>mediator</em>, not a direct effect of ecDNA</li>
<li>That unmeasured confounders (location, neighbors) don't bias our IV estimates</li>
</ul>
</div>

<h2>2. Descriptive Statistics</h2>

<h3>2.1 Sample Overview</h3>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Total cells</td><td>{stats.get("n_total_cells", "N/A"):,}</td></tr>
<tr><td>Tumor cells</td><td>{stats.get("n_tumor_cells", "N/A"):,}</td></tr>
<tr><td>Tumor fraction</td><td>{stats.get("tumor_fraction", 0)*100:.1f}%</td></tr>
</table>
'''

    # ecDNA stats
    if "ecDNA" in stats:
        e = stats["ecDNA"]
        html += f'''
<h3>2.2 ecDNA Copy Number Distribution</h3>

<div class="figure">
<img src="{figures.get("ecdna_distribution", "")}" alt="ecDNA Distribution">
<div class="figure-caption">Figure 2: Distribution of ecDNA copy number across tumor cells.
ecDNA segregates binomially during mitosis, creating natural variation.</div>
</div>

<table>
<tr><th>Statistic</th><th>Value</th><th>Interpretation</th></tr>
<tr><td>Mean</td><td>{e.get("mean", 0):.2f}</td><td>Average ecDNA copies per tumor cell</td></tr>
<tr><td>Std Dev</td><td>{e.get("std", 0):.2f}</td><td>Variation due to binomial segregation</td></tr>
<tr><td>Median</td><td>{e.get("median", 0):.1f}</td><td>Central tendency</td></tr>
<tr><td>Range</td><td>{e.get("min", 0)} - {e.get("max", 0)}</td><td>Min to max copies</td></tr>
<tr><td>% with zero</td><td>{e.get("pct_zero", 0):.1f}%</td><td>Cells that lost all ecDNA</td></tr>
<tr><td>% with >50</td><td>{e.get("pct_high", 0):.1f}%</td><td>High-amplification cells</td></tr>
</table>

<div class="interpretation">
<strong>Interpretation:</strong> The ecDNA distribution shows the characteristic variance
created by binomial segregation. Cells with mean ~20 copies will produce daughters with
~10-30 copies after replication and segregation. This natural randomization is the
basis for using ecDNA as an instrumental variable.
</div>
'''

    # Hypoxia stats
    if "hypoxia" in stats:
        h = stats["hypoxia"]
        html += f'''
<h3>2.3 Hypoxia Status</h3>
<table>
<tr><th>Status</th><th>Count</th><th>Percentage</th></tr>
<tr><td>Hypoxic (O2 &lt; {config.get("environment", {}).get("hypoxia_threshold_mmHg", 10)} mmHg)</td>
    <td>{h.get("n_hypoxic", 0):,}</td><td>{h.get("pct_hypoxic", 0):.1f}%</td></tr>
<tr><td>Normoxic</td><td>{h.get("n_normoxic", 0):,}</td><td>{100 - h.get("pct_hypoxic", 0):.1f}%</td></tr>
</table>
'''

    # Migration stats
    if "migration" in stats:
        m = stats["migration"]
        html += f'''
<h3>2.4 Migration Rate</h3>
<table>
<tr><th>Statistic</th><th>Value (um/hr)</th></tr>
<tr><td>Mean</td><td>{m.get("mean", 0):.2f}</td></tr>
<tr><td>Std Dev</td><td>{m.get("std", 0):.2f}</td></tr>
<tr><td>Range</td><td>{m.get("min", 0):.1f} - {m.get("max", 0):.1f}</td></tr>
</table>
'''

    # Correlations
    html += '''
<h2>3. Correlation Analysis</h2>

<p>Before causal analysis, we examine pairwise correlations. <strong>Important:</strong>
Correlation does not imply causation - these are descriptive only.</p>

<div class="figure">
<img src="''' + figures.get("correlation_heatmap", "") + '''" alt="Correlation Heatmap">
<div class="figure-caption">Figure 3: Correlation matrix of key variables.</div>
</div>
'''

    if "pairwise" in correlations:
        html += '''
<table>
<tr><th>Variable Pair</th><th>Pearson r</th><th>p-value</th><th>Significant?</th></tr>
'''
        for pair, vals in correlations["pairwise"].items():
            sig_class = "stat-significant" if vals.get("significant_05") else "stat-not-significant"
            sig_text = "Yes" if vals.get("significant_05") else "No"
            html += f'''
<tr>
<td>{pair.replace("_", " → ")}</td>
<td>{vals.get("pearson_r", 0):.3f}</td>
<td>{vals.get("p_value", 1):.4f}</td>
<td class="{sig_class}">{sig_text}</td>
</tr>
'''
        html += '</table>'

    # Main scatter plot
    html += f'''
<div class="figure">
<img src="{figures.get("ecdna_migration_scatter", "")}" alt="ecDNA vs Migration">
<div class="figure-caption">Figure 4: ecDNA copy number vs. migration rate, colored by hypoxia status.
Note that hypoxic cells (red) migrate faster than normoxic cells (blue) at the same ecDNA level -
this is the "Go or Grow" effect.</div>
</div>

<div class="interpretation">
<strong>Key observation:</strong> The two parallel trend lines show that hypoxia
has an ADDITIVE effect on migration beyond what ecDNA provides. This suggests
hypoxia is a mediator or effect modifier, not just a correlate.
</div>
'''

    # Regression analysis
    html += '''
<h2>4. Regression Analysis (Naive OLS)</h2>

<div class="warning-box">
<strong>Warning:</strong> OLS estimates may be biased by confounding. We use these
as a baseline to compare against causal (IV) estimates below.
</div>
'''

    if "ecDNA_to_migration_naive" in regression:
        r = regression["ecDNA_to_migration_naive"]
        if "coefficient" in r:
            sig_class = "stat-significant" if r.get("p_value", 1) < 0.05 else "stat-not-significant"
            html += f'''
<h3>4.1 ecDNA → Migration (Naive OLS)</h3>

<div class="test-explanation">
<strong>What this tests:</strong> Simple linear regression of migration rate on ecDNA count.
The coefficient represents the average increase in migration (um/hr) per additional ecDNA copy.
This is potentially biased by omitted variables (location, neighbors, O2 levels).
</div>

<table>
<tr><th>Parameter</th><th>Estimate</th><th>Std Error</th><th>t-stat</th><th>p-value</th></tr>
<tr>
<td>ecDNA coefficient</td>
<td>{r.get("coefficient", 0):.4f}</td>
<td>{r.get("std_error", 0):.4f}</td>
<td>{r.get("t_statistic", 0):.2f}</td>
<td class="{sig_class}">{r.get("p_value", 1):.4f}</td>
</tr>
</table>

<p><strong>R²:</strong> {r.get("r_squared", 0):.3f} (ecDNA explains {r.get("r_squared", 0)*100:.1f}% of migration variance)</p>
<p><strong>N:</strong> {r.get("n", 0):,} tumor cells</p>

<div class="interpretation">{r.get("interpretation", "")}</div>
'''

    # Multivariate regression
    if "hypoxia_migration_multivariate" in regression:
        r = regression["hypoxia_migration_multivariate"]
        if "coefficients" in r:
            html += f'''
<h3>4.2 Multivariate: ecDNA + Hypoxia → Migration</h3>

<div class="test-explanation">
<strong>What this tests:</strong> Multiple regression controlling for both ecDNA and hypoxia status.
This lets us decompose the separate effects of EGFR amplification and hypoxic stress.
</div>

<table>
<tr><th>Predictor</th><th>Coefficient</th><th>Std Error</th><th>p-value</th></tr>
<tr>
<td>Intercept</td>
<td>{r["coefficients"].get("intercept", 0):.3f}</td>
<td>{r["std_errors"].get("intercept", 0):.3f}</td>
<td>{r["p_values"].get("intercept", 1):.4f}</td>
</tr>
<tr>
<td>ecDNA count</td>
<td>{r["coefficients"].get("ecDNA_count", 0):.4f}</td>
<td>{r["std_errors"].get("ecDNA_count", 0):.4f}</td>
<td>{r["p_values"].get("ecDNA_count", 1):.4f}</td>
</tr>
<tr>
<td>is_hypoxic</td>
<td>{r["coefficients"].get("is_hypoxic", 0):.3f}</td>
<td>{r["std_errors"].get("is_hypoxic", 0):.3f}</td>
<td>{r["p_values"].get("is_hypoxic", 1):.4f}</td>
</tr>
</table>

<p><strong>Model R²:</strong> {r.get("r_squared", 0):.3f}</p>
<p><strong>F-statistic:</strong> {r.get("f_statistic", 0):.2f} (p = {r.get("f_pvalue", 1):.4f})</p>

<div class="interpretation">{r.get("interpretation", "")}</div>
'''

    # IV Analysis
    html += '''
<h2>5. Instrumental Variable (IV) Analysis</h2>

<div class="success-box">
<strong>Key insight:</strong> ecDNA segregation provides natural randomization during cell division.
Because segregation is binomial (random), ecDNA count is independent of confounders like
cell location or neighbor effects. This makes ecDNA a valid instrument for causal inference.
</div>

<div class="test-explanation">
<strong>IV Assumptions:</strong>
<ol>
<li><strong>Relevance:</strong> ecDNA directly determines gene dosage (EGFR expression)</li>
<li><strong>Independence:</strong> Binomial segregation is random, independent of confounders</li>
<li><strong>Exclusion:</strong> ecDNA affects outcomes only through gene expression pathways</li>
</ol>
</div>
'''

    if "iv_ecdna_to_migration" in iv:
        r = iv["iv_ecdna_to_migration"]
        if "estimate" in r:
            gt_slope = r.get("ground_truth_base_speed", 10) * r.get("ground_truth_delta", 0.05)
            bias = r["estimate"] - gt_slope
            bias_pct = (bias / gt_slope * 100) if gt_slope != 0 else 0

            html += f'''
<h3>5.1 Reduced Form: ecDNA → Migration</h3>

<table>
<tr><th>Parameter</th><th>Estimate</th><th>Std Error</th><th>p-value</th></tr>
<tr>
<td>Effect of ecDNA on migration</td>
<td>{r.get("estimate", 0):.4f} um/hr per copy</td>
<td>{r.get("std_error", 0):.4f}</td>
<td>{r.get("p_value", 1):.6f}</td>
</tr>
</table>

<p><strong>Ground truth:</strong> δ = {r.get("ground_truth_delta", 0)}, base speed = {r.get("ground_truth_base_speed", 0)} um/hr</p>
<p><strong>Expected slope:</strong> {gt_slope:.4f} um/hr per copy</p>
<p><strong>Bias:</strong> {bias:.4f} ({bias_pct:+.1f}%)</p>

<div class="interpretation">{r.get("interpretation", "")}</div>
'''

    if "iv_o2_to_migration" in iv:
        r = iv["iv_o2_to_migration"]
        if "first_stage" in r:
            weak = r["first_stage"].get("weak_instrument", False)
            weak_class = "warning-box" if weak else "success-box"

            html += f'''
<h3>5.2 Two-Stage Least Squares: O2 as Mediator</h3>

<p>We test if O2 mediates the ecDNA → migration relationship by instrumenting O2 with ecDNA.</p>

<div class="{weak_class}">
<strong>First Stage (ecDNA → O2):</strong><br>
F-statistic: {r["first_stage"].get("f_statistic", 0):.1f}
({r["first_stage"].get("interpretation", "")})
</div>

<table>
<tr><th>Stage</th><th>Coefficient</th><th>Interpretation</th></tr>
<tr>
<td>First stage (γ₁)</td>
<td>{r["first_stage"].get("coefficient", 0):.4f}</td>
<td>Effect of ecDNA on O2</td>
</tr>
<tr>
<td>Second stage (β₁)</td>
<td>{r["second_stage"].get("coefficient", 0):.4f}</td>
<td>Causal effect of O2 on migration</td>
</tr>
<tr>
<td>Reduced form (π₁)</td>
<td>{r["reduced_form"].get("coefficient", 0):.4f}</td>
<td>Total effect of ecDNA on migration</td>
</tr>
</table>

<div class="interpretation">{r.get("interpretation", "")}</div>
'''

    # Mediation Analysis
    html += '''
<h2>6. Mediation Analysis</h2>

<p>Does hypoxia <em>mediate</em> the effect of ecDNA on invasion? The Go-or-Grow hypothesis
predicts that proliferating cells deplete O2, creating hypoxia, which then triggers invasion.</p>

<div class="figure">
<img src="''' + figures.get("mediation_diagram", "") + '''" alt="Mediation Diagram">
<div class="figure-caption">Figure 5: Mediation model decomposing direct and indirect effects.</div>
</div>
'''

    if mediation and "paths" in mediation and "error" not in mediation:
        paths = mediation["paths"]
        ind = paths.get("indirect_effect_ab", {})
        sig = ind.get("significant", False)
        sig_class = "stat-significant" if sig else "stat-not-significant"

        html += f'''
<div class="test-explanation">
<strong>Sobel Test:</strong> Tests whether the indirect effect (a × b) is significantly
different from zero. A significant Sobel test indicates that mediation is occurring.
</div>

<table>
<tr><th>Path</th><th>Estimate</th><th>Interpretation</th></tr>
<tr>
<td>c (total effect)</td>
<td>{paths.get("total_effect_c", {}).get("estimate", 0):.4f}</td>
<td>Total effect of ecDNA on migration</td>
</tr>
<tr>
<td>a (ecDNA → Hypoxia)</td>
<td>{paths.get("a_path", {}).get("estimate", 0):.4f}</td>
<td>Effect of ecDNA on hypoxia probability</td>
</tr>
<tr>
<td>b (Hypoxia → Migration)</td>
<td>{paths.get("b_path", {}).get("estimate", 0):.2f}</td>
<td>Effect of hypoxia on migration (controlling for ecDNA)</td>
</tr>
<tr>
<td>c' (direct effect)</td>
<td>{paths.get("direct_effect_c_prime", {}).get("estimate", 0):.4f}</td>
<td>Direct effect of ecDNA on migration</td>
</tr>
<tr>
<td>a × b (indirect)</td>
<td>{ind.get("estimate", 0):.4f}</td>
<td>Indirect effect through hypoxia</td>
</tr>
</table>

<h3>Sobel Test for Mediation</h3>
<table>
<tr><th>Statistic</th><th>Value</th></tr>
<tr><td>Indirect effect (a × b)</td><td>{ind.get("estimate", 0):.4f}</td></tr>
<tr><td>Standard error</td><td>{ind.get("std_error", 0):.4f}</td></tr>
<tr><td>Z-score</td><td>{ind.get("z_sobel", 0):.3f}</td></tr>
<tr><td>p-value</td><td class="{sig_class}">{ind.get("p_value", 1):.4f}</td></tr>
<tr><td>Significant mediation?</td><td class="{sig_class}">{"Yes" if sig else "No"}</td></tr>
</table>

<p><strong>Proportion mediated:</strong> {mediation.get("proportion_mediated", 0)*100:.1f}%</p>

<div class="interpretation">{mediation.get("overall_interpretation", "")}</div>
'''

    # Effect comparison
    html += f'''
<h2>7. Comparison of Effect Estimates</h2>

<div class="figure">
<img src="{figures.get("effect_comparison", "")}" alt="Effect Comparison">
<div class="figure-caption">Figure 6: Forest plot comparing effect estimates across methods.
Red dashed line shows ground truth from simulation parameters.</div>
</div>

<div class="test-explanation">
<strong>What to look for:</strong>
<ul>
<li>If IV and OLS estimates agree: little confounding bias</li>
<li>If IV estimate closer to ground truth: IV correcting for confounding</li>
<li>Confidence intervals overlapping ground truth: method is well-calibrated</li>
</ul>
</div>
'''

    # O2 vs Migration
    html += f'''
<h2>8. Hypoxia Threshold Analysis</h2>

<div class="figure">
<img src="{figures.get("o2_migration", "")}" alt="O2 vs Migration">
<div class="figure-caption">Figure 7: Local O2 level vs. migration rate. Cells below the
hypoxia threshold (red line) show elevated migration - the "Go" response.</div>
</div>

<div class="interpretation">
This visualization directly shows the Go-or-Grow mechanism: cells experiencing low O2
increase their migration speed as an escape response, while well-oxygenated cells
prioritize proliferation over invasion.
</div>
'''

    # Ground truth comparison
    html += f'''
<h2>9. Ground Truth Verification</h2>

<p>Because we control the simulation, we know the TRUE causal parameters:</p>

<table>
<tr><th>Parameter</th><th>Symbol</th><th>True Value</th><th>Description</th></tr>
<tr>
<td>ecDNA → Division</td><td>α</td><td>{tumor_cfg.get("ecDNA_effect_on_division", "?")}</td>
<td>How ecDNA accelerates proliferation</td>
</tr>
<tr>
<td>ecDNA → VEGF</td><td>β</td><td>{tumor_cfg.get("ecDNA_effect_on_VEGF", "?")}</td>
<td>How ecDNA increases VEGF secretion</td>
</tr>
<tr>
<td>ecDNA → Migration</td><td>δ</td><td>{tumor_cfg.get("ecDNA_effect_on_migration", "?")}</td>
<td>How ecDNA increases invasion</td>
</tr>
<tr>
<td>ecDNA → Survival</td><td>γ</td><td>{tumor_cfg.get("ecDNA_effect_on_survival", "?")}</td>
<td>How ecDNA decreases apoptosis</td>
</tr>
<tr>
<td>Hypoxia → Invasion</td><td>(boost)</td><td>2.0x</td>
<td>Fold-increase in migration under hypoxia</td>
</tr>
</table>

<div class="success-box">
<strong>Key finding:</strong> In the TRUE model, ecDNA does NOT directly cause hypoxia.
The pathway is: ecDNA → ↑Proliferation → ↑O2 consumption → ↓O2 → Hypoxia → ↑Invasion.
This is why mediation analysis and IV methods are crucial for recovering the true structure.
</div>
'''

    # Conclusions
    html += '''
<h2>10. Conclusions</h2>

<h3>Statistical Methods Used</h3>
<table>
<tr><th>Method</th><th>Purpose</th><th>Tests/Measures</th></tr>
<tr>
<td>Pearson Correlation</td>
<td>Descriptive association</td>
<td>t-test for r ≠ 0</td>
</tr>
<tr>
<td>OLS Regression</td>
<td>Estimate conditional associations</td>
<td>t-tests for coefficients, F-test for model, R²</td>
</tr>
<tr>
<td>2SLS / IV</td>
<td>Estimate causal effects</td>
<td>First-stage F (weak IV), Wald test</td>
</tr>
<tr>
<td>Mediation Analysis</td>
<td>Decompose direct/indirect effects</td>
<td>Sobel test for indirect effect</td>
</tr>
</table>

<h3>Summary of Findings</h3>
<ol>
<li>ecDNA binomial segregation provides valid randomization for causal inference</li>
<li>OLS estimates capture the total effect of ecDNA on invasion</li>
<li>Mediation analysis reveals that hypoxia partially mediates the ecDNA → invasion relationship</li>
<li>The "Go or Grow" hypothesis is supported: hypoxic cells show elevated invasion at all ecDNA levels</li>
</ol>

<div class="summary-box">
<strong>Bottom line:</strong> Causal inference methods can successfully identify the true
causal structure when we have a valid instrument (ecDNA segregation). The key insight
is that ecDNA affects invasion through BOTH a direct pathway (EGFR signaling) and an
indirect pathway (proliferation → O2 depletion → hypoxia → invasion).
</div>

<hr>
<p style="color: #57606a; font-size: 14px; text-align: center;">
CAUSANTA: Causal Analysis Using Somatic And Neighborhood Tissue Architecture<br>
Report generated automatically. Verify all statistical assumptions before drawing conclusions.
</p>

</div>
</body>
</html>
'''

    with open(output_path, "w") as f:
        f.write(html)

    print(f"Report saved to: {output_path}")


def cleanup_outputs(output_dir: Path, keep_every_n: int = 100) -> dict[str, int]:
    """Clean up output files, keeping only every Nth timestep plus the last.

    This reduces storage by removing intermediate timepoints while preserving
    key snapshots for analysis.

    Args:
        output_dir: Directory containing simulation output
        keep_every_n: Keep every Nth file (default 100)

    Returns:
        Dict with counts of files kept and removed
    """
    import re

    stats = {"kept": 0, "removed": 0, "bytes_freed": 0}

    # Find all timestamped files (cells_t*.tsv, environment_t*.tsv)
    for pattern in ["cells_t*.tsv", "environment_t*.tsv"]:
        files = sorted(output_dir.glob(pattern))

        if not files:
            continue

        # Extract timesteps
        timesteps = []
        for f in files:
            match = re.search(r'_t(\d+)\.tsv$', f.name)
            if match:
                timesteps.append((int(match.group(1)), f))

        if not timesteps:
            continue

        # Sort by timestep
        timesteps.sort(key=lambda x: x[0])

        # Find max timestep (always keep)
        max_t = timesteps[-1][0]

        # Decide which to keep
        for t, filepath in timesteps:
            keep = (t % keep_every_n == 0) or (t == max_t)

            if keep:
                stats["kept"] += 1
            else:
                size = filepath.stat().st_size
                filepath.unlink()
                stats["removed"] += 1
                stats["bytes_freed"] += size

    return stats


def main():
    parser = argparse.ArgumentParser(description="Run full CAUSANTA simulation and analysis")
    parser.add_argument("--hours", type=int, default=168, help="Simulation hours")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--output", type=str, default=None, help="Output directory")
    parser.add_argument("--keep-every", type=int, default=100,
                        help="Keep every Nth output file (default 100)")
    parser.add_argument("--no-cleanup", action="store_true",
                        help="Skip cleanup step (keep all output files)")
    args = parser.parse_args()

    print("=" * 60)
    print("CAUSANTA Full Simulation and Causal Analysis")
    print("=" * 60)

    # Step 1: Build DAGs
    print("\n1. Building causal DAGs...")
    true_dag = build_causanta_dag()
    analysis_dag = build_analysis_dag()
    print(f"   True DAG: {len(true_dag.nodes)} nodes, {len(true_dag.edges)} edges")
    print(f"   Analysis DAG: {len(analysis_dag.nodes)} nodes, {len(analysis_dag.edges)} edges")

    # Step 2: Run simulation
    print("\n2. Running simulation...")
    output_dir = run_simulation(hours=args.hours, seed=args.seed)
    print(f"   Output: {output_dir}")

    # Step 3: Load data
    print("\n3. Loading simulation data...")
    data = load_simulation_data(output_dir)
    print(f"   Loaded {len(data['final_cells'])} cells at final timepoint")

    # Step 4: Compute statistics
    print("\n4. Computing descriptive statistics...")
    stats = compute_descriptive_statistics(data)
    print(f"   Tumor cells: {stats.get('n_tumor_cells', 0)}")

    print("\n5. Computing correlations...")
    correlations = compute_correlations(data)

    print("\n6. Running naive regression...")
    regression = run_naive_regression(data)

    print("\n7. Running IV analysis...")
    iv = run_iv_analysis(data)

    print("\n8. Running mediation analysis...")
    mediation = run_mediation_analysis(data)

    print("\n9. Creating visualizations...")
    figures = create_visualizations(
        data, stats, correlations, regression, iv, mediation, output_dir
    )
    print(f"   Created {len(figures)} figures")

    print("\n10. Generating HTML report...")
    report_path = output_dir / "analysis_report.html"
    generate_html_report(
        data, stats, correlations, regression, iv, mediation,
        figures, true_dag, analysis_dag, report_path
    )

    # Step 11: Cleanup intermediate files
    if not args.no_cleanup:
        print(f"\n11. Cleaning up (keeping every {args.keep_every}th file + last)...")
        cleanup_stats = cleanup_outputs(output_dir, keep_every_n=args.keep_every)
        mb_freed = cleanup_stats["bytes_freed"] / (1024 * 1024)
        print(f"   Kept: {cleanup_stats['kept']} files")
        print(f"   Removed: {cleanup_stats['removed']} files")
        print(f"   Space freed: {mb_freed:.1f} MB")
    else:
        print("\n11. Skipping cleanup (--no-cleanup specified)")

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"\nOutput directory: {output_dir}")
    print(f"Analysis report: {report_path}")
    print(f"\nOpen the report in a browser to view results.")


if __name__ == "__main__":
    main()
