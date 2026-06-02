"""Parameter sweep system for systematic simulation studies.

Enables:
- Systematic variation of causal effect parameters (α, β, δ, γ)
- Sample size studies (vary simulation duration)
- Domain size studies
- Confounding strength studies
- Multi-replicate runs for statistical power

Designed for reproducible, parallel execution of simulation batteries.
"""

from __future__ import annotations

import itertools
import json
import multiprocessing
import os
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class ParameterSet:
    """A single parameter configuration for simulation."""

    name: str
    params: dict[str, Any]
    seed: int
    replicate: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "params": self.params,
            "seed": self.seed,
            "replicate": self.replicate,
        }


@dataclass
class SweepConfig:
    """Configuration for a parameter sweep study."""

    name: str
    base_config_path: Path
    output_dir: Path
    n_replicates: int = 10
    base_seed: int = 42

    # Parameter variations (each key maps to list of values to try)
    vary_params: dict[str, list[Any]] = field(default_factory=dict)

    # Fixed overrides (applied to all runs)
    fixed_params: dict[str, Any] = field(default_factory=dict)

    def generate_parameter_sets(self) -> list[ParameterSet]:
        """Generate all parameter combinations."""
        param_sets = []

        if not self.vary_params:
            # No variation - just replicates
            for rep in range(self.n_replicates):
                seed = self.base_seed + rep
                param_sets.append(ParameterSet(
                    name=f"baseline_rep{rep:02d}",
                    params=dict(self.fixed_params),
                    seed=seed,
                    replicate=rep,
                ))
            return param_sets

        # Generate all combinations
        param_names = list(self.vary_params.keys())
        param_values = list(self.vary_params.values())

        for combo in itertools.product(*param_values):
            combo_dict = dict(zip(param_names, combo))
            combo_name = "_".join(f"{k}{v}" for k, v in combo_dict.items())

            for rep in range(self.n_replicates):
                seed = self.base_seed + rep
                params = {**self.fixed_params, **combo_dict}

                param_sets.append(ParameterSet(
                    name=f"{combo_name}_rep{rep:02d}",
                    params=params,
                    seed=seed,
                    replicate=rep,
                ))

        return param_sets

    def total_runs(self) -> int:
        """Total number of simulation runs."""
        if not self.vary_params:
            return self.n_replicates

        n_combos = 1
        for values in self.vary_params.values():
            n_combos *= len(values)
        return n_combos * self.n_replicates


def apply_params_to_config(
    base_config: dict,
    params: dict[str, Any],
) -> dict:
    """Apply parameter overrides to base configuration.

    Supports nested paths like "cell_types.6.ecDNA_effect_on_division"
    for modifying specific values.
    """
    config = json.loads(json.dumps(base_config))  # Deep copy

    for key, value in params.items():
        parts = key.split(".")
        obj = config

        # Navigate to parent
        for part in parts[:-1]:
            if part.isdigit():
                part = int(part)
            if isinstance(obj, dict):
                obj = obj.get(part, {})
            elif isinstance(obj, list):
                obj = obj[part]

        # Set value
        final_key = parts[-1]
        if final_key.isdigit():
            final_key = int(final_key)

        if isinstance(obj, dict):
            obj[final_key] = value

    return config


def create_sweep_config_file(
    base_config_path: Path,
    params: dict[str, Any],
    seed: int,
    output_path: Path,
) -> Path:
    """Create a modified config file for a sweep run."""
    with open(base_config_path) as f:
        config = json.load(f)

    # Apply parameter overrides
    config = apply_params_to_config(config, params)

    # Set seed
    config["rng_seed"] = seed

    # Write to output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(config, f, indent=2)

    return output_path


def run_single_simulation(args: tuple) -> dict[str, Any]:
    """Run a single simulation (for parallel execution).

    Args:
        args: Tuple of (param_set, base_config_path, output_base_dir)

    Returns:
        Dictionary with run results
    """
    param_set, base_config_path, output_base_dir = args

    # Import here to avoid circular imports
    from .core import Simulation

    # Create output directory for this run
    run_dir = Path(output_base_dir) / param_set.name
    run_dir.mkdir(parents=True, exist_ok=True)

    # Create config file
    config_path = run_dir / "params.json"
    create_sweep_config_file(
        base_config_path,
        param_set.params,
        param_set.seed,
        config_path,
    )

    # Run simulation
    try:
        sim = Simulation(config_path)
        sim.output_dir = run_dir
        sim.run()

        return {
            "name": param_set.name,
            "status": "success",
            "output_dir": str(run_dir),
            "params": param_set.params,
            "seed": param_set.seed,
        }
    except Exception as e:
        return {
            "name": param_set.name,
            "status": "error",
            "error": str(e),
            "params": param_set.params,
            "seed": param_set.seed,
        }


def run_sweep(
    sweep_config: SweepConfig,
    n_workers: int | None = None,
    progress_callback: Any = None,
) -> list[dict[str, Any]]:
    """Run a complete parameter sweep.

    Args:
        sweep_config: Sweep configuration
        n_workers: Number of parallel workers (None = n_cpus - 1)
        progress_callback: Optional callback for progress updates

    Returns:
        List of run results
    """
    if n_workers is None:
        n_workers = max(1, multiprocessing.cpu_count() - 1)

    # Generate parameter sets
    param_sets = sweep_config.generate_parameter_sets()

    # Create output directory
    sweep_config.output_dir.mkdir(parents=True, exist_ok=True)

    # Save sweep configuration
    sweep_meta = {
        "name": sweep_config.name,
        "n_replicates": sweep_config.n_replicates,
        "base_seed": sweep_config.base_seed,
        "vary_params": sweep_config.vary_params,
        "fixed_params": sweep_config.fixed_params,
        "total_runs": len(param_sets),
    }
    with open(sweep_config.output_dir / "sweep_config.json", "w") as f:
        json.dump(sweep_meta, f, indent=2)

    # Prepare arguments
    args = [
        (ps, sweep_config.base_config_path, sweep_config.output_dir)
        for ps in param_sets
    ]

    # Run simulations
    results = []

    if n_workers == 1:
        # Sequential execution
        for i, arg in enumerate(args):
            result = run_single_simulation(arg)
            results.append(result)
            if progress_callback:
                progress_callback(i + 1, len(args), result)
    else:
        # Parallel execution
        with multiprocessing.Pool(n_workers) as pool:
            for i, result in enumerate(pool.imap_unordered(run_single_simulation, args)):
                results.append(result)
                if progress_callback:
                    progress_callback(i + 1, len(args), result)

    # Save results summary
    with open(sweep_config.output_dir / "sweep_results.json", "w") as f:
        json.dump(results, f, indent=2)

    return results


# Pre-defined sweep configurations for common studies

def baseline_validation_sweep(
    base_config: Path,
    output_dir: Path,
    n_replicates: int = 20,
) -> SweepConfig:
    """Baseline validation: multiple replicates with default parameters."""
    return SweepConfig(
        name="baseline_validation",
        base_config_path=base_config,
        output_dir=output_dir,
        n_replicates=n_replicates,
        vary_params={},
        fixed_params={},
    )


def effect_size_sweep(
    base_config: Path,
    output_dir: Path,
    n_replicates: int = 10,
) -> SweepConfig:
    """Vary all four causal effect parameters."""
    return SweepConfig(
        name="effect_size_sweep",
        base_config_path=base_config,
        output_dir=output_dir,
        n_replicates=n_replicates,
        vary_params={
            "cell_types.6.ecDNA_effect_on_division": [0.15, 0.3, 0.45],
            "cell_types.6.ecDNA_effect_on_VEGF": [0.05, 0.1, 0.15],
            "cell_types.6.ecDNA_effect_on_migration": [0.025, 0.05, 0.075],
            "cell_types.6.ecDNA_effect_on_survival": [0.1, 0.2, 0.3],
        },
    )


def sample_size_sweep(
    base_config: Path,
    output_dir: Path,
    n_replicates: int = 10,
) -> SweepConfig:
    """Vary simulation duration to study sample size effects."""
    return SweepConfig(
        name="sample_size_sweep",
        base_config_path=base_config,
        output_dir=output_dir,
        n_replicates=n_replicates,
        vary_params={
            "time.total_hours": [50, 100, 150, 200, 300, 500],
        },
    )


def domain_size_sweep(
    base_config: Path,
    output_dir: Path,
    n_replicates: int = 10,
) -> SweepConfig:
    """Vary domain size to study scaling behavior."""
    return SweepConfig(
        name="domain_size_sweep",
        base_config_path=base_config,
        output_dir=output_dir,
        n_replicates=n_replicates,
        vary_params={
            "domain.width_um": [500, 1000, 2000],
            "domain.height_um": [500, 1000, 2000],
        },
    )


def initial_ecdna_sweep(
    base_config: Path,
    output_dir: Path,
    n_replicates: int = 10,
) -> SweepConfig:
    """Vary initial ecDNA copy number."""
    return SweepConfig(
        name="initial_ecdna_sweep",
        base_config_path=base_config,
        output_dir=output_dir,
        n_replicates=n_replicates,
        vary_params={
            "tumor_seeds.0.ecDNA_count": [5, 10, 20, 40, 80],
        },
    )


def confounding_strength_sweep(
    base_config: Path,
    output_dir: Path,
    n_replicates: int = 10,
) -> SweepConfig:
    """Vary hypoxia effect (confounding strength)."""
    return SweepConfig(
        name="confounding_sweep",
        base_config_path=base_config,
        output_dir=output_dir,
        n_replicates=n_replicates,
        vary_params={
            # Vary the hypoxia threshold to change confounding severity
            "environment.hypoxia_threshold_mmHg": [5.0, 10.0, 15.0, 20.0],
        },
    )


def aggregate_sweep_results(
    sweep_dir: Path,
) -> dict[str, Any]:
    """Aggregate results from a completed sweep.

    Loads all results and computes summary statistics across
    replicates and parameter combinations.

    Args:
        sweep_dir: Directory containing sweep outputs

    Returns:
        Aggregated results dictionary
    """
    # Load sweep config
    config_path = sweep_dir / "sweep_config.json"
    if not config_path.exists():
        return {"error": "No sweep_config.json found"}

    with open(config_path) as f:
        sweep_config = json.load(f)

    # Load results
    results_path = sweep_dir / "sweep_results.json"
    if not results_path.exists():
        return {"error": "No sweep_results.json found"}

    with open(results_path) as f:
        results = json.load(f)

    # Group by parameter combination
    from collections import defaultdict
    grouped = defaultdict(list)

    for result in results:
        if result["status"] != "success":
            continue

        # Create key from params (excluding replicate-varying items)
        params = result.get("params", {})
        key = tuple(sorted(params.items()))
        grouped[key].append(result)

    # Aggregate
    aggregated = {
        "sweep_config": sweep_config,
        "total_runs": len(results),
        "successful_runs": sum(1 for r in results if r["status"] == "success"),
        "failed_runs": sum(1 for r in results if r["status"] != "success"),
        "parameter_groups": [],
    }

    for params_key, group_results in grouped.items():
        params_dict = dict(params_key)
        group_summary = {
            "params": params_dict,
            "n_replicates": len(group_results),
            "output_dirs": [r["output_dir"] for r in group_results],
        }
        aggregated["parameter_groups"].append(group_summary)

    return aggregated


def collect_analysis_from_sweep(
    sweep_dir: Path,
    analysis_fn: Any,
) -> list[dict[str, Any]]:
    """Run an analysis function on all successful sweep runs.

    Args:
        sweep_dir: Directory containing sweep outputs
        analysis_fn: Function that takes output_dir and returns results

    Returns:
        List of analysis results with params attached
    """
    results_path = sweep_dir / "sweep_results.json"
    if not results_path.exists():
        return []

    with open(results_path) as f:
        results = json.load(f)

    analysis_results = []

    for result in results:
        if result["status"] != "success":
            continue

        output_dir = Path(result["output_dir"])
        try:
            analysis = analysis_fn(output_dir)
            analysis_results.append({
                "name": result["name"],
                "params": result.get("params", {}),
                "seed": result.get("seed"),
                "analysis": analysis,
            })
        except Exception as e:
            analysis_results.append({
                "name": result["name"],
                "params": result.get("params", {}),
                "error": str(e),
            })

    return analysis_results
