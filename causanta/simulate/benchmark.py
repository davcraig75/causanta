#!/usr/bin/env python3
"""Benchmark script for CAUSANTA simulation performance.

Profiles individual components to identify bottlenecks.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from .config import load_config
from .core import Simulation
from .domain import Domain
from .environment import (
    DiffusionSolver,
    EnvironmentFields,
    _thomas_solve_batch,
    accumulate_sources_and_sinks,
)


def benchmark_thomas_algorithm(n_iterations: int = 100) -> dict:
    """Benchmark the Thomas algorithm solver."""
    # Create test data matching typical simulation grid
    ny, nx = 100, 100  # Typical grid size
    n_interior = nx - 2

    # Create batch of tridiagonal systems
    r = 0.1  # Typical diffusion parameter
    a = np.full((ny, n_interior - 1), -r)
    b = np.full((ny, n_interior), 1.0 + 2.0 * r)
    c = np.full((ny, n_interior - 1), -r)
    d = np.random.rand(ny, n_interior)

    # Warmup
    for _ in range(5):
        _thomas_solve_batch(a, b, c, d)

    # Benchmark
    start = time.perf_counter()
    for _ in range(n_iterations):
        _thomas_solve_batch(a, b, c, d)
    elapsed = time.perf_counter() - start

    return {
        "function": "_thomas_solve_batch",
        "grid_size": f"{ny}x{nx}",
        "iterations": n_iterations,
        "total_time_ms": elapsed * 1000,
        "time_per_call_ms": elapsed * 1000 / n_iterations,
    }


def benchmark_diffusion_step(n_iterations: int = 10) -> dict:
    """Benchmark full diffusion solver step."""
    # Create minimal config
    config_path = Path(__file__).parent / "params" / "default.json"
    if not config_path.exists():
        # Create a minimal test
        return {"error": "No default config found"}

    config = load_config(config_path)
    domain = Domain(config.domain)
    env = EnvironmentFields(domain, config.environment)
    solver = DiffusionSolver(domain, config.environment)

    # Create dummy sources/sinks
    ny, nx = domain.ny, domain.nx
    sources = {
        "O2": np.random.rand(ny, nx) * 0.01,
        "glucose": np.random.rand(ny, nx) * 0.01,
        "VEGF": np.random.rand(ny, nx) * 0.001,
        "lactate": np.random.rand(ny, nx) * 0.01,
    }
    sinks = {
        "O2": np.random.rand(ny, nx) * 0.01,
        "glucose": np.random.rand(ny, nx) * 0.01,
        "VEGF": np.zeros((ny, nx)),
        "lactate": np.zeros((ny, nx)),
    }

    dt = config.time.dt_diffusion_hr
    n_substeps = max(1, int(1.0 / dt))

    # Warmup
    for _ in range(2):
        solver.step(env, sources, sinks, dt, n_substeps)

    # Benchmark
    start = time.perf_counter()
    for _ in range(n_iterations):
        solver.step(env, sources, sinks, dt, n_substeps)
    elapsed = time.perf_counter() - start

    return {
        "function": "DiffusionSolver.step",
        "grid_size": f"{domain.ny}x{domain.nx}",
        "substeps_per_call": n_substeps,
        "iterations": n_iterations,
        "total_time_ms": elapsed * 1000,
        "time_per_call_ms": elapsed * 1000 / n_iterations,
    }


def benchmark_full_simulation(n_steps: int = 10) -> dict:
    """Benchmark a full simulation for n steps."""
    config_path = Path(__file__).parent / "params" / "default.json"
    if not config_path.exists():
        return {"error": "No default config found"}

    import tempfile
    import json

    # Create config with limited steps
    with open(config_path) as f:
        data = json.load(f)
    data["time"]["total_hours"] = n_steps
    data["time"]["output_interval_hr"] = n_steps + 1  # No output

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
        json.dump(data, f)
        tmp_path = f.name

    try:
        sim = Simulation(tmp_path)

        # Initialize (don't time this)
        sim.initialize()

        # Time the stepping
        phase_times = {
            "environment": 0.0,
            "cells": 0.0,
            "angiogenesis": 0.0,
            "immune": 0.0,
            "cleanup": 0.0,
        }

        total_start = time.perf_counter()

        for step in range(1, n_steps + 1):
            sim.current_time_hr = float(step)

            t0 = time.perf_counter()
            sim._step_environment()
            phase_times["environment"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            sim._step_cells()
            phase_times["cells"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            sim._step_angiogenesis()
            phase_times["angiogenesis"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            sim._step_immune_recruitment()
            phase_times["immune"] += time.perf_counter() - t0

            t0 = time.perf_counter()
            sim._step_cleanup()
            phase_times["cleanup"] += time.perf_counter() - t0

        total_elapsed = time.perf_counter() - total_start

        # Convert to percentages
        for key in phase_times:
            phase_times[key] = {
                "time_ms": phase_times[key] * 1000,
                "percent": phase_times[key] / total_elapsed * 100,
            }

        return {
            "n_steps": n_steps,
            "total_time_ms": total_elapsed * 1000,
            "time_per_step_ms": total_elapsed * 1000 / n_steps,
            "final_cell_count": sim.population.count(),
            "phases": phase_times,
        }
    finally:
        import os
        os.unlink(tmp_path)


def run_benchmarks():
    """Run all benchmarks and print results."""
    print("=" * 60)
    print("CAUSANTA Simulation Performance Benchmark")
    print("=" * 60)
    print()

    # Thomas algorithm
    print("1. Thomas Algorithm (tridiagonal solver)")
    print("-" * 40)
    result = benchmark_thomas_algorithm(n_iterations=100)
    for k, v in result.items():
        print(f"   {k}: {v}")
    print()

    # Diffusion solver
    print("2. Diffusion Solver (full step)")
    print("-" * 40)
    result = benchmark_diffusion_step(n_iterations=10)
    for k, v in result.items():
        print(f"   {k}: {v}")
    print()

    # Full simulation
    print("3. Full Simulation (per-phase breakdown)")
    print("-" * 40)
    result = benchmark_full_simulation(n_steps=20)
    if "error" not in result:
        print(f"   Total time: {result['total_time_ms']:.1f} ms for {result['n_steps']} steps")
        print(f"   Time per step: {result['time_per_step_ms']:.1f} ms")
        print(f"   Final cells: {result['final_cell_count']}")
        print()
        print("   Phase breakdown:")
        for phase, data in result["phases"].items():
            print(f"      {phase:15s}: {data['time_ms']:8.1f} ms ({data['percent']:5.1f}%)")
    else:
        print(f"   {result['error']}")
    print()
    print("=" * 60)


if __name__ == "__main__":
    run_benchmarks()
