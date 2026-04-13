"""Data loading utilities for CAUSANTA analysis.

Loads simulation output files and prepares data for causal analysis.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator

import numpy as np


@dataclass
class SimulationData:
    """Container for loaded simulation data."""

    output_dir: Path
    config: dict[str, Any]
    cells: dict[int, list[dict]]  # timestep -> list of cell records
    environment: dict[int, list[dict]]  # timestep -> list of env records
    lineage: list[dict]
    summary: list[dict]
    timesteps: list[int]

    @property
    def n_timesteps(self) -> int:
        return len(self.timesteps)

    @property
    def final_cells(self) -> list[dict]:
        """Get cells at final timestep."""
        if not self.timesteps:
            return []
        return self.cells.get(self.timesteps[-1], [])

    @property
    def tumor_cells(self) -> list[dict]:
        """Get tumor cells at final timestep."""
        return [c for c in self.final_cells if c["cell_type"] == 6]

    @property
    def n_divisions(self) -> int:
        return len(self.lineage)

    def get_tumor_ecDNA_series(self) -> dict[int, list[int]]:
        """Get ecDNA counts for tumor cells over time."""
        result = {}
        for t in self.timesteps:
            tumor = [c["ecDNA_count"] for c in self.cells.get(t, []) if c["cell_type"] == 6]
            result[t] = tumor
        return result


def load_config(output_dir: Path) -> dict[str, Any]:
    """Load simulation configuration from params.json."""
    params_path = output_dir / "params.json"
    if not params_path.exists():
        return {}
    with open(params_path) as f:
        return json.load(f)


def load_cells_tsv(filepath: Path) -> list[dict]:
    """Load a cells TSV file."""
    records = []
    if not filepath.exists():
        return records

    with open(filepath, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            records.append({
                "cell_id": int(row["cell_id"]),
                "parent_id": int(row["parent_id"]),
                "cell_type": int(row["cell_type"]),
                "cell_type_name": row["cell_type_name"],
                "x": int(row["x"]),
                "y": int(row["y"]),
                "angle": float(row.get("angle", 0)),
                "ecDNA_count": int(row["ecDNA_count"]),
                "ecDNA_cargo": row.get("ecDNA_cargo", ""),
                "O2_local": float(row["O2_local"]),
                "glucose_local": float(row["glucose_local"]),
                "is_hypoxic": row["is_hypoxic"] == "True",
                "is_quiescent": row.get("is_quiescent", "False") == "True",
                "egfr_expression": float(row.get("egfr_expression", 0)),
                "migration_rate": float(row.get("migration_rate", 0)),
                "VEGF_secretion": float(row.get("VEGF_secretion", 0)),
                "generation": int(row["generation"]),
                "cell_cycle_phase": row.get("cell_cycle_phase", "G0"),
            })
    return records


def load_environment_tsv(filepath: Path) -> list[dict]:
    """Load an environment TSV file."""
    records = []
    if not filepath.exists():
        return records

    with open(filepath, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            records.append({
                "x_grid": int(row["x_grid"]),
                "y_grid": int(row["y_grid"]),
                "x_um": float(row["x_um"]),
                "y_um": float(row["y_um"]),
                "O2": float(row["O2"]),
                "glucose": float(row["glucose"]),
                "VEGF": float(row["VEGF"]),
                "lactate": float(row["lactate"]),
                "pH": float(row.get("pH", 7.4)),
                "vascular_density": float(row.get("vascular_density", 0)),
            })
    return records


def load_lineage_tsv(filepath: Path) -> list[dict]:
    """Load lineage TSV file."""
    records = []
    if not filepath.exists():
        return records

    with open(filepath, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            records.append({
                "time_hr": float(row["time_hr"]),
                "parent_id": int(row["parent_id"]),
                "parent_ecDNA_before": int(row["parent_ecDNA_before"]),
                "parent_ecDNA_after": int(row["parent_ecDNA_after"]),
                "daughter_id": int(row["daughter_id"]),
                "daughter_ecDNA": int(row["daughter_ecDNA"]),
                "x": int(row["x"]),
                "y": int(row["y"]),
                "generation": int(row["generation"]),
            })
    return records


def load_summary_log(filepath: Path) -> list[dict]:
    """Load summary log file."""
    records = []
    if not filepath.exists():
        return records

    with open(filepath, "r") as f:
        lines = f.readlines()
        if not lines:
            return records

        header = lines[0].strip().split("\t")
        for line in lines[1:]:
            values = line.strip().split("\t")
            record = {}
            for i, col in enumerate(header):
                if i < len(values):
                    try:
                        record[col] = float(values[i])
                    except ValueError:
                        record[col] = values[i]
            records.append(record)
    return records


def load_simulation_output(
    output_dir: str | Path,
    load_all_timesteps: bool = False,
    timesteps: list[int] | None = None,
) -> SimulationData:
    """Load complete simulation output.

    Args:
        output_dir: Path to simulation output directory
        load_all_timesteps: If True, load cells/env for all timesteps
        timesteps: Specific timesteps to load (if None and load_all_timesteps=False,
                   loads only first and last)

    Returns:
        SimulationData container with all loaded data.
    """
    output_dir = Path(output_dir)

    # Load config
    config = load_config(output_dir)

    # Find all timestep files
    cells_files = sorted(output_dir.glob("cells_t*.tsv"))
    env_files = sorted(output_dir.glob("environment_t*.tsv"))

    # Extract timesteps from filenames
    all_timesteps = []
    for f in cells_files:
        t = int(f.stem.replace("cells_t", ""))
        all_timesteps.append(t)

    # Determine which timesteps to load
    if timesteps is not None:
        load_timesteps = [t for t in timesteps if t in all_timesteps]
    elif load_all_timesteps:
        load_timesteps = all_timesteps
    else:
        # Load first and last only
        load_timesteps = []
        if all_timesteps:
            load_timesteps.append(all_timesteps[0])
            if len(all_timesteps) > 1:
                load_timesteps.append(all_timesteps[-1])

    # Load cell data
    cells = {}
    for t in load_timesteps:
        cells_path = output_dir / f"cells_t{t:06d}.tsv"
        cells[t] = load_cells_tsv(cells_path)

    # Load environment data
    environment = {}
    for t in load_timesteps:
        env_path = output_dir / f"environment_t{t:06d}.tsv"
        environment[t] = load_environment_tsv(env_path)

    # Load lineage and summary (always load all)
    lineage = load_lineage_tsv(output_dir / "lineage.tsv")
    summary = load_summary_log(output_dir / "summary.log")

    return SimulationData(
        output_dir=output_dir,
        config=config,
        cells=cells,
        environment=environment,
        lineage=lineage,
        summary=summary,
        timesteps=all_timesteps,
    )


def cells_to_array(cells: list[dict], columns: list[str]) -> np.ndarray:
    """Convert cell records to numpy array for analysis.

    Args:
        cells: List of cell dictionaries
        columns: Column names to extract

    Returns:
        2D numpy array with shape (n_cells, n_columns)
    """
    if not cells:
        return np.array([]).reshape(0, len(columns))

    data = []
    for cell in cells:
        row = []
        for col in columns:
            val = cell.get(col, 0)
            if isinstance(val, bool):
                val = 1.0 if val else 0.0
            elif isinstance(val, str):
                val = hash(val) % 1000  # Simple encoding
            row.append(float(val))
        data.append(row)

    return np.array(data)
