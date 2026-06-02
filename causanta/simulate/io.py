"""File I/O for CAUSANTA simulation outputs.

Writes TSV files for cell states, environment fields, lineage records,
and summary statistics as specified in Section 9 of the framework.
"""

from __future__ import annotations

import csv
import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np

from .behaviors import LineageRecord
from .cells import Cell, CellPopulation
from .domain import Domain
from .environment import EnvironmentFields

# TSV column definitions matching spec Section 9.2
CELLS_COLUMNS = [
    "cell_id", "parent_id", "cell_type", "cell_type_name",
    "x", "y", "angle", "shape_path", "size_scale",
    "cell_cycle_phase", "ecDNA_count", "ecDNA_cargo", "egfr_expression",
    "O2_local", "glucose_local", "is_hypoxic", "is_quiescent",
    "generation", "time_born_hr", "migration_rate", "VEGF_secretion",
    # Immune cell state fields
    "activation_level", "exhaustion_level", "kills_performed",
]

# TSV column definitions matching spec Section 9.3
ENVIRONMENT_COLUMNS = [
    "x_grid", "y_grid", "x_um", "y_um",
    "O2", "glucose", "VEGF", "lactate", "pH",
    "ECM_density", "vascular_density",
    "is_occupied", "occupant_cell_id", "occupant_type",
]

# Lineage columns matching spec Section 9.4
LINEAGE_COLUMNS = [
    "time_hr", "parent_id", "parent_ecDNA_before", "parent_ecDNA_after",
    "daughter_id", "daughter_ecDNA", "x", "y", "generation",
]


def create_output_directory(base_path: str | Path = "output", organize: bool = True) -> Path:
    """Create timestamped output directory with organized structure.

    Args:
        base_path: Base output directory
        organize: If True, create subdirectories (data/, figures/, reports/, animations/)

    Returns:
        Path to the run directory
    """
    base = Path(base_path)
    base.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = base / f"run_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=True)

    if organize:
        # Create organized subdirectories
        (run_dir / "data").mkdir(exist_ok=True)
        (run_dir / "figures").mkdir(exist_ok=True)
        (run_dir / "reports").mkdir(exist_ok=True)
        (run_dir / "animations").mkdir(exist_ok=True)

    return run_dir


def write_cells_tsv(population: CellPopulation, filepath: Path) -> None:
    """Write cell state to TSV file."""
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(CELLS_COLUMNS)
        for cell in population.iter_living():
            writer.writerow([
                cell.cell_id,
                cell.parent_id,
                cell.cell_type,
                cell.cell_type_name,
                cell.x,
                cell.y,
                f"{cell.angle:.4f}",
                cell.shape_path,
                f"{cell.size_scale:.4f}",
                cell.cell_cycle_phase,
                cell.ecDNA_count,
                cell.ecDNA_cargo,
                f"{cell.egfr_expression:.4f}",
                f"{cell.O2_local:.2f}",
                f"{cell.glucose_local:.4f}",
                cell.is_hypoxic,
                cell.is_quiescent,
                cell.generation,
                f"{cell.time_born_hr:.1f}",
                f"{cell.migration_rate:.2f}",
                f"{cell.VEGF_secretion:.2f}",
                # Immune cell state
                f"{cell.activation_level:.4f}",
                f"{cell.exhaustion_level:.4f}",
                cell.kills_performed,
            ])


def write_environment_tsv(
    env: EnvironmentFields,
    domain: Domain,
    filepath: Path,
) -> None:
    """Write environment grid to TSV file."""
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(ENVIRONMENT_COLUMNS)
        for gj in range(domain.ny):
            for gi in range(domain.nx):
                x_um, y_um = domain.grid_to_um(gi, gj)
                writer.writerow([
                    gi,
                    gj,
                    f"{x_um:.1f}",
                    f"{y_um:.1f}",
                    f"{env.O2[gj, gi]:.2f}",
                    f"{env.glucose[gj, gi]:.4f}",
                    f"{env.VEGF[gj, gi]:.4f}",
                    f"{env.lactate[gj, gi]:.4f}",
                    f"{env.pH[gj, gi]:.4f}",
                    f"{env.ECM_density[gj, gi]:.4f}",
                    f"{env.vascular_density[gj, gi]:.4f}",
                    env.is_occupied[gj, gi],
                    int(env.occupant_cell_id[gj, gi]),
                    int(env.occupant_type[gj, gi]),
                ])


def write_lineage_tsv(
    records: list[LineageRecord],
    filepath: Path,
    append: bool = False,
) -> None:
    """Write lineage records to TSV file."""
    mode = "a" if append else "w"
    write_header = not append or not filepath.exists()

    with open(filepath, mode, newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        if write_header:
            writer.writerow(LINEAGE_COLUMNS)
        for rec in records:
            writer.writerow([
                f"{rec.time_hr:.1f}",
                rec.parent_id,
                rec.parent_ecDNA_before,
                rec.parent_ecDNA_after,
                rec.daughter_id,
                rec.daughter_ecDNA,
                rec.x,
                rec.y,
                rec.generation,
            ])


def write_summary_log(
    stats: dict,
    filepath: Path,
    step: int,
    append: bool = True,
) -> None:
    """Write per-step summary statistics to log file."""
    mode = "a" if append else "w"
    write_header = not append or not filepath.exists()

    with open(filepath, mode) as f:
        if write_header:
            f.write("step\t" + "\t".join(stats.keys()) + "\n")
        values = "\t".join(
            f"{v:.4f}" if isinstance(v, float) else str(v)
            for v in stats.values()
        )
        f.write(f"{step}\t{values}\n")


def copy_params_json(config_path: str | Path, output_dir: Path) -> None:
    """Copy parameter JSON file to output directory."""
    src = Path(config_path)
    if src.exists():
        # Use data subdirectory if it exists
        data_dir = output_dir / "data"
        dest_dir = data_dir if data_dir.exists() else output_dir
        shutil.copy2(src, dest_dir / "params.json")
