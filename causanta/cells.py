"""Cell agent model and population management for CAUSANTA."""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Iterator

import numpy as np

from .config import CellTypeConfig, DEFAULT_SHAPE_PATHS
from .domain import Domain

# Cell type ID constants
NEURON = 0
ASTROCYTE = 1
OLIGODENDROCYTE = 2
MICROGLIA = 3
ENDOTHELIAL = 4
PERICYTE = 5
TUMOR = 6
RECRUITED_IMMUNE = 7
NECROTIC = 8

TYPE_NAMES = {
    NEURON: "Neuron",
    ASTROCYTE: "Astrocyte",
    OLIGODENDROCYTE: "Oligodendrocyte",
    MICROGLIA: "Microglia",
    ENDOTHELIAL: "Endothelial",
    PERICYTE: "Pericyte",
    TUMOR: "Tumor",
    RECRUITED_IMMUNE: "RecruitedImmune",
    NECROTIC: "Necrotic",
}


@dataclass
class Cell:
    """Complete state vector for a single cell agent."""

    cell_id: int
    parent_id: int
    cell_type: int
    cell_type_name: str
    # Spatial
    x: int
    y: int
    angle: float = 0.0
    # Morphology
    shape_path: str = ""
    size_scale: float = 1.0
    # Proliferation
    cell_cycle_phase: str = "G0"
    cycle_clock_hr: float = 0.0
    total_cycle_time_hr: float = 0.0
    # Genomic / ecDNA
    ecDNA_count: int = 0
    ecDNA_cargo: str = ""
    # Lineage
    generation: int = 0
    time_born_hr: float = 0.0
    # Local environment (sampled each step)
    O2_local: float = 38.0
    glucose_local: float = 5.0
    # Derived state
    is_hypoxic: bool = False
    is_quiescent: bool = False
    # Current effective rates
    migration_rate: float = 0.0
    VEGF_secretion: float = 0.0
    # Internal bookkeeping (not written to TSV)
    necrosis_timer_hr: float = 0.0
    persist_direction: float = 0.0
    is_alive: bool = True
    is_reactive: bool = False
    activation_level: float = 0.0
    diameter: float = 10.0


def create_cell(
    cell_id: int,
    parent_id: int,
    cell_type: int,
    x: int,
    y: int,
    time_born_hr: float,
    type_params: CellTypeConfig,
    rng: np.random.Generator,
    *,
    ecDNA_count: int = 0,
    ecDNA_cargo: str = "",
    generation: int = 0,
    cell_cycle_phase: str = "G0",
) -> Cell:
    """Factory function to create a properly initialized cell."""
    diameter = max(1.0, rng.normal(type_params.diameter_mean_um, type_params.diameter_std_um))
    size_scale = diameter / type_params.diameter_mean_um if type_params.diameter_mean_um > 0 else 1.0
    type_name = TYPE_NAMES.get(cell_type, "Unknown")

    # Get shape path
    shape_path = type_params.shape_path
    if not shape_path:
        shape_path = DEFAULT_SHAPE_PATHS.get(type_name, "")

    # Initialize cell cycle if the cell can divide and is not in G0
    total_cycle_time = 0.0
    if type_params.can_divide and cell_cycle_phase != "G0":
        total_cycle_time = max(1.0, rng.normal(type_params.division_time_mean_hr, type_params.division_time_std_hr))

    return Cell(
        cell_id=cell_id,
        parent_id=parent_id,
        cell_type=cell_type,
        cell_type_name=type_name,
        x=x,
        y=y,
        angle=rng.uniform(0, 2 * math.pi),
        shape_path=shape_path,
        size_scale=size_scale,
        cell_cycle_phase=cell_cycle_phase,
        cycle_clock_hr=0.0,
        total_cycle_time_hr=total_cycle_time,
        ecDNA_count=ecDNA_count,
        ecDNA_cargo=ecDNA_cargo,
        generation=generation,
        time_born_hr=time_born_hr,
        migration_rate=type_params.migration_speed_um_hr,
        VEGF_secretion=type_params.VEGF_secretion_amol_hr,
        persist_direction=rng.uniform(0, 2 * math.pi),
        diameter=diameter,
    )


class CellPopulation:
    """Manages all cell agents with spatial hash indexing for fast neighbor queries."""

    SPATIAL_HASH_RESOLUTION = 20  # um per bucket

    def __init__(self, domain: Domain):
        self._cells: dict[int, Cell] = {}
        self._next_id: int = 0
        self._domain = domain
        self._spatial_grid: dict[tuple[int, int], set[int]] = defaultdict(set)

    @property
    def next_id(self) -> int:
        return self._next_id

    def _hash_key(self, x: int, y: int) -> tuple[int, int]:
        return x // self.SPATIAL_HASH_RESOLUTION, y // self.SPATIAL_HASH_RESOLUTION

    def allocate_id(self) -> int:
        """Allocate and return the next cell ID."""
        cid = self._next_id
        self._next_id += 1
        return cid

    def add_cell(self, cell: Cell) -> None:
        """Add cell to population and spatial index."""
        self._cells[cell.cell_id] = cell
        key = self._hash_key(cell.x, cell.y)
        self._spatial_grid[key].add(cell.cell_id)
        if cell.cell_id >= self._next_id:
            self._next_id = cell.cell_id + 1

    def remove_cell(self, cell_id: int) -> None:
        """Remove cell from population and spatial index."""
        cell = self._cells.get(cell_id)
        if cell is None:
            return
        key = self._hash_key(cell.x, cell.y)
        self._spatial_grid[key].discard(cell_id)
        del self._cells[cell_id]

    def get_cell(self, cell_id: int) -> Cell | None:
        return self._cells.get(cell_id)

    def update_position(self, cell: Cell, old_x: int, old_y: int) -> None:
        """Update spatial index after cell movement."""
        old_key = self._hash_key(old_x, old_y)
        new_key = self._hash_key(cell.x, cell.y)
        if old_key != new_key:
            self._spatial_grid[old_key].discard(cell.cell_id)
            self._spatial_grid[new_key].add(cell.cell_id)

    def get_neighbors(self, x: int, y: int, radius: float) -> list[Cell]:
        """Get all living cells within radius of (x, y)."""
        r_sq = radius * radius
        bucket_range = int(math.ceil(radius / self.SPATIAL_HASH_RESOLUTION)) + 1
        cx, cy = x // self.SPATIAL_HASH_RESOLUTION, y // self.SPATIAL_HASH_RESOLUTION
        result: list[Cell] = []
        for bx in range(cx - bucket_range, cx + bucket_range + 1):
            for by in range(cy - bucket_range, cy + bucket_range + 1):
                for cid in self._spatial_grid.get((bx, by), ()):
                    cell = self._cells.get(cid)
                    if cell is not None and cell.is_alive:
                        dx = cell.x - x
                        dy = cell.y - y
                        if dx * dx + dy * dy <= r_sq:
                            result.append(cell)
        return result

    def has_neighbor_within(self, x: int, y: int, radius: float, exclude_id: int = -1) -> bool:
        """Check if any living cell is within radius of (x, y)."""
        r_sq = radius * radius
        bucket_range = int(math.ceil(radius / self.SPATIAL_HASH_RESOLUTION)) + 1
        cx, cy = x // self.SPATIAL_HASH_RESOLUTION, y // self.SPATIAL_HASH_RESOLUTION
        for bx in range(cx - bucket_range, cx + bucket_range + 1):
            for by in range(cy - bucket_range, cy + bucket_range + 1):
                for cid in self._spatial_grid.get((bx, by), ()):
                    if cid == exclude_id:
                        continue
                    cell = self._cells.get(cid)
                    if cell is not None and cell.is_alive:
                        dx = cell.x - x
                        dy = cell.y - y
                        if dx * dx + dy * dy <= r_sq:
                            return True
        return False

    def get_adjacent_empty_positions(
        self, x: int, y: int, cell_diameter: float, rng: np.random.Generator
    ) -> list[tuple[int, int]]:
        """Find unoccupied positions in 8-connected neighborhood for division.

        Checks positions at distance = cell_diameter in 8 cardinal/diagonal directions.
        Returns shuffled list of positions that are in-bounds and unoccupied.
        """
        step = int(round(cell_diameter))
        offsets = [
            (-step, -step), (0, -step), (step, -step),
            (-step, 0),                  (step, 0),
            (-step, step),  (0, step),   (step, step),
        ]
        candidates = []
        for dx, dy in offsets:
            nx, ny = x + dx, y + dy
            if not self._domain.is_in_bounds(nx, ny):
                continue
            # Check for collisions at the candidate position
            if not self.has_neighbor_within(nx, ny, cell_diameter * 0.8):
                candidates.append((nx, ny))
        rng.shuffle(candidates)
        return candidates

    def find_division_position_with_displacement(
        self,
        x: int,
        y: int,
        cell_diameter: float,
        parent_cell_id: int,
        rng: np.random.Generator,
    ) -> tuple[int, int] | None:
        """Find a position for daughter cell, displacing neighbors if necessary.

        For non-contact-inhibited cells (e.g., tumor), if no empty adjacent
        position exists, push the least resistant neighbor aside.
        Returns position or None if impossible.
        """
        # First try normal empty positions
        positions = self.get_adjacent_empty_positions(x, y, cell_diameter, rng)
        if positions:
            return positions[0]

        # No empty position — try to displace a neighbor
        step = int(round(cell_diameter))
        offsets = [
            (-step, -step), (0, -step), (step, -step),
            (-step, 0),                  (step, 0),
            (-step, step),  (0, step),   (step, step),
        ]
        rng.shuffle(offsets)

        for dx, dy in offsets:
            nx, ny = x + dx, y + dy
            if not self._domain.is_in_bounds(nx, ny):
                continue

            # Find occupants at this position
            neighbors = self.get_neighbors(nx, ny, cell_diameter * 0.8)
            # Filter to displaceable cells (not the parent, not endothelial/vascular)
            displaceable = [
                n for n in neighbors
                if n.cell_id != parent_cell_id
                and n.cell_type not in (ENDOTHELIAL, PERICYTE)
                and n.is_alive
            ]
            if not displaceable:
                continue

            # Displace the first displaceable neighbor by pushing it outward
            target = displaceable[0]
            push_dx = target.x - x
            push_dy = target.y - y
            mag = max(1.0, (push_dx**2 + push_dy**2) ** 0.5)
            push_dist = int(round(cell_diameter))
            new_tx = int(round(target.x + push_dist * push_dx / mag))
            new_ty = int(round(target.y + push_dist * push_dy / mag))
            new_tx, new_ty = self._domain.clamp_position(new_tx, new_ty)

            # Move the displaced cell
            old_tx, old_ty = target.x, target.y
            target.x = new_tx
            target.y = new_ty
            self.update_position(target, old_tx, old_ty)

            return (nx, ny)

    def iter_living(self) -> Iterator[Cell]:
        """Iterate over all living cells."""
        for cell in self._cells.values():
            if cell.is_alive:
                yield cell

    def iter_by_type(self, type_id: int) -> Iterator[Cell]:
        """Iterate over living cells of a specific type."""
        for cell in self._cells.values():
            if cell.is_alive and cell.cell_type == type_id:
                yield cell

    def count(self) -> int:
        """Total number of living cells."""
        return sum(1 for c in self._cells.values() if c.is_alive)

    def count_by_type(self) -> dict[int, int]:
        """Count living cells per type."""
        counts: dict[int, int] = defaultdict(int)
        for cell in self._cells.values():
            if cell.is_alive:
                counts[cell.cell_type] += 1
        return dict(counts)

    def shuffled_living(self, rng: np.random.Generator) -> list[Cell]:
        """Return living cells in random order for fair update scheduling."""
        cells = [c for c in self._cells.values() if c.is_alive]
        rng.shuffle(cells)
        return cells

    def all_cells(self) -> Iterator[Cell]:
        """Iterate over ALL cells (including dead/removed for output purposes)."""
        yield from self._cells.values()
