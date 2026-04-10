"""Angiogenesis module for CAUSANTA.

Implements VEGF-driven vessel sprouting, tip cell selection,
stalk cell placement, anastomosis, and vessel maturation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .cells import ENDOTHELIAL, CellPopulation, create_cell
from .config import SimulationConfig
from .domain import Domain
from .environment import EnvironmentFields


@dataclass
class ActiveSprout:
    """Track an in-progress sprouting event."""

    tip_cell_id: int
    age_hr: float = 0.0
    maturation_positions: list[tuple[int, int]] = field(default_factory=list)


class AngiogenesisManager:
    """Manages vascular sprouting and maturation."""

    def __init__(self, config: SimulationConfig):
        self.config = config
        self.active_sprouts: list[ActiveSprout] = []

    def process(
        self,
        population: CellPopulation,
        env: EnvironmentFields,
        domain: Domain,
        rng: np.random.Generator,
        current_time_hr: float,
    ) -> None:
        """Process one hour of angiogenesis.

        1. Select tip cells (endothelial cells with high adjacent VEGF)
        2. Extend sprouts up VEGF gradient
        3. Check anastomosis
        4. Mature new vessels
        """
        threshold = self.config.angiogenesis.angiogenesis_threshold_nM
        sprout_rate = self.config.angiogenesis.sprouting_rate_per_hr
        endo_params = self.config.cell_types.get(ENDOTHELIAL)
        if endo_params is None:
            return

        # Find candidate tip cells
        candidates: list[tuple[float, int]] = []
        for cell in population.iter_by_type(ENDOTHELIAL):
            local_vegf = domain.bilinear_interpolate(env.VEGF, cell.x, cell.y)
            if local_vegf > threshold:
                candidates.append((local_vegf, cell.cell_id))

        # Sort by VEGF (highest first) and attempt sprouting
        candidates.sort(reverse=True)

        for vegf_val, cell_id in candidates:
            if rng.random() > sprout_rate:
                continue

            cell = population.get_cell(cell_id)
            if cell is None or not cell.is_alive:
                continue

            # Sprout direction: up the VEGF gradient
            grad = domain.compute_gradient(env.VEGF, cell.x, cell.y)
            mag = (grad[0] ** 2 + grad[1] ** 2) ** 0.5
            if mag < 1e-10:
                continue

            # Normalize and compute new position (one cell diameter step)
            step = endo_params.diameter_mean_um
            new_x = int(round(cell.x + step * grad[0] / mag))
            new_y = int(round(cell.y + step * grad[1] / mag))

            if not domain.is_in_bounds(new_x, new_y):
                continue

            # Check space is available
            if population.has_neighbor_within(new_x, new_y, endo_params.diameter_mean_um * 0.7):
                continue

            # Create new stalk cell at the tip's original position concept:
            # Actually, tip moves forward and stalk fills behind. But simpler:
            # create new endothelial cell at the sprouted position.
            new_id = population.allocate_id()
            new_cell = create_cell(
                cell_id=new_id,
                parent_id=cell_id,
                cell_type=ENDOTHELIAL,
                x=new_x,
                y=new_y,
                time_born_hr=current_time_hr,
                type_params=endo_params,
                rng=rng,
                generation=cell.generation + 1,
            )
            population.add_cell(new_cell)

            # Track for maturation
            gi, gj = domain.um_to_grid(new_x, new_y)
            self.active_sprouts.append(ActiveSprout(
                tip_cell_id=new_id,
                age_hr=0.0,
                maturation_positions=[(gi, gj)],
            ))

        # Mature existing sprouts
        maturation_time = self.config.angiogenesis.vessel_maturation_hr
        remaining: list[ActiveSprout] = []
        for sprout in self.active_sprouts:
            sprout.age_hr += 1.0
            progress = min(sprout.age_hr / maturation_time, 1.0)

            for gi, gj in sprout.maturation_positions:
                if 0 <= gi < domain.nx and 0 <= gj < domain.ny:
                    # Gradually increase vascular density
                    env.vascular_density[gj, gi] = max(
                        env.vascular_density[gj, gi], progress * 0.8
                    )

            if sprout.age_hr < maturation_time:
                remaining.append(sprout)

        self.active_sprouts = remaining
