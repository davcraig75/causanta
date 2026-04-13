"""Cell behavioral rules for CAUSANTA.

Implements proliferation (with cell cycle), migration (with chemotaxis),
death (necrosis, apoptosis, immune kill), and state transitions.
All rates are per hour (the simulation time unit).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np

from .cells import (
    ASTROCYTE,
    ENDOTHELIAL,
    MICROGLIA,
    NECROTIC,
    PERICYTE,
    TUMOR,
    Cell,
    CellPopulation,
    create_cell,
)
from .config import CellTypeConfig, SimulationConfig
from .domain import Domain
from .ecdna import (
    modulate_apoptosis_rate,
    modulate_division_time,
    modulate_migration_speed,
    modulate_vegf_secretion,
    segregate_ecdna,
)
from .environment import EnvironmentFields

# Cell cycle phase fractions of total division time
PHASE_FRACTIONS = {"G1": 0.40, "S": 0.30, "G2": 0.15, "M": 0.15}
PHASE_ORDER = ["G1", "S", "G2", "M"]


@dataclass
class LineageRecord:
    """Record of a single division event."""

    time_hr: float
    parent_id: int
    parent_ecDNA_before: int
    parent_ecDNA_after: int
    daughter_id: int
    daughter_ecDNA: int
    x: int
    y: int
    generation: int


def _get_phase_boundary(total_time: float, phase: str) -> float:
    """Cumulative time at end of a phase."""
    cumulative = 0.0
    for p in PHASE_ORDER:
        cumulative += total_time * PHASE_FRACTIONS[p]
        if p == phase:
            return cumulative
    return cumulative


def _get_current_phase(cycle_clock: float, total_time: float) -> str:
    """Determine cell cycle phase from clock position."""
    cumulative = 0.0
    for p in PHASE_ORDER:
        cumulative += total_time * PHASE_FRACTIONS[p]
        if cycle_clock < cumulative:
            return p
    return "M"


def sample_environment(
    cell: Cell,
    env: EnvironmentFields,
    domain: Domain,
    hypoxia_threshold: float,
) -> None:
    """Sample local environment at cell position and update cell state."""
    cell.O2_local = domain.bilinear_interpolate(env.O2, cell.x, cell.y)
    cell.glucose_local = domain.bilinear_interpolate(env.glucose, cell.x, cell.y)
    cell.is_hypoxic = cell.O2_local < hypoxia_threshold


def update_effective_rates(cell: Cell, tp: CellTypeConfig) -> None:
    """Update ecDNA-modulated effective rates on cell."""
    cell.migration_rate = modulate_migration_speed(
        tp.migration_speed_um_hr, cell.ecDNA_count, tp.ecDNA_effect_on_migration
    )
    cell.VEGF_secretion = modulate_vegf_secretion(
        tp.VEGF_secretion_amol_hr, cell.ecDNA_count, tp.ecDNA_effect_on_VEGF
    )


def check_necrosis(cell: Cell, tp: CellTypeConfig) -> bool:
    """Check if cell should become necrotic due to sustained hypoxia.

    Returns True if cell transitions to necrotic.
    """
    if cell.cell_type == NECROTIC:
        return False
    if tp.O2_necrosis_threshold_mmHg <= 0:
        return False

    if cell.O2_local < tp.O2_necrosis_threshold_mmHg:
        cell.necrosis_timer_hr += 1.0
        if cell.necrosis_timer_hr >= tp.necrosis_delay_hr:
            return True
    else:
        cell.necrosis_timer_hr = 0.0
    return False


def apply_necrosis(cell: Cell, config: SimulationConfig) -> None:
    """Transition cell to necrotic state."""
    cell.cell_type = NECROTIC
    cell.cell_type_name = "Necrotic"
    cell.cell_cycle_phase = "G0"
    cell.is_quiescent = True
    necrotic_params = config.cell_types.get(NECROTIC)
    if necrotic_params:
        cell.shape_path = necrotic_params.shape_path
        cell.migration_rate = 0.0
        cell.VEGF_secretion = 0.0


def check_apoptosis(cell: Cell, tp: CellTypeConfig, rng: np.random.Generator) -> bool:
    """Stochastic apoptosis check. Returns True if cell dies."""
    if tp.apoptosis_rate_per_hr <= 0:
        return False
    rate = modulate_apoptosis_rate(
        tp.apoptosis_rate_per_hr, cell.ecDNA_count, tp.ecDNA_effect_on_survival
    )
    return rng.random() < rate


def check_immune_kill(
    tumor_cell: Cell,
    immune_cell: Cell,
    kill_rate: float,
    rng: np.random.Generator,
) -> bool:
    """Check if immune cell kills tumor cell this step."""
    p_kill = kill_rate * max(immune_cell.activation_level, 0.1)
    return rng.random() < p_kill


def can_proliferate(
    cell: Cell,
    tp: CellTypeConfig,
    population: CellPopulation,
) -> bool:
    """Check resource gate for proliferation eligibility.

    Implements homeostatic constraints for adult brain tissue:
    - Microglia: only proliferate when activated (near tumor/inflammation)
    - Astrocytes: only proliferate when reactive (near tumor/injury)
    - Endothelial: only proliferate when sprouting (VEGF-stimulated angiogenesis)
    - Other normal cells: standard resource-gated proliferation
    """
    if not tp.can_divide:
        return False
    if tp.max_generations >= 0 and cell.generation >= tp.max_generations:
        return False
    if cell.O2_local < tp.O2_prolif_threshold_mmHg:
        return False
    if cell.glucose_local < tp.glucose_prolif_threshold_mM:
        return False

    # Homeostatic constraints: normal brain cells are quiescent unless activated
    # Microglia require activation (chemokine/inflammation) to proliferate
    if cell.cell_type == MICROGLIA:
        if not cell.is_reactive and cell.activation_level < 0.3:
            return False

    # Astrocytes require reactive state (tumor proximity) to proliferate
    if cell.cell_type == ASTROCYTE:
        if not cell.is_reactive:
            return False

    # Endothelial cells require VEGF stimulus (is_reactive set by vegf_above_threshold rule)
    if cell.cell_type == ENDOTHELIAL:
        if not cell.is_reactive:
            return False

    return True


def advance_cell_cycle(
    cell: Cell,
    tp: CellTypeConfig,
    population: CellPopulation,
    config: SimulationConfig,
    rng: np.random.Generator,
    current_time_hr: float,
) -> Optional[LineageRecord]:
    """Advance cell cycle by 1 hour. Returns LineageRecord if division occurred."""
    if not tp.can_divide or cell.cell_type == NECROTIC:
        return None

    eligible = can_proliferate(cell, tp, population)

    # Cell in G0: check if conditions allow re-entry
    if cell.cell_cycle_phase == "G0":
        if not eligible:
            cell.is_quiescent = True
            return None
        # Enter G1
        cell.cell_cycle_phase = "G1"
        cell.cycle_clock_hr = 0.0
        base_time = tp.division_time_mean_hr
        effective_time = modulate_division_time(
            base_time, cell.ecDNA_count, tp.ecDNA_effect_on_division
        )
        cell.total_cycle_time_hr = max(
            1.0, rng.normal(effective_time, tp.division_time_std_hr)
        )
        cell.is_quiescent = False
        return None

    # Active cycling: advance clock
    old_phase = cell.cell_cycle_phase

    # If in S/G2/M and resources drop, arrest (clock pauses)
    if old_phase in ("S", "G2", "M") and not eligible:
        return None

    # If in G1 and resources drop, go to G0
    if old_phase == "G1" and not eligible:
        cell.cell_cycle_phase = "G0"
        cell.is_quiescent = True
        return None

    cell.cycle_clock_hr += 1.0
    new_phase = _get_current_phase(cell.cycle_clock_hr, cell.total_cycle_time_hr)

    # Check if M-phase is complete
    if cell.cycle_clock_hr >= cell.total_cycle_time_hr:
        return _execute_division(cell, tp, population, config, rng, current_time_hr)

    cell.cell_cycle_phase = new_phase
    return None


def _execute_division(
    parent: Cell,
    tp: CellTypeConfig,
    population: CellPopulation,
    config: SimulationConfig,
    rng: np.random.Generator,
    current_time_hr: float,
) -> Optional[LineageRecord]:
    """Execute cell division, create daughter, partition ecDNA."""
    # Find empty adjacent position — tumor cells can displace neighbors
    # Use shape-aware collision detection
    if not tp.contact_inhibited:
        result = population.find_division_position_with_displacement(
            parent.x, parent.y, parent.diameter, parent.cell_id, rng,
            shape_path=parent.shape_path, angle=parent.angle,
        )
    else:
        positions = population.get_adjacent_empty_positions(
            parent.x, parent.y, parent.diameter, rng,
            shape_path=parent.shape_path, angle=parent.angle,
        )
        result = positions[0] if positions else None

    if result is None:
        # Contact inhibition: enter G0
        parent.cell_cycle_phase = "G0"
        parent.is_quiescent = True
        return None

    daughter_x, daughter_y = result

    # ecDNA segregation
    parent_ecDNA_before = parent.ecDNA_count
    parent_new, daughter_ecDNA = segregate_ecdna(
        parent.ecDNA_count, tp.ecDNA_segregation_p, rng
    )
    parent.ecDNA_count = parent_new

    # Create daughter cell
    daughter_id = population.allocate_id()
    daughter = create_cell(
        cell_id=daughter_id,
        parent_id=parent.cell_id,
        cell_type=parent.cell_type,
        x=daughter_x,
        y=daughter_y,
        time_born_hr=current_time_hr,
        type_params=tp,
        rng=rng,
        ecDNA_count=daughter_ecDNA,
        ecDNA_cargo=parent.ecDNA_cargo,
        generation=parent.generation + 1,
        cell_cycle_phase="G1",
    )
    population.add_cell(daughter)

    # Reset parent
    parent.generation += 1
    parent.cell_cycle_phase = "G1"
    parent.cycle_clock_hr = 0.0
    base_time = tp.division_time_mean_hr
    parent.total_cycle_time_hr = max(
        1.0,
        rng.normal(
            modulate_division_time(base_time, parent.ecDNA_count, tp.ecDNA_effect_on_division),
            tp.division_time_std_hr,
        ),
    )

    # Update effective rates for both cells
    update_effective_rates(parent, tp)
    update_effective_rates(daughter, tp)

    return LineageRecord(
        time_hr=current_time_hr,
        parent_id=parent.cell_id,
        parent_ecDNA_before=parent_ecDNA_before,
        parent_ecDNA_after=parent_new,
        daughter_id=daughter_id,
        daughter_ecDNA=daughter_ecDNA,
        x=parent.x,
        y=parent.y,
        generation=daughter.generation,
    )


def compute_migration(
    cell: Cell,
    tp: CellTypeConfig,
    env: EnvironmentFields,
    domain: Domain,
    population: CellPopulation,
    rng: np.random.Generator,
) -> None:
    """Compute and apply migration for a single cell."""
    if cell.migration_rate <= 0 or cell.cell_type == NECROTIC:
        return

    # Pericytes: constrained to stay near endothelial cells
    if cell.cell_type == PERICYTE:
        neighbors = population.get_neighbors(cell.x, cell.y, 30)
        has_endo = any(n.cell_type == ENDOTHELIAL for n in neighbors)
        if not has_endo:
            return

    # Compute direction
    # 1. Persistence (directional memory)
    persistence = math.exp(-1.0 / max(tp.migration_persistence_hr, 0.01))
    vx = persistence * math.cos(cell.persist_direction)
    vy = persistence * math.sin(cell.persist_direction)

    # 2. Chemotaxis gradients
    if tp.chemotaxis_O2 != 0:
        grad = domain.compute_gradient(env.O2, cell.x, cell.y)
        vx += tp.chemotaxis_O2 * grad[0]
        vy += tp.chemotaxis_O2 * grad[1]

    if tp.chemotaxis_VEGF != 0:
        grad = domain.compute_gradient(env.VEGF, cell.x, cell.y)
        vx += tp.chemotaxis_VEGF * grad[0]
        vy += tp.chemotaxis_VEGF * grad[1]

    # 3. Haptotaxis (ECM gradient)
    if tp.haptotaxis_ECM != 0:
        grad = domain.compute_gradient(env.ECM_density, cell.x, cell.y)
        vx += tp.haptotaxis_ECM * grad[0]
        vy += tp.haptotaxis_ECM * grad[1]

    # 4. Random noise
    noise_angle = rng.uniform(0, 2 * math.pi)
    noise_strength = 0.3
    vx += noise_strength * math.cos(noise_angle)
    vy += noise_strength * math.sin(noise_angle)

    # Compute direction angle
    direction = math.atan2(vy, vx)

    # Displacement
    speed = cell.migration_rate
    dx = speed * math.cos(direction)
    dy = speed * math.sin(direction)

    new_x = int(round(cell.x + dx))
    new_y = int(round(cell.y + dy))

    # Clamp to domain
    new_x, new_y = domain.clamp_position(new_x, new_y)

    # Volume exclusion: check collision at target using shape-aware detection
    collision = population.has_neighbor_within(
        new_x, new_y,
        cell.diameter * 0.8,  # Search radius
        exclude_id=cell.cell_id,
        query_shape=cell.shape_path,
        query_diameter=cell.diameter,
        query_angle=direction,  # Use movement direction as new angle
    )
    if not collision:
        old_x, old_y = cell.x, cell.y
        cell.x = new_x
        cell.y = new_y
        cell.angle = direction  # Update cell orientation to match movement
        population.update_position(cell, old_x, old_y)

    # Update persistence direction
    cell.persist_direction = direction


def check_state_transitions(
    cell: Cell,
    tp: CellTypeConfig,
    env: EnvironmentFields,
    domain: Domain,
    population: CellPopulation,
    config: SimulationConfig,
    rng: np.random.Generator,
) -> None:
    """Check and apply state transition rules (Table 5.4)."""
    if cell.cell_type == NECROTIC:
        return

    for rule in tp.transition_rules:
        triggered = False

        if rule.condition == "adjacent_to_tumor":
            neighbors = population.get_neighbors(cell.x, cell.y, cell.diameter * 2)
            triggered = any(n.cell_type == TUMOR and n.cell_id != cell.cell_id for n in neighbors)

        elif rule.condition == "vegf_above_threshold":
            local_vegf = domain.bilinear_interpolate(env.VEGF, cell.x, cell.y)
            triggered = local_vegf > config.angiogenesis.angiogenesis_threshold_nM

        elif rule.condition == "chemokine_above_threshold":
            # Approximate chemokine as proportional to nearby tumor/necrotic cell density
            neighbors = population.get_neighbors(cell.x, cell.y, 100)
            tumor_nearby = sum(1 for n in neighbors if n.cell_type in (TUMOR, NECROTIC))
            triggered = tumor_nearby > 2

        elif rule.condition == "tgfb_above_threshold":
            # Approximate TGF-beta from tumor density
            neighbors = population.get_neighbors(cell.x, cell.y, 80)
            tumor_nearby = sum(1 for n in neighbors if n.cell_type == TUMOR)
            triggered = tumor_nearby > 5

        if triggered and rng.random() < rule.rate_per_hr:
            if rule.subtype_flag:
                cell.is_reactive = True
                cell.activation_level = min(cell.activation_level + 0.2, 1.0)
                # Microglia: increase migration when activated
                if cell.cell_type == MICROGLIA and rule.subtype_flag == "activated":
                    cell.migration_rate = tp.migration_speed_um_hr * 1.5


def clear_necrotic_cells(
    population: CellPopulation,
    config: SimulationConfig,
    rng: np.random.Generator,
) -> list[int]:
    """Stochastic lysis of necrotic cells. Returns list of removed cell IDs."""
    necrotic_tp = config.cell_types.get(NECROTIC)
    if necrotic_tp is None:
        return []

    lysis_rate = necrotic_tp.lysis_rate_per_hr
    to_remove: list[int] = []
    for cell in population.iter_by_type(NECROTIC):
        if rng.random() < lysis_rate:
            to_remove.append(cell.cell_id)

    for cid in to_remove:
        population.remove_cell(cid)
    return to_remove
