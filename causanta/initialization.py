"""Tissue initialization for CAUSANTA.

Generates synthetic normal brain tissue with vascular network,
normal cell populations, and tumor seeding.
"""

from __future__ import annotations

import math

import numpy as np

from .cells import (
    ASTROCYTE,
    ENDOTHELIAL,
    MICROGLIA,
    NEURON,
    OLIGODENDROCYTE,
    PERICYTE,
    TUMOR,
    CellPopulation,
    create_cell,
)
from .config import SimulationConfig
from .domain import Domain
from .environment import DiffusionSolver, EnvironmentFields, accumulate_sources_and_sinks


def generate_vascular_network(
    domain: Domain,
    population: CellPopulation,
    env: EnvironmentFields,
    config: SimulationConfig,
    rng: np.random.Generator,
) -> None:
    """Generate a branching capillary network.

    Places endothelial cells along vessel paths and pericytes adjacent to them.
    Sets vascular_density field at vessel positions.
    """
    spacing = config.initialization.vascular_spacing_um
    endo_params = config.cell_types.get(ENDOTHELIAL)
    peri_params = config.cell_types.get(PERICYTE)
    if endo_params is None:
        return

    endo_diameter = endo_params.diameter_mean_um
    cell_spacing = int(endo_diameter)
    margin = 20  # um from domain edge

    # Place horizontal trunk vessels
    y_pos = int(spacing / 2)
    while y_pos < domain.height_um - margin:
        x = margin
        while x < domain.width_um - margin:
            cid = population.allocate_id()
            cell = create_cell(
                cell_id=cid,
                parent_id=-1,
                cell_type=ENDOTHELIAL,
                x=x,
                y=y_pos,
                time_born_hr=0.0,
                type_params=endo_params,
                rng=rng,
            )
            cell.cell_cycle_phase = "G0"
            cell.is_quiescent = True
            population.add_cell(cell)

            # Set vascular density
            gi, gj = domain.um_to_grid(x, y_pos)
            env.vascular_density[gj, gi] = max(env.vascular_density[gj, gi], 0.8)

            x += cell_spacing
        y_pos += int(spacing)

    # Place vertical connecting vessels
    x_pos = int(spacing / 2)
    while x_pos < domain.width_um - margin:
        y = margin
        while y < domain.height_um - margin:
            # Only place if not already occupied
            if not population.has_neighbor_within(x_pos, y, endo_diameter * 0.6):
                cid = population.allocate_id()
                cell = create_cell(
                    cell_id=cid,
                    parent_id=-1,
                    cell_type=ENDOTHELIAL,
                    x=x_pos,
                    y=y,
                    time_born_hr=0.0,
                    type_params=endo_params,
                    rng=rng,
                )
                cell.cell_cycle_phase = "G0"
                cell.is_quiescent = True
                population.add_cell(cell)

                gi, gj = domain.um_to_grid(x_pos, y)
                env.vascular_density[gj, gi] = max(env.vascular_density[gj, gi], 0.8)

            y += cell_spacing
        x_pos += int(spacing)

    # Add random branches
    n_branches = max(1, int(domain.width_um * domain.height_um / (spacing * spacing) * 2))
    for _ in range(n_branches):
        bx = rng.integers(margin, domain.width_um - margin)
        by = rng.integers(margin, domain.height_um - margin)
        angle = rng.uniform(0, 2 * math.pi)
        length = rng.integers(30, int(spacing * 0.6))

        for step in range(0, length, cell_spacing):
            px = int(bx + step * math.cos(angle))
            py = int(by + step * math.sin(angle))
            if not domain.is_in_bounds(px, py):
                break
            if population.has_neighbor_within(px, py, endo_diameter * 0.6):
                continue

            cid = population.allocate_id()
            cell = create_cell(
                cell_id=cid,
                parent_id=-1,
                cell_type=ENDOTHELIAL,
                x=px,
                y=py,
                time_born_hr=0.0,
                type_params=endo_params,
                rng=rng,
            )
            cell.cell_cycle_phase = "G0"
            cell.is_quiescent = True
            population.add_cell(cell)

            gi, gj = domain.um_to_grid(px, py)
            if 0 <= gi < domain.nx and 0 <= gj < domain.ny:
                env.vascular_density[gj, gi] = max(env.vascular_density[gj, gi], 0.6)

    # Place pericytes adjacent to endothelial cells
    if peri_params is not None:
        endo_cells = list(population.iter_by_type(ENDOTHELIAL))
        n_pericytes = max(1, len(endo_cells) // 4)  # ~25% coverage
        selected = rng.choice(len(endo_cells), size=min(n_pericytes, len(endo_cells)), replace=False)
        for idx in selected:
            ec = endo_cells[idx]
            # Place pericyte at small offset from endothelial cell
            offset = int(endo_diameter * 0.6)
            angle = rng.uniform(0, 2 * math.pi)
            px = int(ec.x + offset * math.cos(angle))
            py = int(ec.y + offset * math.sin(angle))
            px, py = domain.clamp_position(px, py)

            if not population.has_neighbor_within(px, py, peri_params.diameter_mean_um * 0.6):
                cid = population.allocate_id()
                pcell = create_cell(
                    cell_id=cid,
                    parent_id=-1,
                    cell_type=PERICYTE,
                    x=px,
                    y=py,
                    time_born_hr=0.0,
                    type_params=peri_params,
                    rng=rng,
                )
                pcell.cell_cycle_phase = "G0"
                pcell.is_quiescent = True
                population.add_cell(pcell)


def populate_normal_tissue(
    domain: Domain,
    population: CellPopulation,
    config: SimulationConfig,
    rng: np.random.Generator,
) -> None:
    """Place normal cell types at biologically realistic proportions.

    Uses rejection sampling to avoid overlap with existing cells.
    """
    total_area = domain.width_um * domain.height_um
    target_total = int(total_area * config.initialization.cell_density_per_um2)

    # Subtract already-placed cells (endothelial, pericytes)
    existing = population.count()
    remaining = max(0, target_total - existing)

    fractions = config.initialization.fractions
    type_map = {
        "neuron": NEURON,
        "astrocyte": ASTROCYTE,
        "oligodendrocyte": OLIGODENDROCYTE,
        "microglia": MICROGLIA,
    }

    # Normalize fractions for the remaining cell types
    frac_sum = sum(fractions.get(name, 0) for name in type_map)
    if frac_sum <= 0:
        return

    max_attempts = 50

    for name, type_id in type_map.items():
        frac = fractions.get(name, 0)
        if frac <= 0:
            continue

        tp = config.cell_types.get(type_id)
        if tp is None:
            continue

        n_cells = int(remaining * frac / frac_sum)
        placed = 0
        min_dist = tp.diameter_mean_um * 0.8

        for _ in range(n_cells):
            for attempt in range(max_attempts):
                x = rng.integers(10, domain.width_um - 10)
                y = rng.integers(10, domain.height_um - 10)
                if not population.has_neighbor_within(x, y, min_dist):
                    cid = population.allocate_id()
                    cell = create_cell(
                        cell_id=cid,
                        parent_id=-1,
                        cell_type=type_id,
                        x=int(x),
                        y=int(y),
                        time_born_hr=0.0,
                        type_params=tp,
                        rng=rng,
                    )
                    # Start normal cells in G0 (quiescent)
                    cell.cell_cycle_phase = "G0"
                    cell.is_quiescent = True
                    population.add_cell(cell)
                    placed += 1
                    break


def seed_tumor(
    population: CellPopulation,
    config: SimulationConfig,
    rng: np.random.Generator,
    current_time_hr: float,
) -> None:
    """Place initial tumor cells at configured positions with ecDNA."""
    tp = config.cell_types.get(TUMOR)
    if tp is None:
        return

    for seed in config.tumor_seeds:
        # Clear any existing cell at tumor position
        cid = population.allocate_id()
        cargo_str = ";".join(seed.ecDNA_cargo)
        cell = create_cell(
            cell_id=cid,
            parent_id=-1,
            cell_type=TUMOR,
            x=seed.x,
            y=seed.y,
            time_born_hr=current_time_hr,
            type_params=tp,
            rng=rng,
            ecDNA_count=seed.ecDNA_count,
            ecDNA_cargo=cargo_str,
            cell_cycle_phase=seed.cell_cycle_phase,
        )
        population.add_cell(cell)


def equilibrate_environment(
    env: EnvironmentFields,
    solver: DiffusionSolver,
    population: CellPopulation,
    domain: Domain,
    config: SimulationConfig,
) -> None:
    """Run diffusion for burn-in period to reach steady-state fields.

    Only diffusion and source/sink; no cell behaviors during burn-in.
    """
    burnin_hours = config.time.burnin_hours
    dt = config.time.dt_diffusion_hr
    n_substeps = max(1, int(1.0 / dt))

    print(f"  Equilibrating environment ({burnin_hours} hours)...")
    sources, sinks = accumulate_sources_and_sinks(population, env, domain, config)

    for hour in range(burnin_hours):
        solver.step(env, sources, sinks, dt, n_substeps)
        if (hour + 1) % 20 == 0:
            mean_o2 = float(env.O2.mean())
            print(f"    Burn-in hour {hour + 1}/{burnin_hours}: mean O2 = {mean_o2:.1f} mmHg")
