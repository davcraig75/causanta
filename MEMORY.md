# CAUSANTA Development Log

## Project Status

Simulation engine v0.1.0 fully implemented with 12 modules.

## Architecture Overview

- **Package**: `causanta/` (Python 3.10+, NumPy)
- **CLI**: `python -m causanta.simulation [config.json] [--hours N] [--seed N]`
- **Default config**: `causanta/params/default.json` (1mm domain, 168 hours, 9 cell types)
- **Output**: Timestamped directories under `output/` with TSV snapshots, lineage records, and Vega.js visualization

## Module Summary

| Module | Responsibility |
|---|---|
| `config.py` | Load and validate JSON parameters into frozen dataclass hierarchy |
| `domain.py` | 2D spatial domain, coordinate transforms, bilinear interpolation |
| `cells.py` | `Cell` dataclass (30+ fields), `CellPopulation` with spatial hash |
| `ecdna.py` | Binomial segregation, phenotype modulation functions |
| `environment.py` | Reaction-diffusion solver (vectorized Thomas algorithm, implicit LOD) |
| `behaviors.py` | Cell cycle, migration, death, state transitions |
| `angiogenesis.py` | VEGF-driven vessel sprouting |
| `initialization.py` | Vascular network, tissue population, tumor seeding, burn-in |
| `simulation.py` | Main loop (6 phases/hour), CLI entry point |
| `io.py` | TSV I/O for cells, environment, lineage, summary |
| `visualization.py` | Vega.js spec and HTML viewer generation |

## Key Design Decisions

1. **Cell storage**: `dict[int, Cell]` + spatial hash grid (20 um buckets) for O(1) neighbor queries
2. **Diffusion**: Vectorized Thomas algorithm (all rows/columns as NumPy batch), implicit LOD (unconditionally stable, no CFL constraint)
3. **Source/sink**: Accumulated once per hour, reused across all diffusion sub-steps
4. **Tumor division**: Non-contact-inhibited cells displace neighbors via `find_division_position_with_displacement()`
5. **Reproducibility**: Single `np.random.Generator` seeded from config, passed through entire simulation

## Conventions

- All rates are per hour (the simulation time unit)
- Environment grid indexing: `field[row, col]` = `field[gj, gi]`
- Cell coordinates: integer micrometers, 1 um precision
- Cell IDs: monotonically increasing integers, allocated via `CellPopulation.allocate_id()`

## Performance Baseline

- 1mm x 1mm domain, ~5,000 initial cells: ~0.8 simulated hours/sec
- 168-hour simulation completes in ~255 seconds
- Bottleneck: diffusion solver (100 sub-steps/hour x 4 substrates)
