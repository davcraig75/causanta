# CAUSANTA

**Causal Analysis Using Somatic And Neighborhood Tissue Architecture**

CAUSANTA is a 2D tissue simulation engine that generates synthetic tissue sections with known ground-truth causal structure, enabling rigorous benchmarking of causal discovery algorithms. It treats heritable somatic variation (particularly extrachromosomal DNA) as Somatic Instrumental Variables (SIVs) within a Structural Causal Model, allowing identification of causal effects from otherwise observational spatial data.

## Features

- **9 cell types** modeled as individual agents: Neuron, Astrocyte, Oligodendrocyte, Microglia, Endothelial, Pericyte, Tumor (Glioma), Recruited Immune, and Necrotic
- **Reaction-diffusion environment** with O2, glucose, VEGF, and lactate fields solved via implicit LOD operator splitting (unconditionally stable)
- **ecDNA segregation** as a natural randomization mechanism: binomial partitioning during mitosis generates inter-cell heterogeneity that satisfies instrumental variable assumptions
- **Ground-truth causal graph** embedded in the simulation: ecDNA copy number causally affects proliferation rate, VEGF secretion, migration speed, and apoptosis resistance
- **Cell behaviors**: 5-phase cell cycle, chemotaxis/haptotaxis migration, necrosis, apoptosis, immune-mediated killing, state transitions (reactive astrocytes, activated microglia)
- **Angiogenesis**: VEGF-driven vessel sprouting with tip cell selection and vessel maturation
- **Vega.js visualization** with interactive time slider, environment heatmap overlay, and cell type filtering

## Quick Start

### Requirements

- Python >= 3.10
- NumPy >= 1.24

### Installation

```bash
pip install -e .
```

### Running a Simulation

```bash
# Default: 1mm x 1mm domain, 168 hours (1 week), single tumor seed
python -m causanta.simulation

# Custom duration
python -m causanta.simulation --hours 72

# Custom configuration
python -m causanta.simulation my_params.json

# Set RNG seed for reproducibility
python -m causanta.simulation --seed 123
```

### Output

Each run creates a timestamped directory under `output/`:

```
output/run_YYYYMMDD_HHMMSS/
    params.json                    # Frozen copy of parameters
    cells_t000000.tsv              # Cell state at each hour
    cells_t000001.tsv
    ...
    environment_t000000.tsv        # Environment grid at each hour
    environment_t000001.tsv
    ...
    lineage.tsv                    # Division events with ecDNA segregation
    summary.log                    # Per-step statistics
    visualization.vg.json          # Vega.js specification
    index.html                     # Self-contained HTML viewer
```

## Configuration

All parameters are specified in a single JSON file. See [`causanta/params/default.json`](causanta/params/default.json) for the complete default configuration.

Key configurable sections:

| Section | Parameters |
|---------|-----------|
| `domain` | `width_um`, `height_um`, `env_grid_um` (up to 12mm x 12mm) |
| `time` | `total_hours`, `dt_diffusion_hr`, `output_interval_hr`, `burnin_hours` |
| `environment` | Substrate diffusion coefficients, decay rates, boundary values |
| `cell_types` | Per-type geometry, proliferation, migration, metabolism, ecDNA effects |
| `tumor_seeds` | Initial tumor positions, ecDNA copy number, cargo genes |
| `angiogenesis` | VEGF threshold, sprouting rate, vessel maturation time |
| `immune_recruitment` | Chemokine threshold, recruitment rate, kill rate |

## Architecture

```
causanta/
    __init__.py
    config.py              # JSON parameter loading and validation
    domain.py              # Grid, coordinate system, bilinear interpolation
    cells.py               # Cell agent dataclass, spatial hash population manager
    ecdna.py               # ecDNA segregation and phenotype modulation
    environment.py         # Reaction-diffusion fields, vectorized Thomas algorithm
    behaviors.py           # Proliferation, migration, death, state transitions
    angiogenesis.py        # VEGF-driven vessel sprouting
    initialization.py      # Synthetic tissue setup (vasculature, cells, tumor)
    simulation.py          # Main loop orchestration and CLI entry point
    io.py                  # TSV/JSON file I/O
    visualization.py       # Vega.js spec and HTML viewer generation
    params/
        default.json       # Default parameter configuration
```

### Simulation Loop (per hour)

1. **Environment diffusion** -- Implicit LOD sub-stepping for O2, glucose, VEGF, lactate
2. **Cell behaviors** -- Death checks, state transitions, cell cycle advancement, migration (random-shuffled order)
3. **Angiogenesis** -- VEGF-driven vessel sprouting
4. **Immune recruitment** -- Spawn immune cells from vasculature near tumor signals
5. **Cleanup** -- Necrotic cell lysis, occupancy grid update
6. **Output** -- Write TSV snapshots at configured interval

### Key Design Decisions

- **Cell storage**: Individual `Cell` dataclass objects in a `dict[int, Cell]` with a spatial hash grid (20 um buckets) for O(1) neighbor queries
- **Diffusion solver**: Vectorized Thomas algorithm processes all grid rows/columns simultaneously via NumPy; implicit LOD scheme is unconditionally stable (no CFL constraint)
- **Source/sink accumulation**: Computed once per hour and reused across all diffusion sub-steps, since cells only act once per hour
- **Tumor displacement**: Non-contact-inhibited cells (tumor) push neighboring cells aside during division, matching invasive growth biology
- **Reproducibility**: Single seeded `np.random.Generator` passed through the entire simulation

## The Causal Inference Framework

CAUSANTA embeds a known Structural Causal Model:

```
ecDNA_count --> EGFR_expression --> proliferation_rate
     |                                      ^
     |                              O2_local (confounder)
     |                                      ^
     +-----> VEGF_secretion --> angiogenesis -+
```

The ground-truth causal effects are:

| Relationship | Equation |
|---|---|
| ecDNA -> division time | T_eff = T_base / (1 + alpha * log2(1 + ecDNA)) |
| ecDNA -> VEGF secretion | S_eff = S_base * (1 + beta * ecDNA) |
| ecDNA -> migration speed | v_eff = v_base * (1 + delta * ecDNA) |
| ecDNA -> apoptosis resistance | a_eff = a_base / (1 + gamma * ecDNA) |

ecDNA segregation during mitosis follows `Binomial(N, 0.5)`, providing the randomization that satisfies instrumental variable assumptions (relevance, independence, exclusion restriction).

## Performance

| Domain | Cells | Speed | 168h Runtime |
|---|---|---|---|
| 1mm x 1mm | ~5,000 | ~0.8 hr/s | ~4 min |

## License

MIT License. See [LICENSE](LICENSE) for details.
