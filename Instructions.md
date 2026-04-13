# CAUSANTA: Complete Rebuild Instructions

## For Reconstructing the Project from Scratch in Any LLM

**Version:** 0.2.0
**Date:** 2026-04-13

---

## 1. What Is CAUSANTA?

**CAUSANTA** (Causal Analysis Using Somatic And Neighborhood Tissue Architecture) is a Python simulation framework that generates synthetic tumor tissue with **known ground-truth causal structure**. Its purpose is benchmarking causal inference methods in spatial biology.

### The Core Scientific Insight

Extrachromosomal DNA (ecDNA) segregates randomly during mitosis via `Binomial(N, 0.5)`. This randomness satisfies the three instrumental variable (IV) assumptions:

1. **Relevance**: ecDNA copy number determines EGFR expression (gene dosage)
2. **Independence**: Segregation is random, uncorrelated with confounders like hypoxia
3. **Exclusion**: ecDNA affects phenotype only through gene expression

This makes ecDNA a **Somatic Instrumental Variable (SIV)** — a natural randomization mechanism embedded in the biology. The simulation generates data where you can validate whether causal inference methods (2SLS, IV regression, propensity matching, etc.) correctly recover the true causal effects.

### The Causal DAG (Ground Truth)

```
ecDNA_count (INSTRUMENT, random segregation)
    |
    v
EGFR_expression (EXPOSURE, gene dosage)
    |
    |---> Proliferation rate    (alpha = 0.3)
    |---> Migration speed       (delta = 0.05)
    |---> VEGF secretion        (beta = 0.1)
    |---> Apoptosis resistance  (gamma = 0.2)
              ^
              |
        O2_local (CONFOUNDER, affects outcomes but NOT ecDNA)
```

The analyst's challenge: recover alpha, beta, delta, gamma from observational data, despite O2 confounding.

---

## 2. Technology Stack

```
Python >= 3.10
numpy >= 1.24          # Core numerics
scipy >= 1.10          # Statistics, optimization
matplotlib >= 3.7      # Visualization
numba >= 0.58          # Optional JIT for diffusion solver (~2x speedup)
```

**Build system:** setuptools >= 68.0 + wheel

**CLI entry points:**
- `causanta-simulate` -> `causanta.simulate.core:main`
- `causanta-analyze` -> `causanta.analyze.cli:main`

---

## 3. Repository Structure

```
causanta/                          # Root (also the git repo name)
├── pyproject.toml                 # Package metadata, deps, entry points
├── README.md                      # Project overview with causal DAG
├── LICENSE                        # MIT
├── .gitignore
│
├── causanta/                      # Python package (flat layout, standard)
│   ├── __init__.py                # Exports: simulate, analyze, graph submodules
│   │
│   ├── simulate/                  # Tumor growth simulation engine
│   │   ├── __init__.py
│   │   ├── core.py                # Main loop, CLI
│   │   ├── config.py              # Frozen dataclass hierarchy from JSON
│   │   ├── domain.py              # 2D spatial domain, interpolation
│   │   ├── cells.py               # Cell dataclass, CellPopulation, spatial hash
│   │   ├── ecdna.py               # Segregation, EGFR expression, phenotype modulation
│   │   ├── environment.py         # Reaction-diffusion PDE solver
│   │   ├── behaviors.py           # Division, migration, death, immune killing
│   │   ├── angiogenesis.py        # VEGF-driven vessel sprouting
│   │   ├── initialization.py      # Vascular network, tissue, tumor seeding
│   │   ├── io.py                  # TSV read/write for cells, environment, lineage
│   │   ├── visualization.py       # Vega-Lite spec generation
│   │   ├── reporting.py           # HTML simulation report
│   │   ├── viewer.py              # Interactive SVG cell viewer
│   │   ├── movie.py               # Animated HTML tumor growth
│   │   ├── shapes.py              # SVG path parsing, collision detection
│   │   ├── cleanup.py             # Reduce intermediate TSV files
│   │   ├── benchmark.py           # Performance testing
│   │   ├── _numba_kernels.py      # JIT-compiled Thomas algorithm
│   │   └── params/
│   │       └── default.json       # All simulation parameters
│   │
│   ├── analyze/                   # Causal inference analysis tools
│   │   ├── __init__.py
│   │   ├── cli.py                 # CLI: causanta-analyze
│   │   ├── loader.py              # Load TSV data -> SimulationData
│   │   ├── effects.py             # Estimate alpha, beta, delta, gamma from data
│   │   ├── iv.py                  # Two-Stage Least Squares (2SLS)
│   │   ├── regression.py          # OLS, polynomial, doubly robust
│   │   ├── discovery.py           # PC algorithm, GES algorithm
│   │   ├── matching.py            # Propensity score matching, IPW
│   │   ├── sem.py                 # Structural equation modeling
│   │   ├── sensitivity.py         # Rosenbaum bounds, E-values, OVB
│   │   └── report.py              # HTML analysis report generator
│   │
│   └── graph/                     # Causal DAG representation
│       ├── __init__.py
│       └── causal_dag.py          # CausalDAG, nodes, edges, Mermaid/Graphviz
│
├── docs/                          # Documentation (Markdown + rendered HTML)
│   ├── tutorial.md / .html        # Causal inference tutorial
│   ├── analysis_guide.md / .html  # Statistical methods guide
│   ├── immune_system.md / .html   # Immune system documentation
│   └── assets/
│       └── style.css              # Unified documentation CSS
│
├── examples/
│   ├── README.md
│   └── analyze_simulation.py      # End-to-end analysis example
│
├── scripts/
│   └── render_docs.py             # Markdown -> HTML renderer
│
├── old/                           # Archived stale files
│
└── output/                        # Simulation output (gitignored)
    └── run_YYYYMMDD_HHMMSS/
        ├── data/                  # TSV snapshots, lineage, params
        ├── figures/               # Vega-Lite specs
        ├── reports/               # HTML reports
        └── animations/            # Tumor growth movie
```

---

## 4. Simulation Engine — Module-by-Module

### 4.1 Configuration System (`config.py`)

All parameters live in a JSON file. The config module loads them into a hierarchy of **frozen dataclasses**.

**Hierarchy:**

```
SimulationConfig
├── rng_seed: int = 42
├── domain: DomainConfig
│   ├── width_um: int = 1000
│   ├── height_um: int = 1000
│   └── env_grid_um: int = 10
├── time: TimeConfig
│   ├── total_hours: int = 168
│   ├── dt_diffusion_hr: float = 0.01
│   ├── output_interval_hr: int = 1
│   └── burnin_hours: int = 20
├── environment: EnvironmentConfig
│   ├── O2_blood_mmHg: float = 60.0
│   ├── q_O2_transfer_per_hr: float = 5.0
│   ├── q_glucose_transfer_per_hr: float = 10.0
│   ├── hypoxia_threshold_mmHg: float = 10.0
│   └── substrates: dict[str, SubstrateConfig]
│       ├── O2:      D=6e6, decay=0,    boundary=38, initial=38, clamp=[0,60]
│       ├── glucose:  D=2.4e6, decay=0,  boundary=5,  initial=5,  clamp=[0,5.5]
│       ├── VEGF:    D=36000, decay=0.6, boundary=0,  initial=0,  clamp=[0,100]
│       └── lactate: D=9e5,  decay=0.06, boundary=1,  initial=1,  clamp=[0,40]
├── angiogenesis: AngiogenesisConfig
│   ├── angiogenesis_threshold_nM: float = 5.0
│   ├── vessel_maturation_hr: float = 48.0
│   ├── sprouting_rate_per_hr: float = 0.01
│   └── max_sprout_length: int = 200
├── immune_recruitment: ImmuneRecruitmentConfig
│   ├── chemokine_threshold: 2.0
│   ├── recruitment_rate_per_hr: 0.003
│   ├── kill_radius_um: 15.0
│   ├── synapse_formation_time_hr: 2.5
│   ├── kill_probability_per_synapse: 0.25
│   ├── min_contact_for_kill_hr: 1.5
│   ├── max_kills_before_exhaustion: 3
│   ├── exhaustion_per_kill: 0.35
│   ├── exhaustion_recovery_rate_per_hr: 0.005
│   ├── exhausted_kill_penalty: 0.8
│   ├── activation_radius_um: 30.0
│   ├── activation_rate_per_hr: 0.1
│   ├── deactivation_rate_per_hr: 0.15
│   └── min_activation_for_kill: 0.5
├── initialization: InitConfig
│   ├── cell_density_per_um2: 0.005
│   ├── fractions: {neuron:0.45, astrocyte:0.25, oligo:0.12, microglia:0.07, endo:0.04, peri:0.02}
│   └── vascular_spacing_um: 150
├── tumor_seeds: list[TumorSeed]
│   └── [{x:525, y:525, ecDNA_count:20, ecDNA_cargo:["EGFR"], phase:"G1"}]
└── cell_types: dict[int, CellTypeConfig]   # 9 types (IDs 0-8)
```

**Cell Type IDs:**
```
0=Neuron, 1=Astrocyte, 2=Oligodendrocyte, 3=Microglia,
4=Endothelial, 5=Pericyte, 6=Tumor, 7=RecruitedImmune, 8=Necrotic
```

**Each CellTypeConfig contains:**
- Morphology: `diameter_mean_um`, `diameter_std_um`, `shape_path` (SVG), `color`
- Division: `can_divide`, `division_time_mean_hr`, `division_time_std_hr`, `max_generations`
- Thresholds: `O2_prolif_threshold_mmHg`, `glucose_prolif_threshold_mM`
- Migration: `migration_speed_um_hr`, `migration_persistence_hr`, `chemotaxis_O2/VEGF/chemokine`, `haptotaxis_ECM`
- Death: `O2_necrosis_threshold_mmHg`, `necrosis_delay_hr`, `apoptosis_rate_per_hr`, `lysis_rate_per_hr`
- Metabolism: `O2_consumption_amol_hr`, `glucose_consumption_amol_hr`, `VEGF_secretion_amol_hr`, `lactate_production_amol_hr`
- ecDNA: `ecDNA_init_count`, `ecDNA_cargo`, `ecDNA_segregation_p=0.5`, `ecDNA_effect_on_division/VEGF/migration/survival`
- Transitions: `transition_rules` (condition -> target_type at rate)
- `contact_inhibited: bool` (normal cells yes, tumor no)

**Key Tumor Cell Parameters (type 6):**
```json
{
    "division_time_mean_hr": 36.0,
    "migration_speed_um_hr": 10.0,
    "O2_consumption_amol_hr": 72000.0,
    "VEGF_secretion_amol_hr": 600.0,
    "ecDNA_init_count": 20,
    "ecDNA_effect_on_division": 0.3,
    "ecDNA_effect_on_VEGF": 0.1,
    "ecDNA_effect_on_migration": 0.05,
    "ecDNA_effect_on_survival": 0.2
}
```

### 4.2 Spatial Domain (`domain.py`)

A 2D rectangular domain with grid-based environment fields.

**Properties:**
- `width_um`, `height_um`: Physical size in micrometers
- `env_grid_um`: Grid cell spacing (e.g., 10 um)
- `nx = width_um // env_grid_um`, `ny = height_um // env_grid_um`

**Coordinate transforms:**
- `um_to_grid(x_um, y_um) -> (gi, gj)`: Clamped integer indices
- `grid_to_um(gi, gj) -> (x_um, y_um)`: Center of grid cell

**Bilinear interpolation** for smooth field sampling at arbitrary cell positions:
```
gx = x_um / grid_um - 0.5
gy = y_um / grid_um - 0.5
value = v00*(1-fx)*(1-fy) + v10*fx*(1-fy) + v01*(1-fx)*fy + v11*fx*fy
```

**Central-difference gradients** for chemotaxis:
```
dfdx = (field[gj, gi+1] - field[gj, gi-1]) / (2 * grid_um)
dfdy = (field[gj+1, gi] - field[gj-1, gi]) / (2 * grid_um)
```

### 4.3 Cell Agents (`cells.py`)

**Cell dataclass** (~30 fields):

| Category | Fields |
|----------|--------|
| Spatial | `x, y: int` (um), `angle: float` (radians) |
| Morphology | `shape_path: str` (SVG), `size_scale: float`, `diameter: float` |
| Cell cycle | `cell_cycle_phase: str` {G0,G1,S,G2,M}, `cycle_clock_hr: float`, `total_cycle_time_hr: float` |
| Genomics | `ecDNA_count: int`, `ecDNA_cargo: str`, `egfr_expression: float` |
| Lineage | `generation: int`, `time_born_hr: float`, `parent_id: int` |
| Environment | `O2_local: float`, `glucose_local: float`, `is_hypoxic: bool` |
| Phenotype | `migration_rate: float`, `VEGF_secretion: float` |
| Immune | `is_reactive: bool`, `activation_level: float`, `contact_target_id: int`, `contact_duration_hr: float`, `kills_performed: int`, `exhaustion_level: float` |
| State | `is_alive: bool`, `necrosis_timer_hr: float`, `persist_direction: float` |

**CellPopulation** stores cells in `dict[int, Cell]` with a spatial hash grid:
- Bucket size: `SPATIAL_HASH_RESOLUTION = 20 um`
- `get_neighbors(x, y, radius)` — O(1) average via hash lookup
- `has_neighbor_within(...)` — Shape-aware collision detection using SVG bounds
- `find_division_position_with_displacement(...)` — PhysiCell-style "budging": tries 8 adjacent positions, then cascades displacement up to 10 cells deep
- `shuffled_living(rng)` — Random iteration order each step (fairness)
- Cell IDs: monotonically increasing integers via `allocate_id()`

### 4.4 ecDNA and Gene Expression (`ecdna.py`)

This is the **causal engine** — where ground truth effects are implemented.

**Constants:**
```python
MAX_ECDNA_COPIES = 100
EGFR_BASE_EXPRESSION = 1.0
EGFR_PER_ECDNA_COPY = 0.5
EGFR_NOISE_CV = 0.1
```

**EGFR Expression (gene dosage):**
```
EGFR = (1.0 + 0.5 * min(ecDNA, 100)) * LogNormal(0, CV=0.1)
```

**ecDNA Segregation at Division:**
```python
def replicate_and_segregate_ecdna(parent_count, p=0.5, rng, replication_fidelity=0.95):
    # S-phase: each copy replicates with 95% fidelity
    replicated = parent_count
    for i in range(parent_count):
        if rng.random() < 0.95:
            replicated += 1
    replicated = min(replicated, MAX_ECDNA_COPIES)

    # M-phase: random binomial segregation
    daughter = rng.binomial(replicated, p)
    parent_new = replicated - daughter
    return parent_new, daughter
```

**The Four Causal Effect Functions:**

| Effect | Formula | Parameter | Default |
|--------|---------|-----------|---------|
| Division time | `T_eff = T_base / (1 + alpha * log2(1 + EGFR))` | alpha | 0.3 |
| VEGF secretion | `S_eff = S_base * (1 + beta * sqrt(EGFR))` | beta | 0.1 |
| Migration speed | `v_eff = v_base * (1 + delta * EGFR) * hypoxia_mult` | delta | 0.05 |
| Apoptosis rate | `a_eff = a_base / (1 + gamma * log2(1 + EGFR))` | gamma | 0.2 |

Where `hypoxia_mult = 2.0` if hypoxic, else `1.0` (the Go-or-Grow switch).

**Numerical examples (EGFR expression by ecDNA count):**
- ecDNA=0: EGFR ~ 1.0
- ecDNA=10: EGFR ~ 6.0 (6x amplification)
- ecDNA=20: EGFR ~ 11.0 (11x)
- ecDNA=50: EGFR ~ 26.0 (26x)

### 4.5 Environment / Diffusion Solver (`environment.py`)

**EnvironmentFields** holds 2D arrays (ny x nx) for each substrate:
- `O2`, `glucose`, `VEGF`, `lactate`
- Plus: `ECM_density`, `vascular_density`, `pH`
- Occupancy tracking: `is_occupied`, `occupant_cell_id`, `occupant_type`

**Reaction-Diffusion PDE:**
```
du/dt = D * laplacian(u) + S(x,t) - lambda * u
```

**Numerical Method: Implicit LOD (Locally One-Dimensional)**

Each diffusion sub-step (dt = 0.01 hr, ~100 sub-steps per hour per substrate):

1. **X-sweep** (implicit tridiagonal):
   ```
   r = D * dt / (2 * dx^2)
   -r * u_{i-1} + (1+2r) * u_i - r * u_{i+1} = RHS
   ```
2. **Y-sweep** (same, transposed)
3. **Reaction step** (explicit): `field += dt * (source - sink - decay * field)`
4. **Clamp** to `[clamp_min, clamp_max]`

**Thomas algorithm** solves each tridiagonal system. Vectorized to process all rows/columns as a NumPy batch. Optional Numba JIT for ~2x speedup.

**Unconditionally stable** — no CFL constraint on dt, allowing large diffusion coefficients.

**Source/Sink Accumulation:**
- **Vascular O2 supply:** `O2_source = q_O2 * vasc_density * (O2_blood - local_O2)`
- **Cell O2 consumption:** `O2_sink[gj,gi] += cell_type.O2_consumption / voxel_area`
- **VEGF secretion (hypoxia-gated):**
  - Hypoxic cells: full VEGF rate (HIF-1alpha active)
  - Normoxic cells: 20% basal rate
- Sources/sinks accumulated once per hour, reused across all sub-steps

### 4.6 Cell Behaviors (`behaviors.py`)

**Cell Cycle Phases and Fractions:**
```
G1: 40%, S: 30%, G2: 15%, M: 15% of total cycle time
```

**Proliferation:**
1. `can_proliferate()` checks: O2 > threshold, glucose > threshold, not at max generation
2. G0 cells try to enter G1 each step; sample `total_cycle_time = N(modulated_mean, std)`
3. Clock advances by 1.0 hr each step; phase determined by fraction of total time
4. Arrest: if resources drop during S/G2/M, clock pauses; if during G1, return to G0
5. Division at clock >= total_time: calls `_execute_division()`

**Division Execution:**
- Non-contact-inhibited (tumor): `find_division_position_with_displacement()` — pushes neighbors aside
- Contact-inhibited (normal): `get_adjacent_empty_positions()` — only empty spots
- ecDNA segregation: `replicate_and_segregate_ecdna(parent_count, p, rng)`
- Creates daughter cell at found position; parent resets to G1
- Records `LineageRecord` with before/after ecDNA

**Migration:**
```python
# Persistence (biased random walk)
vx += exp(-1/persistence_hr) * cos(persist_direction)

# Chemotaxis (gradient-following)
vx += chemotaxis_O2 * grad_O2_x
vx += chemotaxis_VEGF * grad_VEGF_x

# Random noise
noise_angle = rng.uniform(0, 2*pi)
vx += 0.3 * cos(noise_angle)

# Move
new_x = x + migration_rate * cos(direction)  # Rounded to integer um

# Volume exclusion: reject move if collision detected
```

**Death:**
- **Necrosis:** If `O2_local < necrosis_threshold` for `necrosis_delay_hr` hours -> transition to NECROTIC
- **Apoptosis:** Stochastic with rate `modulate_apoptosis_rate(base, EGFR, gamma)`
- **Necrotic lysis:** Stochastic removal at `lysis_rate_per_hr`

**Immune Killing:**
```python
# Contact tracking
if same_target: contact_duration += dt
else: reset contact, start new

# Kill probability
if contact_duration < min_contact_for_kill: p_kill = 0
else:
    synapse_progress = contact_duration / synapse_formation_time
    base_kill = kill_probability_per_synapse * synapse_progress
    exhaustion_penalty = exhaustion_level * exhausted_kill_penalty
    p_kill = base_kill * (1 - exhaustion_penalty) * activation_level * dt
```

**Immune Activation:**
- Tumor within `activation_radius_um`: `activation += activation_rate * dt`, set `is_reactive = True`
- No tumor nearby: `activation -= deactivation_rate * dt`
- Recovery when not in contact: `exhaustion -= recovery_rate * dt`

**State Transitions:**
- Condition-based rules (e.g., "adjacent_to_tumor" -> astrocyte becomes reactive)
- Applied stochastically at configured `rate_per_hr`

### 4.7 Initialization (`initialization.py`)

**Execution order:**
1. **Vascular network**: Horizontal + vertical lines at `vascular_spacing_um` intervals; random branches; pericyte coverage (~25%); sets `vascular_density` on environment grid
2. **Normal tissue**: Fill domain to `cell_density_per_um2` with type fractions; rejection sampling with collision detection; all start in G0
3. **Environment burn-in**: Run diffusion for `burnin_hours` (e.g., 20h) WITHOUT cell behaviors to reach O2/glucose steady state
4. **Tumor seeding**: Place tumor cells per `tumor_seeds` config

### 4.8 Main Simulation Loop (`core.py`)

**Six phases per hour:**
1. **Environment diffusion** — solve reaction-diffusion PDEs
2. **Cell behaviors** — iterate all cells in random order: sample environment, update rates, check death, advance cycle, compute migration, check transitions
3. **Angiogenesis** — VEGF-driven vessel sprouting
4. **Immune recruitment** — spawn new immune cells at vasculature if chemokine > threshold
5. **Cleanup** — lyse necrotic cells, update occupancy
6. **Output** — write TSV snapshots at `output_interval_hr`

**CLI:**
```bash
python -m causanta.simulate.core [config.json] [--hours N] [--seed N]
```

### 4.9 Angiogenesis (`angiogenesis.py`)

- Endothelial cells with local `VEGF > threshold` can sprout
- Sprouting probability: `sprouting_rate_per_hr` per eligible cell
- New endothelial cell placed up the VEGF gradient
- Sprout tracks age; `vascular_density` increases as `min(age/maturation_time, 1.0) * 0.8`
- Mature when `age >= vessel_maturation_hr` (48h)

### 4.10 I/O (`io.py`)

**Output directory structure:**
```
output/run_YYYYMMDD_HHMMSS/
├── data/
│   ├── cells_t000000.tsv      # Cell states per timestep
│   ├── environment_t000000.tsv # Grid fields per timestep
│   ├── lineage.tsv            # All division events (appended)
│   ├── summary.log            # Per-step statistics
│   └── params.json            # Copy of parameters used
├── figures/
│   └── simulation.vl.json     # Vega-Lite spec
├── reports/
│   ├── report.html            # Simulation report
│   └── viewer.html            # Interactive cell viewer
└── animations/
    └── tumor_growth.html      # Animated tumor movie
```

**TSV Columns:**

Cells: `cell_id, parent_id, cell_type, cell_type_name, x, y, angle, shape_path, size_scale, cell_cycle_phase, ecDNA_count, ecDNA_cargo, egfr_expression, O2_local, glucose_local, is_hypoxic, is_quiescent, generation, time_born_hr, migration_rate, VEGF_secretion, activation_level, exhaustion_level, kills_performed`

Environment: `x_grid, y_grid, x_um, y_um, O2, glucose, VEGF, lactate, pH, ECM_density, vascular_density, is_occupied, occupant_cell_id, occupant_type`

Lineage: `time_hr, parent_id, parent_ecDNA_before, parent_ecDNA_after, daughter_id, daughter_ecDNA, x, y, generation`

### 4.11 Movie Renderer (`movie.py`)

- Reads all cell TSV files for a run
- `render_mode='tumor_focus'`: Only renders tumor (type 6), immune (7), necrotic (8) cells + 10% random sample of normal tissue
- Frame subsampling: `max_frames` parameter limits total frames
- Generates self-contained HTML with SVG frames and JavaScript playback controls (play/pause, step, speed slider, frame scrubber)

### 4.12 Cleanup Utility (`cleanup.py`)

```python
cleanup_simulation_outputs(output_dir, keep_fraction=0.01)
```
- Keeps first timestep, last timestep, and `keep_fraction` uniformly sampled from the rest
- Reports files deleted and disk space freed
- Dry-run mode available

```python
organize_output_directory(source_dir)
```
- Moves files into `data/`, `figures/`, `reports/`, `animations/` subdirectories

### 4.13 Shape Collision Detection (`shapes.py`)

- Each cell type has an SVG path defining its shape
- `parse_svg_path_bounds(path)` -> axis-aligned bounding box (cached with `@lru_cache`)
- Collision: rotated AABB intersection test
- Default shapes: neurons have branched paths, astrocytes are star-like, tumor cells are irregular spheres, immune cells are circles, etc.

---

## 5. Analysis Module — Module-by-Module

### 5.1 Data Loading (`loader.py`)

`SimulationData` dataclass holds all loaded data:
```python
@dataclass
class SimulationData:
    output_dir: Path
    config: dict
    cells: dict[int, list[dict]]       # timestep -> cell records
    environment: dict[int, list[dict]]  # timestep -> grid records
    lineage: list[dict]                 # all division events
    summary: list[dict]                 # per-step stats
    timesteps: list[int]
```

Key properties: `final_cells`, `tumor_cells` (type 6 at final step), `n_divisions`.

By default loads only first + last timesteps (memory efficiency).

### 5.2 Ground Truth Effect Estimation (`effects.py`)

Estimates the four causal parameters from simulation output:

- **alpha (division):** From lineage data, uses `T_eff = T_base / (1 + alpha * log2(1 + ecDNA))` relationship
- **beta (VEGF):** Linear regression of VEGF_secretion on ecDNA_count for tumor cells; `beta = slope / intercept`
- **delta (migration):** Linear regression of migration_rate on ecDNA_count for tumor cells
- **gamma (survival):** Enrichment analysis comparing mean ecDNA of surviving cells vs all cells born
- **Segregation analysis:** Tests if daughter fractions have mean 0.5 (one-sample t-test)

### 5.3 Instrumental Variable Regression (`iv.py`)

**Two-Stage Least Squares (2SLS):**
```
Stage 1: D = gamma_0 + gamma_1*Z + eta           (instrument -> treatment)
Stage 2: Y = beta_0  + beta_1*D_hat + epsilon     (predicted treatment -> outcome)
```

- Reports first-stage F-statistic (weak instrument test: F > 10 required)
- SE correction using original (not predicted) residuals
- Flags `weak_instrument` if F < 10
- Computes OLS estimate for comparison (bias ratio = OLS/IV)
- Default: instrument=ecDNA, outcomes=[VEGF_secretion, migration_rate], covariates=[O2_local, glucose_local]

**Segregation IV:** Uses daughter ecDNA (random) as instrument for parent ecDNA -> generation (division rate proxy)

### 5.4 OLS Regression (`regression.py`)

**Linear regression with adjustment:**
```
Y = b0 + b1*Treatment + sum(bj*Covariate_j) + epsilon
```
- Computes coefficients via `(X'X)^{-1} X'y`
- Standard errors, t-stats, p-values, R-squared, adjusted R-squared
- Default covariates: O2_local, glucose_local

**Also implements:**
- Polynomial regression (non-linear treatment effects)
- Interaction regression (effect modification)
- Doubly robust estimation (outcome model + propensity weighting + bootstrap SE, n_boot=200)

### 5.5 Causal Discovery (`discovery.py`)

**PC Algorithm (constraint-based):**
1. Start with complete graph
2. For increasing conditioning set sizes (0 to `max_cond_size`):
   - Test conditional independence via partial correlation + Fisher's z-transform
   - Remove edges where independence holds at significance `alpha`
   - Record separating sets
3. Orient v-structures: if i--k--j, i and j not adjacent, k not in sep(i,j) -> i->k<-j

**GES Algorithm (score-based):**
1. Forward: greedily add edges maximizing BIC
2. Backward: greedily remove edges maximizing BIC
3. BIC = sum_j [LL_j - penalty * k * log(n) / 2]

**Comparison to ground truth:** Precision, recall, F1, structural Hamming distance (SHD)

### 5.6 Propensity Score Matching (`matching.py`)

- Binary treatment: ecDNA > threshold (default 10)
- Propensity scores via logistic regression (gradient descent, no external deps)
- Nearest-neighbor matching with caliper (default 0.1)
- Reports: ATT, ATE, standardized mean differences before/after matching
- Also: Inverse Propensity Weighting (IPW) with bootstrap SE

### 5.7 Structural Equation Modeling (`sem.py`)

**Path analysis:**
```
Default model:
  VEGF = f(ecDNA, O2)
  migration = f(ecDNA, VEGF, O2)
```
- Decomposes total effect into direct + indirect (via mediator)
- Sobel test for indirect effect significance: `SE_indirect = sqrt(a^2*SE_b^2 + b^2*SE_a^2)`
- Reports path coefficients, z-stats, p-values

### 5.8 Sensitivity Analysis (`sensitivity.py`)

**Rosenbaum bounds:** How much hidden bias (odds ratio Gamma) would be needed to explain away the treatment effect? Reports breakdown point (first Gamma where p > 0.05).

**E-value:** Minimum strength of association an unmeasured confounder needs with both treatment and outcome to explain the observed effect. `E = RR + sqrt(RR*(RR-1))`

**Omitted variable bias** (Cinelli & Hazlett 2020): What fraction of residual variance would an omitted confounder need to explain?

**Placebo tests:** Regress treatment on outcomes it shouldn't affect (e.g., cell position x, y). Significant effects suggest confounding.

### 5.9 Causal DAG (`graph/causal_dag.py`)

**Data structures:**
- `CausalNode`: name, description, NodeType (INSTRUMENT/EXPOSURE/MEDIATOR/OUTCOME/CONFOUNDER)
- `CausalEdge`: source, target, EffectType (LINEAR/LOG/SQRT/MULTIPLICATIVE/INVERSE/THRESHOLD), parameter symbol (Greek letter), equation template
- `CausalDAG`: collection of nodes and edges with graph operations

**`build_causanta_dag(params)`** constructs the full ground-truth DAG with 9 nodes and 12 edges:

| Edge | Type | Parameter | Symbol |
|------|------|-----------|--------|
| ecDNA -> EGFR_expression | LINEAR | EGFR_per_ecDNA_copy | kappa |
| EGFR -> Proliferation | LOG | ecDNA_effect_on_division | alpha |
| EGFR -> Invasion | MULTIPLICATIVE | ecDNA_effect_on_migration | delta |
| EGFR -> VEGF | SQRT | ecDNA_effect_on_VEGF | beta |
| EGFR -> Survival | LOG | ecDNA_effect_on_survival | gamma |
| Proliferation -> O2_local | LINEAR | O2_consumption | omega |
| O2_local -> Hypoxia | THRESHOLD | hypoxia_threshold | tau |
| Hypoxia -> Invasion | MULTIPLICATIVE | hypoxia_invasion_boost | psi |
| Hypoxia -> Proliferation | INVERSE | hypoxia_prolif_suppression | eta |
| Hypoxia -> VEGF | MULTIPLICATIVE | hypoxia_VEGF_multiplier | phi |
| VEGF -> Vasculature | THRESHOLD | angiogenesis_threshold | theta |
| Vasculature -> O2_local | LINEAR | q_O2_transfer_per_hr | q |

**Output formats:** Mermaid diagram, Graphviz DOT, adjacency matrix

**`compare_dags()`** compares discovered graph against ground truth -> precision, recall, F1, SHD

### 5.10 Report Generation (`report.py`)

`generate_analysis_report()` runs all analysis methods and generates a comprehensive HTML report:
- Dark-themed UI with collapsible sections
- Interactive Mermaid DAG editor with presets (simple, full, default)
- Effect comparison charts (matplotlib -> base64 PNG)
- Tables for each analysis method's results
- Summary statistics

CLI: `causanta-analyze output_dir [--methods iv,discovery,...] [--format html]`

---

## 6. Key Design Decisions

1. **Cell storage**: `dict[int, Cell]` + spatial hash grid (20 um buckets) for O(1) neighbor queries. Not an array — cells are created/destroyed frequently.

2. **Diffusion**: Implicit LOD with vectorized Thomas algorithm. Unconditionally stable, no CFL constraint, allowing large diffusion coefficients (O2: D = 6e6 um^2/hr).

3. **Source/sink accumulation**: Computed once per hour, reused across ~100 diffusion sub-steps. Major performance optimization.

4. **Tumor division**: Non-contact-inhibited. Uses cascade displacement ("budging") up to 10 cells deep, mimicking PhysiCell's approach.

5. **Immune system**: GBM-realistic immunosuppressive defaults. Kill probability 25% (not 70%), exhaustion after ~3 kills, high activation threshold. Tumors should "win."

6. **Reproducibility**: Single `np.random.Generator` seeded from config, passed through entire simulation. Same seed -> same output.

7. **All rates are per hour** — the simulation's fundamental time unit.

8. **Grid indexing convention**: `field[row, col]` = `field[gj, gi]` (row-major, y-first).

9. **Cell coordinates**: Integer micrometers, 1 um precision.

10. **No external ML dependencies**: All statistical methods (logistic regression, PC algorithm, SEM) implemented from scratch using only numpy/scipy.

---

## 7. Performance Characteristics

| Domain | Cells | Speed | 168h Runtime |
|--------|-------|-------|--------------|
| 1mm x 1mm | ~5,000 | ~3 hr/s | ~1 min |
| 2mm x 2mm | ~20,000 | ~1 hr/s | ~3 min |

**Bottleneck**: Diffusion solver (100 sub-steps/hour x 4 substrates). With Numba: ~2x faster.

---

## 8. GBM Immunosuppression Rationale

The default immune parameters model glioblastoma's immunosuppressive tumor microenvironment:

| Mechanism | How Modeled |
|-----------|-------------|
| Blood-brain barrier | Low `recruitment_rate_per_hr` (0.003) |
| PD-L1 / TGF-beta suppression | Low `kill_probability_per_synapse` (0.25) |
| Rapid T cell exhaustion | `exhaustion_per_kill` = 0.35 (exhausted after ~3 kills) |
| Difficult activation | `min_activation_for_kill` = 0.5, slow activation rate |
| Fast deactivation | `deactivation_rate_per_hr` > `activation_rate_per_hr` |

Result: tumors grow despite immune pressure, matching clinical GBM observations.

---

## 9. Documentation System

Three documents, each in Markdown (source) and HTML (rendered):

1. **tutorial.md** — Teaches causal inference: correlation vs causation, IV methods, 2SLS, hands-on Python examples
2. **analysis_guide.md** — Statistical methods reference: OLS, 2SLS, sibling comparison, fixed effects, diagnostics, pitfalls
3. **immune_system.md** — Immune system module: synapse formation, exhaustion, GBM biology, parameter tuning

All HTML files link to `docs/assets/style.css` (Aptos/Calibri font, blue-themed headers, styled tables with alternating row colors, dark code blocks).

Rendering: `python scripts/render_docs.py` (requires `pip install markdown`).

---

## 10. Build Steps to Recreate from Scratch

### Phase 1: Foundation
1. Create `pyproject.toml` with dependencies and entry points
2. Implement `config.py` — frozen dataclass hierarchy, JSON loader
3. Implement `domain.py` — coordinate transforms, interpolation, gradients
4. Create `default.json` with all parameters

### Phase 2: Cell System
5. Implement `cells.py` — Cell dataclass, CellPopulation, spatial hash
6. Implement `shapes.py` — SVG path parsing, AABB collision detection
7. Implement `ecdna.py` — segregation, EGFR expression, four modulation functions

### Phase 3: Environment
8. Implement `environment.py` — EnvironmentFields, DiffusionSolver with Thomas algorithm
9. Implement `_numba_kernels.py` — optional JIT Thomas solver

### Phase 4: Behaviors
10. Implement `behaviors.py` — cell cycle, division, migration, death, immune killing, state transitions
11. Implement `angiogenesis.py` — VEGF-driven sprouting and vessel maturation
12. Implement `initialization.py` — vascular network, tissue population, burn-in, tumor seeding

### Phase 5: Simulation Loop
13. Implement `core.py` — 6-phase loop, CLI entry point
14. Implement `io.py` — TSV writers for cells, environment, lineage, summary

### Phase 6: Visualization
15. Implement `visualization.py` — Vega-Lite spec
16. Implement `viewer.py` — standalone HTML viewer
17. Implement `movie.py` — animated SVG movie
18. Implement `reporting.py` — HTML simulation report

### Phase 7: Analysis
19. Implement `loader.py` — SimulationData, TSV parsing
20. Implement `effects.py` — ground truth estimation
21. Implement `iv.py` — 2SLS
22. Implement `regression.py` — OLS, polynomial, interaction, doubly robust
23. Implement `discovery.py` — PC and GES algorithms
24. Implement `matching.py` — propensity scores, nearest-neighbor, IPW
25. Implement `sem.py` — path analysis, mediation, Sobel test
26. Implement `sensitivity.py` — Rosenbaum, E-value, OVB, placebo
27. Implement `report.py` — HTML analysis report
28. Implement `cli.py` — CLI for analysis
29. Implement `graph/causal_dag.py` — DAG construction and comparison

### Phase 8: Utilities and Docs
30. Implement `cleanup.py` — TSV reduction
31. Write `tutorial.md`, `analysis_guide.md`, `immune_system.md`
32. Create `docs/assets/style.css`
33. Write `scripts/render_docs.py`
34. Write `examples/analyze_simulation.py`
35. Write `README.md`

---

## 11. Testing a Rebuild

A correct rebuild should produce simulations where:

1. **Tumor grows** despite immune pressure (GBM-realistic)
2. **ecDNA segregation** has mean daughter fraction = 0.5 (binomial)
3. **2SLS recovers true effects** within ~10% error (OLS is biased by ~100%+)
4. **First-stage F-statistic** >> 10 (strong instrument)
5. **ecDNA does NOT correlate** with O2_local (instrument independence)
6. **Cells with more ecDNA**: divide faster, migrate faster, secrete more VEGF, resist apoptosis
7. **Hypoxic cells**: migrate faster (Go-or-Grow), don't proliferate
8. **168-hour simulation**: completes in ~1-4 minutes on 1mm domain

---

## 12. Summary of All Mathematical Formulas

| Formula | Module | Purpose |
|---------|--------|---------|
| `EGFR = 1.0 + 0.5*ecDNA + LogNormal(0,0.1)` | ecdna.py | Gene dosage expression |
| `T_div = T_base / (1 + 0.3*log2(1+EGFR))` | ecdna.py | Proliferation acceleration |
| `V_VEGF = V_base * (1 + 0.1*sqrt(EGFR))` | ecdna.py | Angiogenic signaling |
| `v_mig = v_base * (1 + 0.05*EGFR) * hypoxia_boost` | ecdna.py | Invasion capacity |
| `a_apo = a_base / (1 + 0.2*log2(1+EGFR))` | ecdna.py | Apoptosis resistance |
| `r = D*dt / (2*dx^2)` | environment.py | Diffusion number (LOD) |
| `-r*u_{i-1} + (1+2r)*u_i - r*u_{i+1} = u_i^n` | environment.py | Implicit x-sweep |
| `O2_src = q_O2 * vasc_dens * (O2_blood - local_O2)` | environment.py | Vascular O2 supply |
| `P_kill = kill_prob * synapse_progress * (1-exhaust) * activation * dt` | behaviors.py | Immune killing |
| `daughter_ecDNA ~ Binomial(replicated, p=0.5)` | ecdna.py | ecDNA segregation |
| `beta_IV = Cov(Y,Z) / Cov(X,Z)` | iv.py | 2SLS estimator |
| `F = (r^2/1) / ((1-r^2)/(n-2))` | iv.py | First-stage F-stat |
| `E = RR + sqrt(RR*(RR-1))` | sensitivity.py | E-value |
