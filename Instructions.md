# CAUSANTA: Causal Analysis Using Somatic And Neighborhood Tissue Architecture
 
## Theoretical Framework & Simulation Engine Specification
 
---
 
## 1. Motivation and Scientific Context
 
Spatial transcriptomics and proteomics capture high-dimensional tissue organization, but current computational frameworks are limited to descriptive analysis. Existing methods do not distinguish causal drivers from downstream effects, nor do they predict tissue responses to molecular perturbations. These limitations stem from the observational nature of spatial data, which is heavily confounded by microenvironmental "neighborhood effects" that mask true causal drivers. Fundamental questions remain unanswered: *Which computational approaches can reliably identify causal relationships? How can we reconstruct tissue evolution from static data? How can we predict tissue responses to molecular perturbations?*
 
**CAUSANTA addresses these gaps by treating somatic stochasticity as endogenous biological randomization, enabling causal and temporal inference from multimodal spatial data.** The key innovation is to treat heritable cell-to-cell genomic variation as Somatic Instrumental Variables (SIVs) within a Structural Causal Model. This enables identification of causal effects from otherwise observational spatial data. Extrachromosomal DNA (ecDNA) segregation serves as an ideal validation system: its randomized inheritance acts as a natural experiment, and its genomic ground truth enables rigorous benchmarking. While ecDNA anchors validation, the framework generalizes to other forms of heritable somatic variation, enabling causal inference across diverse diseases and biological contexts.
 
The CAUSANTA simulation engine generates synthetic tissue sections with known ground-truth causal structure, enabling rigorous benchmarking of causal discovery algorithms. By simulating evolutionary trajectories within tissues — treating somatic variation as a biological clock — CAUSANTA supports inference of temporal progression, identification of lineage-specific drivers of clonal expansion, and in silico molecular intervention experiments that explore how specific perturbations alter the tissue microenvironment.
 
---
 
## 2. Spatial Domain and Coordinate System
 
### 2.1 Domain
 
CAUSANTA simulates a 2D tissue section representing a histological slice. The domain size is configurable up to **12 mm × 12 mm** (12,000 μm × 12,000 μm).
 
### 2.2 Coordinate Precision
 
The **minimal spatial unit is 1 μm**. All cell positions (x, y) are recorded in micrometer coordinates at integer precision. This provides sub-cellular coordinate resolution for accurate spatial registration with FISH, PhenoCycler, and H&E imaging data.
 
### 2.3 Environment Grid
 
Environment fields (O₂, glucose, etc.) are computed on a **coarser regular grid** to maintain computational tractability. The environment grid resolution is a configurable parameter (`env_grid_um`, default 10 μm). For a 1 mm × 1 mm domain at 10 μm resolution, the environment grid is 100 × 100 = 10,000 voxels. For a 12 mm domain at 10 μm resolution, it is 1,200 × 1,200 = 1,440,000 voxels.
 
Cell agents exist in continuous (integer μm) space and sample the environment grid via bilinear interpolation at their position.
 
### 2.4 Time
 
The **time step is 1 hour**. All rate parameters are expressed per hour. Internal diffusion sub-stepping occurs within each 1-hour step as needed for numerical stability.
 
---
 
## 3. Cell Agent Model
 
### 3.1 Cell Types
 
| Type ID | Cell Type | Typical Diameter (μm) | Role |
|---------|-----------|----------------------|------|
| 0 | Neuron | 10–20 | Parenchyma; post-mitotic; displaced/killed by tumor |
| 1 | Astrocyte | 10–15 (soma) | Glial support; tiles in non-overlapping ~50–100 μm domains; reactive astrogliosis at tumor margin |
| 2 | Oligodendrocyte | 8–12 | Myelinating glia; displaced by tumor; some gliomas have oligo-lineage features |
| 3 | Microglia | 8–10 | Resident immune; 5–10% of glia; activates and migrates toward tumor |
| 4 | Endothelial | 10–15 | Vasculature; O₂/glucose source; angiogenesis |
| 5 | Pericyte | 5–8 | Wraps endothelial cells; blood-brain barrier integrity; vascular remodeling |
| 6 | Tumor (Glioma) | 15–25 | Proliferative, invasive, ecDNA-bearing |
| 7 | Recruited Immune | 10–12 | T cells, macrophages entering from vasculature |
| 8 | Necrotic | variable | Dead cell; releases DAMPs; occupies space until cleared |
 
### 3.2 Cell State Vector
 
Every cell agent c carries a complete state record at each time step:
 
```
cell_id           int       Unique identifier (monotonically increasing)
parent_id         int       cell_id of mitotic parent (-1 for initial cells)
cell_type         int       Type code (0–8, see 3.1)
cell_type_name    string    Human-readable name
 
# Spatial
x                 int       X position in μm (1 μm precision)
y                 int       Y position in μm
angle             float     Orientation angle in radians [0, 2π)
 
# Morphology
shape_path        string    SVG path string defining cell outline (Vega mark path format)
size_scale        float     Scale factor applied to shape_path (1.0 = default for type)
 
# Proliferation
cell_cycle_phase  string    {G0, G1, S, G2, M}
cycle_clock_hr    float     Hours elapsed in current phase
 
# Genomic / ecDNA
ecDNA_count       int       Copy number of extrachromosomal DNA elements
ecDNA_cargo       string    Semicolon-delimited amplified genes (e.g., "EGFR;MYC")
 
# Lineage
generation        int       Division count from founding ancestor
time_born_hr      float     Simulation hour when cell was created
 
# Local environment (sampled)
O2_local          float     Local oxygen (mmHg) at cell position
glucose_local     float     Local glucose (mM) at cell position
 
# Derived state
is_hypoxic        bool      O₂ < hypoxia threshold
is_quiescent      bool      In G0 (contact inhibited or resource-starved)
 
# Current rates (effective, after ecDNA modulation)
migration_rate    float     Effective migration speed (μm/hr)
VEGF_secretion    float     Effective VEGF output (amol/hr)
```
 
### 3.3 Cell Shape Representation
 
Cell shapes are defined as **SVG path strings**, following the same convention as Vega's `path` mark type. Each cell type has a default shape path normalized to a unit bounding box (roughly -1 to 1 in both axes). The `size_scale` factor scales the path uniformly, and `angle` rotates it.
 
**Default shape paths by cell type:**
 
```json
{
  "Neuron":          "M 0 -1 C 0.5 -0.5 0.3 0.2 0.8 1 L 0 0.6 L -0.8 1 C -0.3 0.2 -0.5 -0.5 0 -1 Z",
  "Astrocyte":       "M 0 -1 L 0.4 -0.3 L 1 -0.2 L 0.4 0.2 L 0.5 1 L 0 0.5 L -0.5 1 L -0.4 0.2 L -1 -0.2 L -0.4 -0.3 Z",
  "Oligodendrocyte": "M 0 -0.8 A 0.8 0.8 0 1 1 0 0.8 A 0.8 0.8 0 1 1 0 -0.8 Z",
  "Microglia":       "M 0 -0.6 L 0.3 -0.3 L 0.8 -0.5 L 0.5 0 L 0.8 0.5 L 0.3 0.3 L 0 0.6 L -0.3 0.3 L -0.8 0.5 L -0.5 0 L -0.8 -0.5 L -0.3 -0.3 Z",
  "Endothelial":     "M -1 -0.3 Q 0 -0.5 1 -0.3 L 1 0.3 Q 0 0.5 -1 0.3 Z",
  "Pericyte":        "M -0.6 -0.4 Q 0 -0.6 0.6 -0.4 L 0.6 0.4 Q 0 0.6 -0.6 0.4 Z",
  "Tumor":           "M 0 -0.9 C 0.6 -0.9 0.9 -0.4 0.9 0 C 0.9 0.5 0.5 0.9 0 0.9 C -0.5 0.9 -0.9 0.5 -0.9 0 C -0.9 -0.4 -0.6 -0.9 0 -0.9 Z",
  "RecruitedImmune": "M 0 -0.7 A 0.7 0.7 0 1 1 0 0.7 A 0.7 0.7 0 1 1 0 -0.7 Z",
  "Necrotic":        "M -0.5 -0.5 L 0.5 -0.5 L 0.5 0.5 L -0.5 0.5 Z"
}
```
 
The rendering pipeline (Vega.js) transforms each cell's path by: scale by `size_scale × (diameter / 2)`, rotate by `angle`, translate to `(x, y)`.
 
### 3.4 Cell Type Parameters
 
Each cell type τ carries a parameter set **P_τ**. All rates are **per hour** (the simulation time unit).
 
```
# Identity
type_id                     int
type_name                   string
shape_path                  string     SVG path for default shape
color                       string     Hex color for visualization
 
# Geometry
diameter_mean_um            float      Mean cell body diameter (μm)
diameter_std_um             float      Standard deviation (μm)
 
# Proliferation
can_divide                  bool       Whether this type can undergo mitosis
division_time_mean_hr       float      Mean total cell cycle duration (hours)
division_time_std_hr        float      Stochastic variation (hours)
O2_prolif_threshold_mmHg    float      Minimum O₂ for proliferation
glucose_prolif_threshold_mM float      Minimum glucose for proliferation
max_generations             int        Division limit (-1 = unlimited)
 
# Migration
migration_speed_um_hr       float      Base migration speed (μm per hour)
migration_persistence_hr    float      Directional memory time constant (hours)
chemotaxis_O2               float      Coefficient: bias toward O₂ gradient
chemotaxis_VEGF             float      Coefficient: bias toward VEGF gradient
chemotaxis_chemokine        float      Coefficient: bias toward chemokine gradient
haptotaxis_ECM              float      Coefficient: bias along ECM gradient
contact_inhibited           bool       Whether migration stops at confluent neighbors
 
# Death
O2_necrosis_threshold_mmHg  float      O₂ below which necrosis initiates
necrosis_delay_hr           float      Hours below threshold before death
apoptosis_rate_per_hr       float      Spontaneous death probability per hour
lysis_rate_per_hr           float      Rate at which necrotic cells are cleared
 
# Metabolism & Secretion
O2_consumption_amol_hr      float      Oxygen consumption (attomoles per hour)
glucose_consumption_amol_hr float      Glucose consumption (attomoles per hour)
VEGF_secretion_amol_hr      float      VEGF production when hypoxic
lactate_production_amol_hr  float      Lactate output (Warburg effect)
 
# State transition probabilities (per hour, conditional on trigger)
transition_rules            list       Each entry: {condition, target_type, rate_per_hr}
 
# ecDNA (tumor cells only)
ecDNA_init_count            int        Starting ecDNA copy number
ecDNA_cargo                 list       Gene names on ecDNA
ecDNA_segregation_p         float      Per-copy probability of going to daughter (0.5 = random)
ecDNA_effect_on_division    float      α: proliferation acceleration per log₂(ecDNA)
ecDNA_effect_on_VEGF        float      β: VEGF secretion multiplier per ecDNA copy
ecDNA_effect_on_migration   float      δ: migration speed multiplier per ecDNA copy
ecDNA_effect_on_survival    float      γ: apoptosis resistance per ecDNA copy
```
 
#### Reference Values
 
| Parameter | Neuron | Astrocyte | Oligo | Microglia | Endothelial | Pericyte | Tumor | Recruited Immune |
|-----------|--------|-----------|-------|-----------|-------------|----------|-------|-----------------|
| diameter (μm) | 15 ± 3 | 12 ± 2 | 10 ± 2 | 9 ± 1 | 12 ± 2 | 6 ± 1 | 18 ± 3 | 11 ± 1 |
| can_divide | no | yes | yes | yes | yes | yes | yes | yes |
| division_time (hr) | — | 168 ± 24 | 168 ± 24 | 72 ± 12 | 60 ± 12 | 96 ± 12 | 36 ± 8 | 36 ± 8 |
| migration (μm/hr) | 0 | 3 ± 1 | 1 ± 0.5 | 30 ± 10 | 10 ± 3 | 5 ± 2 | 10 ± 3 | 35 ± 10 |
| O₂ consumption (amol/hr) | 30000 | 12000 | 9000 | 6000 | 3000 | 2000 | 72000 | 12000 |
| O₂ prolif threshold (mmHg) | — | 10 | 10 | 5 | 8 | 8 | 8 | 5 |
| O₂ necrosis threshold (mmHg) | 2.5 | 2.5 | 2.5 | 1.0 | 2.5 | 2.5 | 2.5 | 2.5 |
| VEGF secretion (amol/hr, hypoxic) | 0 | 30 | 0 | 0 | 0 | 0 | 600 | 0 |
| apoptosis rate (/hr) | 1e-4 | 1e-4 | 1e-4 | 1e-4 | 1e-4 | 1e-4 | 5e-5 | 1e-3 |
| color | #4477AA | #66CCEE | #228833 | #EE6677 | #AA3377 | #BBBBBB | #CCBB44 | #EE8866 |
 
---
 
## 4. Environment Field Model
 
### 4.1 Substrates and Governing Equation
 
Each environment field E_k is defined on the environment grid (resolution `env_grid_um`) and evolves via a reaction-diffusion partial differential equation:
 
```
∂E_k/∂t = D_k ∇²E_k  −  λ_k E_k  +  sources(x,t)  −  sinks(x,t)
```
 
where:
- `D_k` = diffusion coefficient (μm²/hr)
- `λ_k` = natural decay rate (1/hr)
- `sources` = supply from vasculature and cell secretion
- `sinks` = cellular consumption
 
### 4.2 Primary Fields
 
| Field | Symbol | Units | D_k (μm²/hr) | λ_k (/hr) | Normal Value | Role |
|-------|--------|-------|---------------|-----------|-------------|------|
| Oxygen | O₂ | mmHg | 6.0 × 10⁶ | ~0 | 38 mmHg (~5%) | Proliferation, hypoxia, necrosis, HIF |
| Glucose | Gluc | mM | 2.4 × 10⁶ | ~0 | 5.0 mM | Metabolic fuel; Warburg effect |
| VEGF | VEGF | nM | 3.6 × 10⁴ | 0.6 | 0 nM | Angiogenesis trigger |
| Lactate | Lac | mM | 9.0 × 10⁵ | 0.06 | 1.0 mM | pH proxy; immunosuppressive |
| ECM density | ECM | [0–1] | 0 (static) | 0 | 0.5 | Migration scaffold; tumor remodeling |
| Vascular density | Vasc | [0–1] | 0 (static) | 0 | varies | O₂/glucose source; immune entry |
 
Note: Diffusion coefficients converted to μm²/hr from standard μm²/min values (×60).
 
### 4.3 Environment Grid Record
 
Each environment grid voxel at position (i, j) stores:
 
```
x_grid            int       Grid column index
y_grid            int       Grid row index
x_um              float     Physical X position (μm)
y_um              float     Physical Y position (μm)
O2                float     Oxygen partial pressure (mmHg) [0–60]
glucose           float     Glucose concentration (mM) [0–5.5]
VEGF              float     VEGF concentration (nM) [0–100]
lactate           float     Lactate concentration (mM) [0–40]
ECM_density       float     Extracellular matrix density [0–1]
vascular_density  float     Local vessel density [0–1]
pH                float     Derived: 7.4 − 0.02 × lactate
is_occupied       bool      Whether a cell occupies this voxel center
occupant_cell_id  int       cell_id of occupant (-1 if empty)
occupant_type     int       Cell type of occupant (-1 if empty)
```
 
### 4.4 Source and Sink Terms
 
**Oxygen supply from vasculature:**
```
source_O2(x) = q_O2 × vasc_density(x) × (O2_blood − O2(x))
```
where `O2_blood` = 60 mmHg, `q_O2` is a transfer coefficient (hr⁻¹).
 
**Oxygen consumption by cells:**
```
sink_O2(x) = Σ_c  O2_consumption_c × δ(x − x_c) / V_voxel
```
where the sum runs over cells c located within the voxel, and V_voxel = (env_grid_um)² is the voxel area.
 
**VEGF secretion (hypoxia-dependent):**
```
source_VEGF(x) = Σ_c  VEGF_rate_c × H(θ_hypoxia − O2(x_c)) × δ(x − x_c) / V_voxel
```
where H is the Heaviside function.
 
**Glucose** follows the same source/sink pattern as oxygen. **Lactate** is produced by tumor cells proportional to glucose consumption (Warburg effect).
 
### 4.5 Boundary Conditions
 
Dirichlet (fixed value) at domain edges, representing surrounding normal tissue:
- O₂ = 38 mmHg
- Glucose = 5.0 mM
- VEGF = 0 nM
- Lactate = 1.0 mM
 
### 4.6 Numerical Diffusion Scheme
 
**Implicit Locally One-Dimensional (LOD) operator splitting**, solving each spatial direction via the Thomas algorithm (tridiagonal matrix solve):
 
For each substrate, within each 1-hour time step, diffusion is sub-stepped at `dt_diffusion` (configurable, default 6 seconds = 0.1 min):
 
1. X-sweep: solve tridiagonal system implicit in x for each row
2. Y-sweep: solve tridiagonal system implicit in y for each column
3. Apply source/sink terms as a separate operator step
4. Enforce boundary conditions
5. Clamp values to physical ranges
 
The implicit scheme is **unconditionally stable** — no CFL constraint on `dt_diffusion`, allowing large sub-steps even with oxygen's high diffusion coefficient.
 
---
 
## 5. Cell Behavioral Rules
 
All rates in this section are **per hour** unless otherwise noted.
 
### 5.1 Proliferation (Mitosis)
 
For each cell c whose type permits division, at each time step:
 
**1. Resource gate:**
```
can_proliferate = (O2_local > O2_prolif_threshold)
                  AND (glucose_local > glucose_prolif_threshold)
                  AND (generation < max_generations OR max_generations == -1)
                  AND (at least one adjacent position unoccupied)
```
 
**2. Cell cycle progression:**
 
Phase durations are fractions of the total division time T (drawn once per cycle from Normal(mean, std)):
- G1: 0.40 × T
- S:  0.30 × T
- G2: 0.15 × T
- M:  0.15 × T
 
`cycle_clock_hr` advances by 1 (the time step) each hour. Phase transitions occur when accumulated time exceeds phase duration. If resources drop below threshold during S/G2/M, the cell arrests (clock pauses, does not revert).
 
**3. Division (M phase completion):**
 
- Select a random unoccupied position adjacent to the parent (8-connected neighborhood)
- If no position available → cell enters G0 (contact inhibition), clock pauses
- If position available:
  - Create daughter cell with new `cell_id`, `parent_id` = parent's `cell_id`
  - Daughter placed at selected position
  - Both cells reset `cycle_clock_hr = 0`, enter G1
  - `generation` incremented for both cells
 
**ecDNA segregation (the SIV mechanism):**
```
For parent with ecDNA_count = N:
    daughter_ecDNA = Binomial(N, p)    where p = ecDNA_segregation_p (default 0.5)
    parent_ecDNA   = N − daughter_ecDNA
```
 
This binomial partitioning is the core mechanism that generates inter-cell ecDNA heterogeneity — creating the "natural experiment" that enables causal inference. The randomness of segregation satisfies the independence assumption required for instrumental variable analysis.
 
**ecDNA modulation of cell behavior:**
```
division_time_effective = division_time_base / (1 + α × log₂(1 + ecDNA_count))
VEGF_secretion_effective = VEGF_secretion_base × (1 + β × ecDNA_count)
migration_effective = migration_base × (1 + δ × ecDNA_count)
apoptosis_effective = apoptosis_base / (1 + γ × ecDNA_count)
```
 
These are the **ground-truth causal effects** that CAUSANTA embeds in the simulation. The causal discovery algorithms being benchmarked should recover these relationships from the simulated spatial data.
 
### 5.2 Migration
 
Each hour, motile cells compute a migration vector and move:
 
**1. Direction:**
```
v = v_persist + χ_O2 ∇O₂ + χ_VEGF ∇VEGF + χ_chemo ∇chemokine + χ_ECM ∇ECM + v_noise
```
 
- `v_persist` = previous direction × exp(−1/persistence_time), providing directional memory
- Gradient terms: environment fields sampled at cell position, finite-difference gradient over neighboring voxels
- `v_noise` = uniform random angle perturbation
 
**2. Displacement:**
```
Δx = migration_speed × cos(θ_v) × 1 hr
Δy = migration_speed × sin(θ_v) × 1 hr
```
 
Position updated: `x += round(Δx)`, `y += round(Δy)` (integer μm).
 
**3. Collision:** If target position is within one cell diameter of another cell, the move is rejected (volume exclusion). The cell may attempt an alternate direction or stay.
 
**Type-specific behaviors:**
- **Neurons, Oligodendrocytes:** migration_speed = 0 (sessile)
- **Astrocytes:** Low motility; increased near tumor margin (reactive astrogliosis)
- **Microglia:** High motility; strong chemotaxis toward DAMPs from necrotic regions
- **Tumor:** Moderate motility; biased toward O₂ gradients, along ECM fibers; optional Go-or-Grow coupling: `speed_effective = speed_base × (1 − proliferation_activity)`
- **Recruited Immune:** Highest motility; chemotaxis toward tumor-secreted chemokines
- **Pericytes:** Low motility; constrained to remain adjacent to endothelial cells
 
### 5.3 Cell Death
 
**Necrosis:** When `O2_local < O2_necrosis_threshold` for `necrosis_delay_hr` consecutive hours:
- `cell_type` → 8 (Necrotic)
- Cell remains at position, blocks space, ceases metabolism
- Cleared stochastically: `P(lysis) = lysis_rate × 1 hr` per time step
 
**Apoptosis:** At each step:
- `P(apoptosis) = apoptosis_rate_per_hr × 1 hr`
- Apoptotic cells removed from grid immediately
 
**Immune-mediated kill:** When a recruited immune cell or activated microglia is within kill_radius of a tumor cell:
- `P(kill) = kill_rate × activation_level × 1 hr`
 
### 5.4 Cell Type Transitions
 
| From | Trigger Condition | To | Rate (per hr) |
|------|------------------|----|---------------|
| Astrocyte | Adjacent to tumor cell | Reactive Astrocyte (subtype flag) | 0.01 |
| Microglia | Local chemokine > threshold | Activated Microglia (subtype flag) | 0.05 |
| Activated Microglia | Local TGF-β > threshold | Immunosuppressed Microglia | 0.02 |
| Endothelial | Local VEGF > angio_threshold | Sprouting Endothelial (subtype flag) | 0.01 |
| Any living cell | O₂ < necrosis threshold (sustained) | Necrotic | deterministic |
| (vascular source) | Local chemokine > recruit_threshold | Recruited Immune (new cell spawned) | 0.005 |
 
### 5.5 Angiogenesis
 
When VEGF exceeds `angiogenesis_threshold_nM` at a position adjacent to an existing endothelial cell:
 
1. **Tip cell selection:** endothelial cell with highest local VEGF
2. **Sprouting:** tip cell extends one position per hour up the VEGF gradient; a new endothelial cell is placed at the vacated position (stalk cell)
3. **Anastomosis:** if tip reaches another vessel, vessels connect
4. **Maturation:** new vessel positions gradually increase `vascular_density` over `vessel_maturation_hr` hours
5. **Effect:** new vascular positions become O₂/glucose sources
 
---
 
## 6. ecDNA as Somatic Instrumental Variable (SIV)
 
### 6.1 The Causal Inference Problem
 
In observational spatial data, associations between molecular features (e.g., gene expression) and cellular phenotypes (e.g., proliferation) are confounded by shared microenvironment. A cell near a vessel has both higher O₂ *and* different gene expression; separating cause from correlation is impossible without randomization.
 
### 6.2 ecDNA Segregation as Natural Randomization
 
ecDNA elements segregate during mitosis via a mechanism that is:
 
- **Random:** Each ecDNA copy is independently assigned to a daughter cell with probability p ≈ 0.5 (binomial process)
- **Heritable:** Daughter cells inherit the randomly assigned count, which then influences their phenotype
- **Independent of confounders:** The physical segregation mechanism is not driven by the local microenvironment
 
This satisfies the three instrumental variable assumptions:
1. **Relevance:** ecDNA copy number affects the instrumented variable (e.g., EGFR expression)
2. **Independence:** ecDNA segregation is independent of confounders (microenvironment)
3. **Exclusion restriction:** ecDNA affects outcome only through the instrumented pathway
 
### 6.3 Ground-Truth Causal Graph
 
CAUSANTA embeds a known Structural Causal Model:
 
```
ecDNA_count ──→ EGFR_expression ──→ proliferation_rate
     │                                      ↑
     │                              O2_local (confounder)
     │                                      ↑
     └──→ VEGF_secretion ──→ angiogenesis ──┘
```
 
The simulation knows the true causal effect sizes (α, β, γ, δ from Section 5.1). Causal discovery algorithms operating on the output spatial data should recover this structure. The degree to which they succeed — measured against these known ground-truth edges — constitutes the benchmark.
 
### 6.4 Counterfactual Trajectories
 
CAUSANTA supports **intervention experiments**: at any time point, the user can modify ecDNA status at specific positions and re-run the simulation forward to observe how the tissue evolves differently. Comparing the factual trajectory (observed ecDNA distribution) against counterfactual trajectories (modified ecDNA) quantifies causal effects spatially.
 
---
 
## 7. Initialization
 
### 7.1 Synthetic Initialization (Default Mode)
 
Normal brain tissue is populated according to biologically realistic proportions for cortical gray matter:
 
| Cell Type | Fraction | Placement Rule |
|-----------|----------|---------------|
| Neuron | 40–50% | Distributed across parenchyma |
| Astrocyte | 20–30% | Tiled in non-overlapping ~50–100 μm domains |
| Oligodendrocyte | 10–15% | Interspersed, denser near white matter tracts |
| Microglia | 5–10% | Scattered uniformly; 5–10% of all glial cells |
| Endothelial | ~3–5% | Organized as capillary network segments |
| Pericyte | ~1–2% | Adjacent to endothelial cells |
 
**Vascular network:** Placed first as a branching capillary pattern with inter-capillary distance of ~100–200 μm. Pericytes assigned to positions adjacent to endothelial cells.
 
**Tumor seeding:** One or more tumor cells placed at specified coordinates with initial ecDNA_count, ecDNA_cargo, and cell cycle phase.
 
**Environment equilibration:** After placing cells and vasculature, the diffusion solver runs for a burn-in period (default 100 hours equivalent) to reach steady-state O₂/glucose fields before the first recorded time step.
 
### 7.2 Data-Driven Initialization (for Real Data)
 
Load real spatial data and register to the CAUSANTA coordinate grid:
 
- **H&E:** Tissue architecture → ECM_density map, cell segmentation → positions and morphology
- **FISH:** Per-cell ecDNA copy number → `ecDNA_count` per cell
- **PhenoCycler / CODEX:** Protein panel → cell type classification, positions, protein expression
- **Visium:** 55 μm spots deconvolved to cell neighborhoods → transcriptomic features per cell
- **AlphaCycler:** Additional spatial proteomic channels
 
Environment fields initialized from tissue context (vessel positions from CD31 staining, ECM from collagen markers) rather than synthetic baseline.
 
---
 
## 8. Simulation Main Loop
 
```
INITIALIZE:
    1. Load parameters from JSON configuration file
    2. Create environment grid at specified resolution
    3. Initialize environment fields to baseline values
    4. Place vascular network (synthetic or from data)
    5. Populate with normal cell types at specified proportions
    6. Seed tumor cell(s) at specified position(s)
    7. Equilibrate environment fields (burn-in diffusion)
    8. Write initial state: cells_t000000.tsv, environment_t000000.tsv
    9. Copy parameter JSON to output directory
 
FOR each time step t = 1, 2, ..., T_total:
    
    PHASE 1: ENVIRONMENT DIFFUSION
        For each diffusion sub-step within the 1-hour main step:
            For each substrate (O₂, glucose, VEGF, lactate):
                a. Compute source terms (vascular supply, cell secretion)
                b. Compute sink terms (cell consumption)
                c. Solve diffusion via implicit LOD (Thomas algorithm)
                d. Apply natural decay (λ_k)
                e. Enforce boundary conditions and physical range clamps
        Compute derived fields (pH from lactate, hypoxia flags)
    
    PHASE 2: CELL BEHAVIOR (cells processed in random-shuffled order)
        For each living cell c:
            a. Sample local environment via bilinear interpolation
            b. Update is_hypoxic, is_quiescent flags
            c. Death check:
                - Necrosis: O₂ below threshold for required duration?
                - Apoptosis: stochastic roll
            d. State transition check (Table 5.4)
            e. Proliferation:
                - Advance cell cycle clock
                - If M-phase complete and space available:
                    Execute division, partition ecDNA binomially
                    Create daughter with new cell_id, parent_id
            f. Migration:
                - Compute direction (gradients + persistence + noise)
                - Compute displacement (speed × 1 hr)
                - Attempt move (volume exclusion check)
    
    PHASE 3: ANGIOGENESIS
        For each endothelial cell with adjacent VEGF > threshold:
            Attempt sprouting (Section 5.5)
    
    PHASE 4: IMMUNE RECRUITMENT
        For each vascular position with chemokine > threshold:
            Stochastic spawn of recruited immune cell
    
    PHASE 5: CLEANUP
        Stochastic lysis of necrotic cells
        Update environment grid occupancy
    
    PHASE 6: OUTPUT (at configured interval)
        Write cells_t{step:06d}.tsv
        Write environment_t{step:06d}.tsv
        Write summary statistics to log
```
 
---
 
## 9. Output Specification
 
### 9.1 Directory Structure
 
Each simulation run produces a timestamped output directory:
 
```
output/
  run_20260410_143022/
    params.json                    # Complete parameter file (frozen copy)
    cells_t000000.tsv              # Initial cell state
    cells_t000001.tsv              # t = 1 hour
    cells_t000002.tsv              # t = 2 hours
    ...
    environment_t000000.tsv        # Initial environment
    environment_t000001.tsv
    ...
    summary.log                    # Per-step statistics (cell counts, mean O₂, etc.)
    lineage.tsv                    # Complete parent-child lineage table
```
 
### 9.2 Cells File Format (TSV)
 
**Filename:** `cells_t{step:06d}.tsv`
 
| Column | Type | Description |
|--------|------|-------------|
| cell_id | int | Unique cell identifier |
| parent_id | int | Parent cell ID (-1 for initial cells) |
| cell_type | int | Type code (0–8) |
| cell_type_name | string | Human-readable type |
| x | int | X position in μm |
| y | int | Y position in μm |
| angle | float | Orientation (radians) |
| shape_path | string | SVG path string |
| size_scale | float | Scale factor for shape rendering |
| cell_cycle_phase | string | G0, G1, S, G2, or M |
| ecDNA_count | int | ecDNA copy number |
| ecDNA_cargo | string | Semicolon-delimited gene list |
| O2_local | float | Local oxygen (mmHg) |
| glucose_local | float | Local glucose (mM) |
| is_hypoxic | bool | True/False |
| is_quiescent | bool | True/False |
| generation | int | Divisions from ancestor |
| time_born_hr | float | Hour when cell was created |
| migration_rate | float | Current effective speed (μm/hr) |
| VEGF_secretion | float | Current VEGF output (amol/hr) |
 
### 9.3 Environment File Format (TSV)
 
**Filename:** `environment_t{step:06d}.tsv`
 
| Column | Type | Description |
|--------|------|-------------|
| x_grid | int | Grid column index |
| y_grid | int | Grid row index |
| x_um | float | X in μm |
| y_um | float | Y in μm |
| O2 | float | Oxygen (mmHg) |
| glucose | float | Glucose (mM) |
| VEGF | float | VEGF (nM) |
| lactate | float | Lactate (mM) |
| pH | float | Derived pH |
| ECM_density | float | ECM [0–1] |
| vascular_density | float | Vessel density [0–1] |
| is_occupied | bool | Cell present? |
| occupant_cell_id | int | cell_id or -1 |
| occupant_type | int | Cell type or -1 |
 
### 9.4 Lineage File
 
**Filename:** `lineage.tsv`
 
Append-only file recording every division event:
 
| Column | Type | Description |
|--------|------|-------------|
| time_hr | float | Hour of division |
| parent_id | int | Dividing cell's ID |
| parent_ecDNA_before | int | Parent ecDNA count before division |
| parent_ecDNA_after | int | Parent ecDNA count after segregation |
| daughter_id | int | New daughter cell's ID |
| daughter_ecDNA | int | Daughter ecDNA count |
| x | int | Position of division (μm) |
| y | int | Position of division (μm) |
| generation | int | Generation number of daughter |
 
---
 
## 10. Visualization (Vega.js)
 
### 10.1 Rendering Pipeline
 
The cells TSV is loaded directly into a Vega specification where each cell is rendered as a `path` mark:
 
1. Read `shape_path` as the Vega path datum
2. Apply transform: scale by `size_scale × (diameter/2)`, rotate by `angle`, translate to `(x, y)`
3. Color by `cell_type` using the type color palette
4. Opacity modulated by state (e.g., necrotic cells at 0.3 opacity, hypoxic at 0.7)
 
### 10.2 Environment Overlay
 
The environment TSV renders as a `rect` grid heatmap behind the cell layer:
- Channel selectable: O₂, glucose, VEGF, lactate, pH, ECM, vascular density
- Color scale: sequential (e.g., viridis for O₂, magma for VEGF)
- Opacity adjustable to see cells through environment
 
### 10.3 Interactive Features
 
- Time slider: scrub through simulation steps
- Cell selection: click cell to see full state vector, lineage tree
- Layer toggles: show/hide cell types, environment channels
- Zoom: pan and zoom across the tissue section
- ecDNA highlight mode: cells colored by ecDNA_count (continuous scale)
 
---
 
## 11. Core Mathematical Relationships
 
| Process | Equation | Key Parameters |
|---------|----------|---------------|
| Diffusion | ∂E/∂t = D∇²E − λE + S | D, λ per substrate |
| O₂ vascular supply | S = q × Vasc × (O₂_blood − O₂) | q (transfer coeff), O₂_blood = 60 mmHg |
| O₂-dependent proliferation | r = r_max × max(0, (O₂ − θ_H)/(O₂_max − θ_H)) | r_max, θ_H |
| ecDNA segregation | daughter_ecDNA ~ Binom(N, p) | N = parent count, p = 0.5 |
| ecDNA → division | T_eff = T_base / (1 + α log₂(1 + ecDNA)) | α (effect size) |
| ecDNA → VEGF | S_eff = S_base × (1 + β × ecDNA) | β |
| ecDNA → apoptosis resistance | a_eff = a_base / (1 + γ × ecDNA) | γ |
| ecDNA → migration | v_eff = v_base × (1 + δ × ecDNA) | δ |
| VEGF secretion | S_VEGF = r × H(θ_hyp − O₂) | r, θ_hyp |
| Migration | v = v_persist + χ∇O₂ + χ_ECM∇ECM + noise | χ coefficients |
| Necrosis | type → Necrotic if O₂ < θ_nec for > T_nec hours | θ_nec, T_nec |
| Immune kill | P(kill) = k × activation × Δt | k (kill rate) |
| pH from lactate | pH = 7.4 − 0.02 × [lactate] | empirical approximation |
 
---
 
## 12. Validation Targets
 
The simulation should reproduce these known biological phenomena:
 
1. **Avascular spheroid zonation:** Proliferating rim (~100–200 μm from vessels), quiescent intermediate zone, necrotic core — matching in vitro spheroid data
2. **Oxygen gradient:** Hypoxia onset at ~100–150 μm from nearest vessel
3. **Pseudopalisading necrosis:** Migrating waves of hypoxic cells around necrotic foci (GBM pathognomonic feature)
4. **ecDNA heterogeneity dynamics:** Inter-cell variance in ecDNA copy number increases with division count, following Var = N × p × (1−p) per generation under binomial segregation
5. **Angiogenic switch:** VEGF accumulation triggers vessel sprouting as tumor mass exceeds ~1–2 mm
6. **Immune exclusion:** High tumor density + immunosuppressive signaling creates immune-cold core with immune cells restricted to tumor periphery
7. **Gompertzian growth:** Tumor growth curve transitions from exponential → linear → plateau as carrying capacity limits and necrosis increase
8. **Causal recoverability:** Causal discovery algorithms applied to output data should recover the ground-truth ecDNA → phenotype causal edges (Section 6.3) at rates significantly above chance
 
---
 
## 13. Implementation Architecture
 
```
causanta/
    __init__.py
    config.py              # JSON parameter loading and validation
    domain.py              # Grid, coordinate system, boundary conditions
    cells.py               # Cell agent class, state vector, type parameters
    environment.py         # Environment fields, diffusion solver (LOD/Thomas)
    behaviors.py           # Proliferation, migration, death, transitions
    angiogenesis.py        # Vascular sprouting logic
    ecdna.py               # ecDNA segregation, modulation effects
    initialization.py      # Synthetic and data-driven tissue setup
    simulation.py          # Main loop orchestration
    io.py                  # TSV/JSON reading and writing
    visualization.py       # Vega.js spec generation
    README.md              # Project documentation
    memory.md              # Development log and decisions
    params/
        default.json       # Default parameter configuration
    output/                # Simulation run outputs (timestamped subdirs)
```
 
**Performance strategy:** Core simulation loop in Python (NumPy for environment field operations). If profiling reveals bottlenecks, migrate diffusion solver and cell iteration to Cython, Rust (via PyO3), or C extensions. The 100×100 environment grid with ~5,000–50,000 cells is well within pure-Python/NumPy performance for exploratory work.