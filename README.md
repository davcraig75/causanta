# CAUSANTA

**Causal Analysis Using Somatic And Neighborhood Tissue Architecture**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Executive Summary

CAUSANTA is a **simulation framework for benchmarking causal inference methods** in spatial biology. It generates synthetic tumor tissue with **known ground-truth causal structure**, enabling rigorous validation of causal discovery algorithms.

**The key innovation**: Extrachromosomal DNA (ecDNA) segregation acts as a **natural randomization mechanism**, satisfying instrumental variable (IV) assumptions and allowing causal effect identification from observational spatial data.

**Why this matters**: In biology, we rarely have both (1) ground truth causal effects AND (2) realistic confounding. CAUSANTA provides both, letting you test whether methods like 2SLS, propensity matching, or causal discovery algorithms actually recover true effects.

---

## Table of Contents

1. [Scientific Background](#scientific-background)
2. [The Core Innovation: Somatic Instrumental Variables](#the-core-innovation-somatic-instrumental-variables)
3. [Architecture Overview](#architecture-overview)
4. [Module Reference](#module-reference)
5. [Data Flow](#data-flow)
6. [Ground Truth Parameters](#ground-truth-parameters)
7. [Quick Start](#quick-start)
8. [Analysis Capabilities](#analysis-capabilities)
9. [Documentation](#documentation)
10. [Scope and Limitations](#scope-and-limitations)
11. [Development Notes](#development-notes)
12. [Citation](#citation)

---

## Scientific Background

### The Problem: Causal Inference in Biology

Causal inference in biology is fundamentally hard because:

1. **Randomized experiments are often infeasible** - You can't randomly assign cancer mutations to patients
2. **Observational data has confounding** - Hypoxia affects both gene expression AND cell behavior
3. **Ground truth is rarely known** - You don't know the true causal effect to validate against

### ecDNA Biology

**Extrachromosomal DNA (ecDNA)** are circular DNA molecules found in ~30% of cancers that:

- **Lack centromeres** → Cannot attach to spindle fibers during mitosis
- **Segregate randomly** → Each daughter cell gets ~Binomial(N, 0.5) copies
- **Carry oncogenes** → Often EGFR, MYC, or other drivers
- **Create heterogeneity** → Sisters can have very different copy numbers

This random segregation is the biological foundation for our instrumental variable approach.

### The Confounding Problem

In tumor tissue, **hypoxia (low oxygen)** creates confounding:

```
                 O2 (Hypoxia)
                ↙            ↘
      EGFR expression  →  Cell Proliferation
```

Hypoxic cells upregulate EGFR AND proliferate faster. Naive regression (OLS) conflates the direct EGFR→Proliferation effect with the confounded path through hypoxia.

---

## The Core Innovation: Somatic Instrumental Variables

### The SIV Framework

We use ecDNA copy number as an **instrumental variable** to identify causal effects:

```
    ecDNA (Instrument Z)
          ↓
    EGFR expression (Treatment D)  ←--  O2 (Confounder U)
          ↓                              ↓
    Proliferation (Outcome Y)  ←---------┘
```

### IV Assumptions (and how ecDNA satisfies them)

| Assumption | Requirement | How ecDNA Satisfies It |
|------------|-------------|------------------------|
| **Relevance** | Z affects D | ecDNA carries EGFR → gene dosage effect |
| **Independence** | Z ⊥ U | Binomial(N, 0.5) segregation is random |
| **Exclusion** | Z → Y only through D | ecDNA affects phenotypes via gene expression |

### Why This Works

During cell division:
1. ecDNA replicates (approximately doubles)
2. ecDNA segregates randomly: `N_daughter ~ Binomial(N_parent, 0.5)`
3. This creates **within-lineage variation** independent of environment
4. Sister cells in the same microenvironment have different ecDNA → different EGFR
5. This mimics **randomization within strata**

### Two-Stage Least Squares (2SLS)

**Stage 1** (First Stage): Predict treatment from instrument
```
EGFR_i = π₀ + π₁·ecDNA_i + η_i
```

**Stage 2** (Second Stage): Regress outcome on predicted treatment
```
Proliferation_i = β₀ + β₁·EGFR̂_i + ε_i
```

The coefficient β₁ is the **causal effect** of EGFR on proliferation, purged of confounding.

---

## Architecture Overview

### System Design

CAUSANTA has three main subsystems:

```
┌─────────────────────────────────────────────────────────────────┐
│                        CAUSANTA                                  │
├─────────────────┬─────────────────────┬─────────────────────────┤
│   SIMULATE      │      ANALYZE        │      VISUALIZE          │
├─────────────────┼─────────────────────┼─────────────────────────┤
│ Agent-based     │ IV/2SLS estimation  │ Nature Methods style    │
│ tumor growth    │ Bootstrap CI        │ Publication figures     │
│ ecDNA dynamics  │ Power analysis      │ Interactive viewers     │
│ Environment PDE │ Heterogeneity       │                         │
│ Immune system   │ Sensitivity         │                         │
└─────────────────┴─────────────────────┴─────────────────────────┘
```

### Simulation Architecture

The simulation runs a **6-phase loop** each timestep (0.5 hours):

```
┌──────────────────────────────────────────────────────────────┐
│                     SIMULATION LOOP                          │
├──────────────────────────────────────────────────────────────┤
│  1. ENVIRONMENT DIFFUSION                                    │
│     - Solve reaction-diffusion PDEs for O2, glucose, VEGF   │
│     - Vessels are sources, cells are sinks                  │
│                                                              │
│  2. CELL BEHAVIORS (random order)                           │
│     For each cell:                                           │
│     - Sample local environment                               │
│     - Update effective rates (ecDNA modulation)             │
│     - Check division → ecDNA segregation (Binomial)         │
│     - Check migration → move up gradients                   │
│     - Check apoptosis/necrosis                               │
│                                                              │
│  3. ANGIOGENESIS                                             │
│     - VEGF-driven vessel sprouting                          │
│                                                              │
│  4. IMMUNE RECRUITMENT                                       │
│     - T-cell recruitment based on tumor burden              │
│     - Immune killing of tumor cells                          │
│                                                              │
│  5. CLEANUP                                                  │
│     - Remove dead cells                                      │
│     - Clear necrotic debris                                  │
│                                                              │
│  6. OUTPUT                                                   │
│     - Write cell states, environment, lineage               │
└──────────────────────────────────────────────────────────────┘
```

### Cell Types

| Type | Division | Migration | ecDNA | Role |
|------|----------|-----------|-------|------|
| Tumor | Yes (36h) | Yes (10 μm/hr) | Yes | Primary subject |
| Neuron | No | No | No | Normal tissue |
| Astrocyte | Slow (100h) | Slow | No | Normal tissue |
| Microglia | No | Yes | No | Resident immune |
| Recruited Immune | No | Fast (30 μm/hr) | No | T-cells |
| Endothelial | Slow | No | No | Blood vessels |
| Necrotic | No | No | No | Dead cells |

### Environment Fields

The simulation tracks 3D concentration fields solved via finite differences:

| Field | Diffusion | Sources | Sinks |
|-------|-----------|---------|-------|
| O2 (mmHg) | 2000 μm²/s | Vessels | All cells |
| Glucose (mM) | 600 μm²/s | Vessels | All cells |
| VEGF (nM) | 100 μm²/s | Hypoxic tumor cells | Decay |

---

## Module Reference

### `causanta/simulate/` - Tumor Growth Simulation

| Module | Purpose | Key Functions/Classes |
|--------|---------|----------------------|
| `core.py` | Main simulation loop | `Simulation`, `run()` |
| `cells.py` | Cell agents and population | `CellPopulation`, `create_cell()` |
| `behaviors.py` | Cell behaviors | `advance_cell_cycle()`, `compute_migration()` |
| `ecdna.py` | **SIV mechanism** | `segregate_ecdna()`, `compute_egfr_expression()` |
| `environment.py` | Reaction-diffusion | `EnvironmentFields`, `DiffusionSolver` |
| `angiogenesis.py` | Vessel dynamics | `AngiogenesisManager` |
| `sweep.py` | Parameter sweeps | `ParameterSweep`, `run_sweep()` |
| `config.py` | Configuration | `SimulationConfig`, `load_config()` |

### `causanta/analyze/` - Causal Inference Tools

| Module | Purpose | Key Functions/Classes |
|--------|---------|----------------------|
| `iv.py` | **2SLS estimation** | `two_stage_least_squares()`, `IVEstimate`, `IVDiagnostics` |
| `effects.py` | Effect estimation | `estimate_effects_iv()`, `compare_iv_ols_performance()` |
| `bootstrap.py` | BCa confidence intervals | `bootstrap_iv_estimate()`, `bca_confidence_interval()` |
| `power.py` | Power analysis | `compute_power()`, `required_sample_size()`, `minimum_detectable_effect()` |
| `heterogeneity.py` | Stratified analysis | `HeterogeneityAnalysis`, `analyze_heterogeneity()` |
| `loader.py` | Data loading | `SimulationData`, `load_simulation()` |
| `discovery.py` | Causal structure learning | PC algorithm, GES |

### `causanta/visualize/` - Publication Figures

| Module | Purpose | Key Functions |
|--------|---------|---------------|
| `nature_style.py` | Journal styling | `apply_nature_style()`, `NATURE_COLORS` |
| `causal_figures.py` | Main figures | `create_main_figure_1()` through `create_main_figure_6()` |

### `causanta/graph/` - Causal DAG

| Module | Purpose | Key Classes |
|--------|---------|-------------|
| `causal_dag.py` | DAG representation | `CausalDAG`, `draw_dag()` |

---

## Data Flow

### From Simulation to Analysis

```
┌─────────────┐     ┌──────────────┐     ┌─────────────┐
│  SIMULATE   │────▶│   OUTPUT     │────▶│   ANALYZE   │
│             │     │              │     │             │
│ - Cells     │     │ cells_t*.tsv │     │ - Load data │
│ - ecDNA     │     │ lineage.tsv  │     │ - Run 2SLS  │
│ - Environ   │     │ environ*.tsv │     │ - Bootstrap │
│ - Division  │     │ params.json  │     │ - Power     │
└─────────────┘     └──────────────┘     └─────────────┘
```

### Key Data Files

| File | Contents | Used For |
|------|----------|----------|
| `cells_t{TIME}.tsv` | Cell states per timestep | Cross-sectional IV analysis |
| `lineage.tsv` | Division events with ecDNA segregation | Validating Binomial segregation |
| `environment_t{TIME}.tsv` | O2, glucose, VEGF fields | Understanding confounding |
| `params.json` | Ground truth parameters | Validation benchmarks |

### Variables in Cell Data

| Variable | Description | Role in IV |
|----------|-------------|------------|
| `ecDNA_count` | ecDNA copies | **Instrument (Z)** |
| `EGFR_expression` | EGFR protein level | **Treatment (D)** |
| `division_rate` | Proliferation rate | **Outcome (Y)** |
| `O2_local` | Local oxygen | Confounder (U) |
| `x`, `y` | Position | Spatial analysis |
| `cell_type` | Type code | Filtering |

---

## Ground Truth Parameters

These are the **true causal effects** embedded in the simulation. Your analysis should recover these.

### Effect Parameters

| Parameter | Symbol | Default | Equation | Interpretation |
|-----------|--------|---------|----------|----------------|
| `ecDNA_effect_on_division` | α | 0.30 | T_eff = T_base / (1 + α·log₂(1+ecDNA)) | 30% faster division per ecDNA doubling |
| `ecDNA_effect_on_VEGF` | β | 0.10 | S_eff = S_base · (1 + β·ecDNA) | 10% more VEGF per ecDNA copy |
| `ecDNA_effect_on_migration` | δ | 0.05 | v_eff = v_base · (1 + δ·ecDNA) | 5% faster migration per ecDNA copy |
| `ecDNA_effect_on_survival` | γ | 0.20 | a_eff = a_base / (1 + γ·ecDNA) | 20% survival boost per ecDNA copy |

### ecDNA Segregation Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `segregation_p` | 0.5 | Binomial probability (random) |
| `replication_fidelity` | 0.95 | ecDNA replication rate |
| `max_copies` | 100 | Biological cap |

### Success Criteria

For publication-ready validation:

| Criterion | Target | How to Measure |
|-----------|--------|----------------|
| IV accuracy | <15% bias | \|β̂_IV - β_true\| / β_true |
| CI coverage | ≥90% | Fraction of CIs containing true value |
| Power | ≥80% | At effect size 0.05, n=1000 |
| Segregation | p > 0.05 | KS test vs Binomial(N, 0.5) |
| First-stage F | >10 | Staiger & Stock rule |
| E-value | >2.0 | Robustness to confounding |

---

## Quick Start

### Installation

```bash
pip install -e .
```

### Run a Simulation

```bash
# Default: 1mm x 1mm domain, 168 hours
python -m causanta.simulate.core

# Custom duration
python -m causanta.simulate.core --hours 500

# Custom parameters
python -m causanta.simulate.core my_params.json
```

### Analyze Results

```bash
# Comprehensive analysis
python scripts/run_comprehensive_analysis.py output/run_*/data/

# Full validation study (baseline, effect sizes, power)
python scripts/run_simulation_study.py --study full --output-dir results/validation
```

### Generate Publication Figures

```python
from causanta.visualize import (
    apply_nature_style,
    create_main_figure_1,  # Conceptual framework (DAG)
    create_main_figure_2,  # Instrument validation (segregation)
    create_main_figure_3,  # Causal estimation (IV vs OLS)
    create_main_figure_4,  # Power analysis
    create_main_figure_5,  # Heterogeneity
    create_main_figure_6,  # Robustness (sensitivity)
)

apply_nature_style()
create_main_figure_1(output_dir="figures/")
```

---

## Analysis Capabilities

### IV Diagnostics

| Diagnostic | Purpose | Module |
|------------|---------|--------|
| F-statistic | Instrument strength | `iv.py` |
| Wu-Hausman test | Endogeneity detection | `iv.py` |
| Anderson-Rubin CI | Weak-IV robust inference | `iv.py` |
| Rosenbaum bounds | Sensitivity to hidden confounding | `iv.py` |
| E-value | Minimum confounding to explain away | `iv.py` |

### Statistical Methods

| Method | Purpose | Module |
|--------|---------|--------|
| 2SLS | Causal effect estimation | `iv.py` |
| BCa Bootstrap | Bias-corrected CIs | `bootstrap.py` |
| Power curves | Sample size planning | `power.py` |
| Stratified IV | Effect heterogeneity | `heterogeneity.py` |
| PC/GES algorithms | Causal discovery | `discovery.py` |

### Heterogeneity Strata

Effects can be estimated separately by:

| Stratum | Definition | Biological Rationale |
|---------|------------|---------------------|
| **Region** | Core (<100μm), Margin (100-300μm), Infiltrating (>300μm) | Different microenvironments |
| **Hypoxia** | Normoxic (>20 mmHg), Mild (10-20), Severe (<10) | O2-dependent pathways |
| **ecDNA burden** | Low (<10), Medium (10-30), High (>30) | Nonlinear dose-response |

---

## Documentation

All documentation is in the [`docs/`](docs/) folder.

### Manuscript

| Document | Description |
|----------|-------------|
| **[Main Manuscript](docs/manuscript/manuscript.md)** | Full paper: *Extrachromosomal DNA as a Causal Instrument for Spatial Multi-Omics* |
| **[Manuscript Guide](docs/manuscript/README.md)** | Figures, viewing instructions, citation info |
| **[Supplementary Methods](docs/supplementary_materials.md)** | Mathematical framework, 2SLS derivation, bootstrap, power analysis |

### Guides

| Document | Description |
|----------|-------------|
| **[Tutorial: Causal Inference](docs/tutorial.md)** | Step-by-step guide to causal effect estimation |
| **[Analysis Guide](docs/analysis_guide.md)** | Detailed statistical methods and interpretation |
| **[Immune System](docs/immune_system.md)** | Tumor-immune interaction modeling |

### Reference

| Document | Description |
|----------|-------------|
| **[Parameters](causanta/simulate/params/default.json)** | Full simulation configuration options |

---

## Repository Structure

```
causanta/
├── README.md                 # This file
├── pyproject.toml            # Package configuration
├── LICENSE
│
├── causanta/                 # Python package
│   ├── simulate/             # Tumor growth simulation
│   │   ├── core.py           # Main simulation loop (6 phases)
│   │   ├── cells.py          # Cell agents (9 types)
│   │   ├── behaviors.py      # Division, migration, death
│   │   ├── environment.py    # Reaction-diffusion solver
│   │   ├── ecdna.py          # ecDNA segregation (SIV mechanism)
│   │   ├── angiogenesis.py   # VEGF-driven vessel sprouting
│   │   ├── sweep.py          # Parameter sweep system
│   │   └── params/
│   │       └── default.json  # Default parameters
│   │
│   ├── analyze/              # Causal inference tools
│   │   ├── iv.py             # 2SLS estimation + diagnostics
│   │   ├── effects.py        # Effect estimation wrapper
│   │   ├── bootstrap.py      # BCa bootstrap CIs
│   │   ├── power.py          # Power analysis
│   │   ├── heterogeneity.py  # Stratified IV analysis
│   │   ├── loader.py         # Data loading utilities
│   │   └── discovery.py      # PC/GES causal discovery
│   │
│   ├── visualize/            # Publication figures
│   │   ├── nature_style.py   # Nature Methods styling
│   │   └── causal_figures.py # Figure 1-6 generation
│   │
│   └── graph/                # Causal DAG
│       └── causal_dag.py     # DAG construction
│
├── docs/                     # Documentation
│   ├── tutorial.md
│   ├── analysis_guide.md
│   ├── supplementary_materials.md
│   ├── immune_system.md
│   └── manuscript/           # Publication manuscript
│       ├── manuscript.md
│       └── figures/
│
├── scripts/                  # Analysis scripts
│   ├── run_comprehensive_analysis.py
│   ├── run_simulation_study.py
│   └── scenarios/            # Biological validation scenarios
│       ├── margin_invasion.py
│       ├── drug_response.py
│       ├── driver_vs_passenger.py
│       └── immune_selection.py
│
└── tests/                    # Test suite
    └── test_*.py
```

---

## Scope and Limitations

CAUSANTA is a **causal inference benchmark prototype**, not a validated GBM tissue simulator.

### What it is good for

- Testing whether causal inference methods (2SLS, IV, propensity matching) recover known effects
- Teaching causal inference concepts with biologically motivated examples
- Generating spatially structured synthetic data with embedded ground truth
- Benchmarking causal discovery algorithms against known DAG

### Known limitations

| Limitation | Details |
|------------|---------|
| **Phenomenological physics** | Environment PDEs mix units; parameters tuned for behavior, not accuracy |
| **Single tumor state** | Real GBM has OPC/NPC/AC/MES programs; we use one generic tumor type |
| **Single ecDNA axis** | Real ecDNA carries multiple oncogenes; we model one EGFR-like axis |
| **Simplified tissue physics** | 2D, endpoint collision detection, no force-based mechanics |
| **No calibration** | Parameters hand-set, not calibrated to experimental data |

### Future directions

- GBM state programs and niche-aware differentiation
- Multi-oncogene ecDNA modeling
- Calibration against spatial transcriptomics data
- 3D simulation and sectioning logic

---

## Development Notes

### For Claude Code Users

When working on this codebase:

1. **Ground truth is king** - All analysis validates against known effects in `default.json`
2. **The SIV mechanism is in `ecdna.py`** - This is the core innovation
3. **2SLS implementation is in `iv.py`** - Full diagnostics included
4. **Confounding comes from O2** - Hypoxia affects both EGFR and outcomes
5. **Segregation must be random** - Always validate with KS test vs Binomial
6. **Tests are in `tests/`** - Run with `pytest`

### Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| Binomial segregation | Matches biological mechanism; creates valid instrument |
| Log-scale effects | `log₂(1+ecDNA)` models diminishing returns |
| Hypoxia as confounder | Biologically realistic; creates identifiable bias |
| 2SLS over matching | More efficient; works with continuous treatment |
| BCa bootstrap | Handles skewness in IV estimates |

### Running Tests

```bash
pytest tests/ -v
```

### Output Management

Simulations can generate 100s of MB. To reduce:

```bash
# Keep only 2% of timesteps
python -m causanta.simulate.cleanup output/run_*/data --keep 0.02
```

---

## Session Finalization

### End-of-Session Prompt

**Copy and paste this prompt at the end of your Claude Code session to finalize changes:**

```
Please finalize my session with the following steps:

1. **Verify all code changes:**
   - Run `pytest tests/ -v` and fix any failures
   - Check that all modified modules have consistent imports
   - Verify no syntax errors in changed files

2. **Update documentation:**
   - If I modified any analysis methods, update docs/supplementary_materials.md
   - If I changed simulation behavior, update the relevant section in README.md
   - If I added new functions/classes, add them to the Module Reference tables in README.md

3. **Check manuscript consistency:**
   - If my changes affect the scientific claims, flag sections of docs/manuscript/manuscript.md that may need updates
   - Verify figure references still match if I changed visualization code

4. **Update the changelog:**
   - Add an entry to the "Change Log" section in README.md with today's date
   - Summarize what was changed, added, fixed, or removed
   - Include file paths for major changes

5. **Commit with descriptive message:**
   - Stage all related changes
   - Write a commit message summarizing the session's work

6. **Final status report:**
   - List any TODOs or follow-up items for next session
   - Note any known issues introduced or discovered
```

### What the Finalization Covers

| Check | Purpose |
|-------|---------|
| Test suite | Catch regressions |
| Module reference | Keep README accurate |
| Supplementary materials | Maintain methods documentation |
| Manuscript flags | Prevent paper/code divergence |
| Changelog | Track project history |
| Commit | Clean git history |

---

## Change Log

### 2024-04-23

**Major documentation overhaul and infrastructure additions**

- **Added:** Comprehensive causal analysis infrastructure
  - `causanta/analyze/bootstrap.py` - BCa bootstrap confidence intervals
  - `causanta/analyze/power.py` - Power analysis and sample size calculation
  - `causanta/analyze/heterogeneity.py` - Stratified IV analysis by region/hypoxia
  - `causanta/simulate/sweep.py` - Parameter sweep system
  - `causanta/visualize/` - Publication figure generation with Nature Methods styling

- **Added:** Biological validation scenarios in `scripts/scenarios/`
  - `margin_invasion.py` - Does ecDNA drive invasion?
  - `drug_response.py` - Does ecDNA predict drug response?
  - `driver_vs_passenger.py` - Can IV distinguish drivers from passengers?
  - `immune_selection.py` - Does immune pressure select for/against ecDNA?

- **Moved:** `manuscript/` → `docs/manuscript/`
  - Converted Word document to GitHub-renderable Markdown
  - Added CSS styling for academic paper appearance
  - Created publication figures (DAG, segregation, IV estimates, spatial)

- **Updated:** README.md expanded for collaborator onboarding
  - Added scientific background and SIV framework explanation
  - Added architecture overview with 6-phase simulation loop
  - Added module reference tables for all subsystems
  - Added ground truth parameters with equations
  - Added session finalization prompt and changelog

- **Fixed:** `causanta/analyze/iv.py` - `se_first` variable was undefined
- **Fixed:** `causanta/analyze/heterogeneity.py` - Missing `.get()` defaults for environment fields
- **Removed:** Deprecated `old/` directory contents
- **Removed:** Generated HTML docs (now using Markdown)

---

## Citation

If you use CAUSANTA in your research, please cite:

```bibtex
@article{craig2024ecdna,
  title = {Extrachromosomal DNA as a Causal Instrument for Spatial Multi-Omics},
  author = {Craig, David W. and Rodin, Andrei S.},
  year = {2024},
  institution = {City of Hope}
}
```

See the [full manuscript](docs/manuscript/manuscript.md) for details.

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

## Authors

- **David W. Craig** - Department of Integrative Translational Science, City of Hope
- **Andrei S. Rodin** - Department of Computational Medicine, City of Hope

*Correspondence: dacraig@coh.org*
