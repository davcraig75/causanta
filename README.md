# CAUSANTA

**Causal Analysis Using Somatic And Neighborhood Tissue Architecture**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-0.2.0-green.svg)](#)

---

## Executive Summary

CAUSANTA is a **mathematical demonstration**. It states a proposition, proves it, and verifies it numerically:

> **Proposition.** When extrachromosomal DNA (ecDNA) segregates randomly at cell division — `N_daughter ~ Binomial(N_parent, 0.5)` — ecDNA copy number is a *valid instrumental variable* for the causal effect of EGFR expression on cell phenotype.

The demonstration has two parts, and the repository exists to make both checkable:

1. **A structural proof.** The random-segregation mechanism satisfies the three defining conditions of an instrumental variable (relevance, independence, exclusion) by construction. Given those conditions, the classical IV identification theorem says the IV estimand equals the structural causal effect. This part is mathematics: it follows from the generative model, not from any dataset.

2. **A constructive numerical verification.** Mathematics tells us the estimand is correct *if* the assumptions hold; a simulation lets us confirm the estimator recovers the right number *when we know what the right number is*. The accompanying agent-based simulation generates synthetic tumor tissue with **known ground-truth causal effects** embedded by the model, alongside realistic confounding (hypoxia). We then check that the IV estimator recovers the planted effects while naive regression does not.

**Why both parts are needed.** In biology one rarely has simultaneously (1) a known true causal effect and (2) realistic confounding to defeat naive estimators. The simulation supplies both *by construction*, turning "ecDNA is a valid instrument" from an assertion into a statement that is *proved* (the segregation mechanism meets the IV conditions analytically) and then *verified* (2SLS recovers the embedded effects to within a stated tolerance; OLS is biased by the confounder). This is a proof-and-verification artifact, not a method offered for any particular practical use.

---

## Table of Contents

1. [Scientific Background](#scientific-background)
2. [The Demonstration: Somatic Instrumental Variables](#the-demonstration-somatic-instrumental-variables)
3. [Architecture Overview](#architecture-overview)
4. [Module Reference](#module-reference)
5. [Data Flow](#data-flow)
6. [Ground Truth Parameters](#ground-truth-parameters)
7. [Quick Start](#quick-start)
8. [Analysis Capabilities](#analysis-capabilities)
9. [Documentation](#documentation)
10. [Repository Structure](#repository-structure)
11. [Scope and Limitations](#scope-and-limitations)
12. [Change Log](#change-log)
13. [Citation](#citation)

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

## The Demonstration: Somatic Instrumental Variables

This section states the proposition, gives the structural proof, states the identification result it implies, and describes the numerical verification that closes the argument.

### Setup (the generative model)

Define, for a cell *i*:

- `Z` = ecDNA copy number — the candidate **instrument**
- `D` = EGFR expression — the **treatment**
- `Y` = phenotype (proliferation, migration, VEGF secretion) — the **outcome**
- `U` = local microenvironment, principally oxygen / hypoxia — the **confounder**

The model generates data from the following structural equations (an unconfounded edge `Z → D`, a confounded edge `D → Y`, and a back-door path through `U`):

```
    ecDNA (Instrument Z)
          ↓
    EGFR expression (Treatment D)  ←--  O2 (Confounder U)
          ↓                              ↓
    Proliferation (Outcome Y)  ←---------┘
```

```
    D = f(Z, U, ε_D)      with  ∂f/∂Z ≠ 0          (treatment depends on instrument and confounder)
    Y = g(D, U, ε_Y)                                (outcome depends on treatment and confounder, NOT on Z)
    Z' ~ Binomial(2Z, 1/2)  drawn independently of U at each division   (random segregation)
```

### Proposition (Somatic Instrumental Variable)

> Under the generative model above, `Z` (ecDNA copy number) is a **valid instrumental variable** for the effect of `D` (EGFR) on `Y` (phenotype). Consequently the IV estimand identifies the structural causal effect.

### Proof — the three IV conditions hold by construction

| Condition | Formal requirement | Why it holds in the model |
|-----------|--------------------|---------------------------|
| **Relevance** | `Cov(Z, D) ≠ 0` | Copy number enters EGFR through gene dosage, `X = (a + b·Z)·(…)` with `b > 0`, so `∂f/∂Z ≠ 0` and `Cov(Z, D) ≠ 0`. |
| **Independence (exogeneity)** | `Z ⊥ U` (conditional on lineage state) | The transmitted count is a binomial draw `Z' ~ Binomial(2Z, 1/2)` whose randomness uses **no** environmental information. The coin flips are independent of oxygen, so `Z` is independent of `U`. |
| **Exclusion** | `Z` affects `Y` only through `D` | In `Y = g(D, U, ε_Y)` the instrument `Z` does not appear: there is no direct edge `Z → Y`. Copy number influences phenotype *solely* via EGFR expression. |

The mechanism that makes independence true is the same one biology supplies: ecDNA lacks centromeres, so it cannot attach to the mitotic spindle and partitions by chance. Sister cells sharing one microenvironment receive different copy numbers — a coin flip *within* each confounding stratum, which is precisely randomization conditional on `U`. ∎

### Identification result (the theorem this invokes)

Given relevance, independence, and exclusion, the standard instrumental-variable identification theorem gives, for the linear case,

```
    β  =  Cov(Z, Y) / Cov(Z, D)
```

and `β` equals the **structural causal effect** of `D` on `Y`, with the confounding through `U` removed. The sample analog is **two-stage least squares**:

```
    Stage 1 (first stage):   D_i = π₀ + π₁·Z_i + η_i          (predict treatment from instrument)
    Stage 2 (second stage):  Y_i = β₀ + β₁·D̂_i + ε_i          (regress outcome on predicted treatment)
```

The coefficient `β₁` is the causal effect, purged of confounding. This is a mathematical consequence of the three conditions — it is not specific to ecDNA, biology, or this codebase.

### Numerical verification (closing the argument)

A proof guarantees the estimand is correct *when the assumptions hold*. The simulation lets us check the *estimator* on data where the true `β` is known by construction:

1. ecDNA replicates (approximately doubles) at division;
2. it segregates as `N_daughter ~ Binomial(N_parent, 0.5)`, independent of environment;
3. this produces within-lineage copy-number variation orthogonal to the microenvironment;
4. the planted structural effects (α, β, δ, γ — see [Ground Truth Parameters](#ground-truth-parameters)) are the numbers the estimator should return.

The verification succeeds when 2SLS recovers the planted effects to within the stated tolerance while OLS, which ignores the back-door path through `U`, is biased — confirming numerically what the proof establishes structurally.

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

All nine types are instantiated in every run: the parenchyma is populated with
neurons and glia at biologically realistic fractions, the vasculature is built as
an endothelial/pericyte network, the tumor is seeded with ecDNA, and recruited
immune + necrotic populations appear dynamically during the run.

| Type | Division | Migration (μm/hr) | ecDNA | Role |
|------|----------|-------------------|-------|------|
| Neuron | No | 0 | No | Normal parenchyma (excitable) |
| Astrocyte | Yes (168h) | 3 | No | Normal parenchyma (reactive near tumor) |
| Oligodendrocyte | No | 1 | No | Normal parenchyma (myelinating) |
| Microglia | No | 30 | No | Resident immune |
| Endothelial | Yes (60h) | 10 | No | Blood vessels |
| Pericyte | No | 5 | No | Vessel support |
| Tumor | Yes (24h) | 10 | Yes | Primary subject (ecDNA → EGFR) |
| Recruited Immune | Yes (36h) | 35 | No | Infiltrating CTL/NK/macrophage |
| Necrotic | No | 0 | No | Dead cells (hypoxic core) |

Tissue is filled to `cell_density_per_um2 = 0.002`; parenchymal/vascular metabolism
is scaled 0.5× so a populated tissue and a growing tumor coexist.

### Environment Fields

The simulation tracks 2D concentration fields solved via finite differences:

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
| `core.py` | Main simulation loop + CLI | `Simulation`, `run()`, `main()` |
| `cells.py` | Cell agents and population | `CellPopulation`, `create_cell()` |
| `behaviors.py` | Cell behaviors | `advance_cell_cycle()`, `compute_migration()` |
| `ecdna.py` | **SIV mechanism** | `replicate_and_segregate_ecdna()`, `compute_egfr_expression()` |
| `environment.py` | Reaction-diffusion | `EnvironmentFields`, `DiffusionSolver` |
| `angiogenesis.py` | Vessel dynamics | `AngiogenesisManager` |
| `initialization.py` | Initial tissue layout | `populate_normal_tissue()`, vascular network |
| `sweep.py` | Parameter sweeps | `SweepConfig`, `run_sweep()` |
| `config.py` | Configuration | `SimulationConfig`, `load_config()` |
| `io.py` | TSV read/write | Cell, environment, lineage I/O |
| `visualization.py` | Vega-Lite specs | `write_html_viewer()` |
| `reporting.py` | HTML run reports | `generate_report()` |
| `viewer.py` | Interactive SVG viewer | Standalone in-browser viewer |
| `movie.py` | Animated HTML | Time-lapse of cell positions |
| `cleanup.py` | Output pruning | `cleanup_data_keep_final_only()` |
| `benchmark.py` | Performance testing | Profiling utilities |
| `_numba_kernels.py` | JIT diffusion | Thomas-algorithm kernels |

### `causanta/analyze/` - Causal Inference Tools

| Module | Purpose | Key Functions/Classes |
|--------|---------|----------------------|
| `iv.py` | **2SLS estimation** | `two_stage_least_squares()`, `IVEstimate`, `compute_e_value()` |
| `effects.py` | Effect estimation | `run_causal_analysis()`, `analyze_ecDNA_segregation()` |
| `bootstrap.py` | BCa confidence intervals | `bootstrap_bca()`, `bootstrap_iv_estimate()`, `bootstrap_segregation_test()` |
| `power.py` | Power analysis | `iv_power_analytical()`, `sample_size_for_power()`, `minimum_detectable_effect()`, `power_curve()` |
| `heterogeneity.py` | Stratified analysis | `stratified_iv_by_region()`, `stratified_iv_by_hypoxia()`, `stratified_iv_by_ecDNA_level()` |
| `matching.py` | Propensity score methods | `propensity_score_matching()`, `inverse_propensity_weighting()` |
| `regression.py` | Regression adjustment | `linear_regression_adjustment()`, `doubly_robust_estimation()` |
| `sem.py` | Structural equation models | `mediation_analysis()`, `full_sem()`, `sobel_test()` |
| `sensitivity.py` | Robustness checks | `rosenbaum_bounds()`, `e_value_analysis()`, `placebo_test()` |
| `discovery.py` | Causal structure learning | `pc_algorithm()`, `ges_algorithm()`, `compare_to_ground_truth()` |
| `loader.py` | Data loading | `SimulationData`, `load_simulation_output()` |
| `report.py` | Combined HTML report | `generate_analysis_report()` |
| `cli.py` | `causanta-analyze` entry | CLI dispatch |

### `causanta/visualize/` - Publication Figures

| Module | Purpose | Key Functions |
|--------|---------|---------------|
| `nature_style.py` | Journal styling | `apply_nature_style()`, `NATURE_COLORS` |
| `causal_figures.py` | Main figures | `create_main_figure_1()` through `create_main_figure_6()` |

### `causanta/graph/` - Causal DAG

| Module | Purpose | Key Classes |
|--------|---------|-------------|
| `causal_dag.py` | DAG representation | `CausalDAG`, `CausalNode`, `CausalEdge`, `build_causanta_dag()` |

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

Per-cell snapshot columns (`cells_t{TIME}.tsv`):

| Variable | Description | Role in IV |
|----------|-------------|------------|
| `ecDNA_count` | ecDNA copies | **Instrument (Z)** |
| `egfr_expression` | EGFR protein level | **Treatment (D)** |
| `migration_rate` | Effective migration speed (μm/hr) | **Outcome (Y) — δ** |
| `VEGF_secretion` | Effective VEGF secretion rate | **Outcome (Y) — β** |
| `O2_local` | Local oxygen | Confounder (U) |
| `is_hypoxic` | Binary indicator (O2 < 18 mmHg) | Confounder (HIF-mediated) |
| `cell_cycle_phase` | G0/G1/S/G2/M | Division-rate proxy |
| `x`, `y` | Position | Spatial analysis |
| `cell_type` | Type code | Filtering |

The α (division) outcome is recovered from inter-division intervals in `lineage.tsv`,
not from a per-cell column — see `scripts/recover_structural_effects.py`.

---

## Ground Truth Parameters

These are the **true causal effects** embedded in the simulation. Your analysis should recover these.

### Effect Parameters

| Parameter | Symbol | Default | Equation | Interpretation |
|-----------|--------|---------|----------|----------------|
| `ecDNA_effect_on_division` | α | 0.30 | T_eff = T_base / (1 + α·log₂(1+EGFR)) | faster division with EGFR |
| `ecDNA_effect_on_VEGF` | β | 0.10 | S_eff = S_base · (1 + β·√EGFR) · (0.2 + 0.8·M) | more VEGF with EGFR (hypoxia-gated) |
| `ecDNA_effect_on_migration` | δ | 0.05 | v_eff = v_base · (1 + δ·EGFR) · (2× if hypoxic) | faster migration with EGFR |
| `ecDNA_effect_on_survival` | γ | 0.50 | a_eff = a_base / (1 + γ·log₂(1+EGFR)) | apoptosis resistance with EGFR |

The modulators act on EGFR expression X = (2.89 + 1.21·ecDNA)·(1 + κ_hyp·M)·noise, not directly on ecDNA; α/β/δ/γ are defined on the EGFR scale. M = is_hypoxic.

### ecDNA Segregation Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `segregation_p` | 0.5 | Binomial probability (random) |
| `replication_fidelity` | 0.95 | ecDNA replication rate |
| `max_copies` | 100 | Biological cap |

### Identifiability Status

| Parameter | Status | Reason |
|-----------|--------|--------|
| κ (first-stage slope) | Identified | Recovered from segregation regression |
| β (VEGF) | Identified | Inverted on κ_hyp = 0 runs |
| δ (migration) | Identified | Inverted on κ_hyp = 0 runs |
| α (division) | Identified (v21) | Nonlinear fit of inter-division intervals |
| γ (survival) | **Not identified** | Baseline apoptosis 5×10⁻⁵/hr → no detectable selection footprint |

### Verification Criteria

The numerical verification is considered to pass when the estimator meets these targets — i.e. when the simulation confirms what the proof predicts:

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

# With test dependencies
pip install -e ".[dev]"
```

Hard dependencies: `numpy>=1.24`, `scipy>=1.10`, `matplotlib>=3.7`, `numba>=0.58`.
Numba JITs the diffusion solver (~2× speedup) and is installed by default.

This registers two console entry points:

- `causanta-simulate` → `causanta.simulate.core:main`
- `causanta-analyze` → `causanta.analyze.cli:main`

### Run a Simulation

```bash
# Default: 1mm x 1mm domain, 168 hours
causanta-simulate

# Custom duration / seed
causanta-simulate --hours 500 --seed 43

# Custom parameters (any file under causanta/simulate/params/ works)
causanta-simulate causanta/simulate/params/large_baseline.json

# Equivalent module form
python -m causanta.simulate.core causanta/simulate/params/large_baseline.json
```

### Analyze Results

```bash
# Causanta-analyze entry point (HTML + JSON report)
causanta-analyze output/<run-dir>/data/

# Comprehensive analysis script
python scripts/run_comprehensive_analysis.py output/<run-dir>/data/

# Full validation study (baseline, effect sizes, power)
python scripts/run_simulation_study.py --study full --output-dir results/validation
```

### Generate Publication Figures

```python
from pathlib import Path
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
out = Path("figures")
out.mkdir(exist_ok=True)
create_main_figure_1(output_dir=out)
```

### Running Tests

```bash
pytest tests/ -v
```

---

## Analysis Capabilities

These are the estimators and diagnostics used to *verify* the proposition numerically — the apparatus that checks the estimator against the known ground truth. They are standard, published causal-inference methods, applied here only as instruments of verification.

### IV Diagnostics

| Diagnostic | Purpose | Module |
|------------|---------|--------|
| F-statistic | Instrument strength | `iv.py` |
| Wu-Hausman test | Endogeneity detection | `iv.py` |
| Anderson-Rubin CI | Weak-IV robust inference | `iv.py` |
| Rosenbaum bounds | Sensitivity to hidden confounding | `sensitivity.py` |
| E-value | Minimum confounding to explain away | `iv.py` / `sensitivity.py` |

### Statistical Methods

| Method | Purpose | Module |
|--------|---------|--------|
| 2SLS | Causal effect estimation | `iv.py` |
| BCa Bootstrap | Bias-corrected CIs | `bootstrap.py` |
| Power curves | Sample size planning | `power.py` |
| Stratified IV | Effect heterogeneity | `heterogeneity.py` |
| Propensity matching / IPW | Selection-on-observables | `matching.py` |
| Regression adjustment / DR | Covariate control | `regression.py` |
| Mediation / full SEM | Path decomposition | `sem.py` |
| PC / GES algorithms | Causal discovery | `discovery.py` |

### Heterogeneity Strata

Effects can be estimated separately by:

| Stratum | Definition | Biological Rationale |
|---------|------------|---------------------|
| **Region** | Core (<100μm), Margin (100-300μm), Infiltrating (>300μm) | Different microenvironments |
| **Hypoxia** | Normoxic (>20 mmHg), Mild (10-20), Severe (<10) | O2-dependent pathways |
| **ecDNA burden** | Low (<10), Medium (10-30), High (>30) | Nonlinear dose-response |

---

## Documentation

All documentation lives under [`docs/manuscript/`](docs/manuscript/).

### Theory & Specification

| Document | Description |
|----------|-------------|
| **[Theoretical Framework](docs/manuscript/theory.md)** | Complete specification: domain, cell model, environment PDEs, behavioral rules, SIV mechanism, causal inference methods |

### Manuscript

| Document | Description |
|----------|-------------|
| **[Main Manuscript](docs/manuscript/manuscript.md)** | Full paper: *Extrachromosomal DNA as a Causal Instrument for Spatial Multi-Omics* |
| **[Manuscript Guide](docs/manuscript/README.md)** | Figure index, viewing instructions, citation info |
| **[Supplementary Methods](docs/manuscript/supplementary_materials.md)** | Statistical methods for publication (2SLS, bootstrap, power) |
| **[Figure Legends](docs/manuscript/figures/figure_legends.md)** | Captions for all manuscript figures |

### Guides

| Document | Description |
|----------|-------------|
| **[Tutorial: Causal Inference](docs/manuscript/tutorial.md)** | Step-by-step guide to causal effect estimation |
| **[Analysis Guide](docs/manuscript/analysis_guide.md)** | Detailed statistical methods and interpretation |
| **[Immune System](docs/manuscript/immune_system.md)** | Tumor-immune interaction modeling |

### Interactive

| Document | Description |
|----------|-------------|
| **[Causal DAG (interactive)](docs/manuscript/causal_dag_interactive.html)** | Draggable structural causal model |
| **[Parameter Explorer](docs/manuscript/simulation_parameter_explorer.html)** | Drag α/β/δ/γ/κ/κ_hyp; live in-browser OLS-vs-2SLS Monte-Carlo |
| **[Multi-scale Explorer](docs/manuscript/simulation_multiscale_explorer.html)** | Section→Tissue→Cell reactive simulation diagram |

### Reference

| Document | Description |
|----------|-------------|
| **[Default Parameters](causanta/simulate/params/default.json)** | Full simulation configuration options |

---

## Repository Structure

```
causanta/
├── README.md                        # This file
├── pyproject.toml                   # Package metadata, deps, entry points
├── LICENSE                          # MIT
│
├── causanta/                        # Python package
│   ├── __init__.py
│   │
│   ├── simulate/                    # Tumor growth simulation
│   │   ├── core.py                  # Main 6-phase loop + CLI
│   │   ├── config.py                # Dataclass hierarchy from JSON
│   │   ├── domain.py                # 2D spatial domain
│   │   ├── cells.py                 # Cell agents (9 types)
│   │   ├── ecdna.py                 # SIV mechanism
│   │   ├── environment.py           # Reaction-diffusion PDE solver
│   │   ├── behaviors.py             # Division, migration, death, immune kill
│   │   ├── angiogenesis.py          # VEGF-driven vessel sprouting
│   │   ├── initialization.py        # Vascular network, tissue, tumor seeding
│   │   ├── io.py                    # TSV read/write
│   │   ├── visualization.py         # Vega-Lite spec generation
│   │   ├── reporting.py             # Per-run HTML report
│   │   ├── viewer.py                # Interactive SVG cell viewer
│   │   ├── movie.py                 # Animated HTML tumor growth
│   │   ├── shapes.py                # SVG path parsing
│   │   ├── cleanup.py               # Reduce intermediate TSV files
│   │   ├── benchmark.py             # Performance testing
│   │   ├── _numba_kernels.py        # JIT-compiled Thomas algorithm
│   │   └── params/                  # JSON configs (default + paper runs)
│   │
│   ├── analyze/                     # Causal inference tools
│   │   ├── cli.py                   # causanta-analyze entry point
│   │   ├── iv.py                    # 2SLS + diagnostics + E-value
│   │   ├── effects.py               # Effect estimation wrapper
│   │   ├── bootstrap.py             # BCa bootstrap CIs
│   │   ├── power.py                 # Power analysis
│   │   ├── heterogeneity.py         # Stratified IV analysis
│   │   ├── matching.py              # Propensity score methods
│   │   ├── regression.py            # Regression adjustment / DR
│   │   ├── sem.py                   # Structural equation models
│   │   ├── sensitivity.py           # Rosenbaum / E-value / placebo
│   │   ├── discovery.py             # PC / GES causal discovery
│   │   ├── loader.py                # Data loading utilities
│   │   └── report.py                # Combined HTML analysis report
│   │
│   ├── visualize/                   # Publication figures
│   │   ├── nature_style.py          # Nature Methods styling
│   │   └── causal_figures.py        # Figure 1–6 generation
│   │
│   └── graph/                       # Causal DAG
│       └── causal_dag.py            # DAG construction
│
├── docs/manuscript/                 # Documentation + manuscript
│   ├── manuscript.{md,docx}         # Main paper
│   ├── theory.{md,docx}             # Theoretical framework
│   ├── tutorial.{md,docx}           # Tutorial
│   ├── analysis_guide.{md,docx}     # Statistical methods guide
│   ├── supplementary_materials.{md,docx}
│   ├── immune_system.{md,docx}
│   ├── figures/                     # All manuscript figures (PNG + PDF)
│   ├── causal_dag_interactive.html
│   ├── simulation_parameter_explorer.html
│   └── simulation_multiscale_explorer.html
│
├── scripts/                         # Analysis + figure scripts
│   ├── run_comprehensive_analysis.py
│   ├── run_simulation_study.py
│   ├── run_all_sims_parallel.sh
│   ├── aggregate_multiseed_results.py
│   ├── multiseed_timeseries_analysis.py
│   ├── recover_structural_effects.py
│   ├── discovery_sweep.py
│   ├── discovery_multiseed_aggregate.py
│   ├── generate_figures.py
│   ├── generate_simulation_figure.py
│   ├── make_*.py                    # Per-figure regeneration scripts
│   └── scenarios/                   # Biological validation scenarios
│       ├── margin_invasion.py
│       ├── drug_response.py
│       ├── driver_vs_passenger.py
│       └── immune_selection.py
│
└── tests/                           # Scientific-contract test suite
    └── test_scientific_contracts.py
```

---

## Scope and Limitations

CAUSANTA is a **mathematical demonstration and causal-inference benchmark** — a proposition with a structural proof and a numerical verification. It is *not* a validated GBM tissue simulator, and it is not offered as a method or product for any clinical or experimental application; its purpose is to establish and check a statement about identifiability.

### What it is good for

- Demonstrating, with proof and numerical confirmation, that random ecDNA segregation makes copy number a valid instrument
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
| **γ not identified** | Baseline apoptosis is small enough that survival selection leaves no footprint in the configured regime |

### Future directions

- GBM state programs and niche-aware differentiation
- Multi-oncogene ecDNA modeling
- Calibration against spatial transcriptomics data
- 3D simulation and sectioning logic

---

## Change Log

Full history available via `git log`. Highlights from recent releases:

### v0.2.0 — 2026-06-01 (release)

- **Documentation consolidated** under `docs/manuscript/` (theory, tutorial, analysis guide, immune system, supplementary materials, all interactive HTML explorers, manuscript itself).
- **Repository cleanup** for release: removed legacy `.old/` artifacts, one-off scripts (presentation generator, diagram generator, MD→HTML renderer, param-celltype repair, ad-hoc batch runners), and `Instructions.md` rebuild seed.
- **Module/CLI alignment** with `pyproject.toml`: `causanta-simulate` and `causanta-analyze` are the canonical entry points.

### v0.2.0-pre — 2026-05-27 (manuscript v21)

- **α (division) identified** via nonlinear fit of inter-division intervals to `T_div = T_base / (1 + α·log₂(1+EGFR))` on the κ_hyp = 0 runs. Recovered α = 0.323 ± 0.010 (2 mm) and 0.307 ± 0.008 (6 mm) vs configured 0.30.
- **γ (survival) reframed as non-identified** — baseline apoptosis (5×10⁻⁵/hr) is too small to produce detectable selection. Deprecated `survival_enrichment` proxies in `analyze/effects.py`.
- **New figures**: `figure_recovery_identifiability.{png,pdf}` (Fig 9); interactive `simulation_parameter_explorer.html` (live OLS-vs-2SLS Monte-Carlo) and `simulation_multiscale_explorer.html` (Section→Tissue→Cell reactive diagram).
- Corrected division-event count for 6 mm baseline seed-43 lineage to 19,534 (was mislabeled 9,768 in earlier docs).

### v0.1.9 — 2026-05-22 (full tissue population fix)

- **Missing cell types restored.** Production param files defined `cell_types` for only `{Neuron, Endothelial, Pericyte, Tumor, Necrotic}`; Astrocyte/Oligodendrocyte/Microglia/RecruitedImmune were undefined, so `populate_normal_tissue` silently skipped them and immune recruitment never fired. Added the four definitions across all 31 production param files.
- **Tissue density fixed**: raised `cell_density_per_um2` from 0.001 to 0.002 (the prior value sat below the vascular network's own cell count, leaving the parenchyma empty); halved O₂/glucose consumption for parenchymal/vascular types so a populated tissue and a growing tumor coexist.
- All nine cell types now appear in every run; immune infiltration is scale-dependent (~390 cells at 2 mm; single digits at 6 mm — reproducing the immune-cold large-GBM phenotype).
- Causal benchmark preserved: OLS-vs-IV bias +8.7% (2 mm) / +12.5% (6 mm) baseline; first-stage F = 29k–150k; κ̂ ≈ 1.22; β = 0.100 and δ = 0.050 recovered exactly; PC F1 = 0.790–1.000; E-values 3.07 (migration) / 1.34 (VEGF).

### v0.1.8 — 2026-05-19 (E-value fix + canonical multi-seed reconciliation)

- **E-value corrected** to VanderWeele & Ding (2017): `E = RR + √(RR(RR−1))`, consolidated in a single canonical implementation (`causanta.analyze.iv.compute_e_value`).
- **Spatial heterogeneity** mode added (`ecDNA_migration_spatial`, default off); per-zone 2SLS recovers planted gradient (0.0302 / 0.0699 / 0.0500); the constant-δ run stays flat at 0.050 as a negative control.
- PC discovery switched from continuous `O2_local` to binary `is_hypoxic` (signal lives in the HIF-threshold step, not the linear oxygen field); F1 jumped to 1.000 at both baseline cells.

### Earlier history

See `git log` for the full record, including: per-hour snapshot retention default change (2026-05-18), multi-scale robustness analysis with the HIF effect knob (2026-04-30), and the initial causal analysis infrastructure (`bootstrap.py`, `power.py`, `heterogeneity.py`, `sweep.py`, `visualize/`) plus biological validation scenarios in `scripts/scenarios/` (2024-04-23).

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

- **David W. Craig** — Department of Integrative Translational Science, City of Hope
- **Andrei S. Rodin** — Department of Computational Medicine, City of Hope

*Correspondence: dacraig@coh.org*
