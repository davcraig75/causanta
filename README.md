# CAUSANTA

**Causal Analysis Using Somatic And Neighborhood Tissue Architecture**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

CAUSANTA is a simulation framework for benchmarking causal inference methods in spatial biology. It generates synthetic tumor tissue with **known ground-truth causal structure**, enabling rigorous validation of causal discovery algorithms.

The key innovation: **extrachromosomal DNA (ecDNA) segregation acts as a natural randomization mechanism**, satisfying instrumental variable assumptions and allowing causal effect identification from observational spatial data.

## Why CAUSANTA?

Causal inference in biology is hard because:
- Randomized experiments are often infeasible
- Observational data has confounding
- Ground truth is rarely known

CAUSANTA solves this by:
1. **Embedding a known causal DAG** into a realistic tissue simulation
2. **Using ecDNA as a Somatic Instrumental Variable (SIV)** - it segregates randomly during mitosis, creating natural variation that mimics randomization
3. **Generating rich spatial data** where you can validate whether your causal inference method recovers the true effects

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
# Command-line analysis
python -m causanta.analyze.cli output/run_*/data/

# Or use Python
python examples/analyze_simulation.py output/run_*/data/
```

---

## Documentation

| Document | Description |
|----------|-------------|
| **[Tutorial: Causal Inference](docs/tutorial.html)** | Step-by-step guide to causal effect estimation with CAUSANTA |
| **[Analysis Guide](docs/analysis_guide.html)** | Detailed statistical methods and interpretation |
| **[Immune System](docs/immune_system.html)** | Tumor-immune interaction modeling |
| **[Parameters Reference](causanta/simulate/params/default.json)** | Full configuration options |

---

## The Causal Structure

CAUSANTA embeds this causal DAG:

```
                    ┌─────────────────────────────────────┐
                    │         CAUSAL DAG                  │
                    └─────────────────────────────────────┘

    ┌──────────────┐         ┌──────────────────┐
    │  ecDNA_count │────────▶│ EGFR_expression  │
    │  (INSTRUMENT)│         │   (EXPOSURE)     │
    └──────────────┘         └────────┬─────────┘
           │                          │
           │                          ▼
           │                 ┌──────────────────┐
           │                 │  Proliferation   │◀──┐
           │                 │    Rate          │   │
           │                 └──────────────────┘   │
           │                          │             │
           │                          ▼             │
           │                 ┌──────────────────┐   │
           └────────────────▶│    Migration     │   │
                             │     Speed        │   │ O2_local
                             └──────────────────┘   │ (CONFOUNDER)
                                      │             │
                                      ▼             │
                             ┌──────────────────┐   │
                             │ VEGF_secretion   │───┘
                             │ → Angiogenesis   │
                             └──────────────────┘
```

### Why ecDNA is a Valid Instrument

1. **Relevance**: ecDNA copy number directly affects EGFR expression (gene dosage)
2. **Independence**: ecDNA segregates via `Binomial(N, 0.5)` during mitosis - essentially random
3. **Exclusion**: ecDNA affects phenotypes *only through* gene expression, not directly

This allows Two-Stage Least Squares (2SLS) regression to recover causal effects from observational data.

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
│   │   ├── core.py           # Main simulation loop
│   │   ├── cells.py          # Cell agents
│   │   ├── behaviors.py      # Cell behaviors (division, migration, death)
│   │   ├── environment.py    # Reaction-diffusion fields
│   │   ├── ecdna.py          # ecDNA segregation and effects
│   │   └── params/
│   │       └── default.json  # Default parameters
│   │
│   ├── analyze/              # Causal inference tools
│   │   ├── effects.py        # Causal effect estimation
│   │   ├── iv.py             # Instrumental variable regression
│   │   ├── regression.py     # OLS and 2SLS
│   │   └── discovery.py      # Causal structure learning
│   │
│   └── graph/                # Causal DAG representation
│       └── causal_dag.py     # DAG construction and visualization
│
├── docs/                     # Documentation
│   ├── tutorial.md           # Causal inference tutorial
│   ├── analysis_guide.md     # Statistical methods guide
│   ├── immune_system.md      # Immune modeling docs
│   └── assets/
│       └── style.css         # Documentation styling
│
├── examples/                 # Example scripts
│   └── analyze_simulation.py
│
└── scripts/                  # Utility scripts
    └── run_enhanced_analysis.py
```

---

## Configuration

All parameters are in [`causanta/simulate/params/default.json`](causanta/simulate/params/default.json).

### Key Sections

| Section | What It Controls |
|---------|------------------|
| `domain` | Simulation size (up to 12mm x 12mm) |
| `time` | Duration, output frequency |
| `cell_types` | 9 cell types with individual behaviors |
| `tumor_seeds` | Initial tumor location and ecDNA |
| `immune_recruitment` | Tumor-immune dynamics |
| `angiogenesis` | VEGF-driven vessel formation |

### Causal Effect Parameters

These are the **ground truth** effects you're trying to recover:

| Parameter | Effect | Equation |
|-----------|--------|----------|
| `ecDNA_effect_on_division` (α) | Faster division | T_eff = T_base / (1 + α·log₂(1+ecDNA)) |
| `ecDNA_effect_on_VEGF` (β) | More VEGF | S_eff = S_base · (1 + β·ecDNA) |
| `ecDNA_effect_on_migration` (δ) | Faster migration | v_eff = v_base · (1 + δ·ecDNA) |
| `ecDNA_effect_on_survival` (γ) | Apoptosis resistance | a_eff = a_base / (1 + γ·ecDNA) |

---

## Output Structure

Each simulation creates an organized output directory:

```
output/run_YYYYMMDD_HHMMSS/
├── data/
│   ├── cells_t000000.tsv      # Cell states per timestep
│   ├── environment_t000000.tsv # Environment fields
│   ├── lineage.tsv            # Division events with ecDNA segregation
│   ├── summary.log            # Per-step statistics
│   └── params.json            # Parameters used
│
├── figures/
│   └── simulation.vl.json     # Vega-Lite visualization spec
│
├── reports/
│   ├── report.html            # Comprehensive analysis report
│   └── viewer.html            # Interactive cell viewer
│
└── animations/
    └── tumor_growth.html      # Animated tumor growth visualization
```

### Cleanup

To reduce disk usage (simulations can generate 100s of MB):

```bash
# Keep only 2% of timesteps (first, last, and samples)
python -m causanta.simulate.cleanup output/run_*/data --keep 0.02
```

---

## Performance

| Domain | Cells | Speed | 168h Runtime |
|--------|-------|-------|--------------|
| 1mm x 1mm | ~5,000 | ~3 hr/s | ~1 min |
| 2mm x 2mm | ~20,000 | ~1 hr/s | ~3 min |

---

## Scope and Limitations

CAUSANTA is a **causal inference benchmark prototype**, not a validated GBM tissue simulator.

**What it is good for:**
- Testing whether causal inference methods (2SLS, IV, propensity matching) recover known effects
- Teaching causal inference concepts with biologically motivated examples
- Generating spatially structured synthetic data with embedded ground truth

**Known limitations:**
- **Phenomenological physics**: Environment reaction terms mix units (mmHg, mM, amol/hr) without conversion. Parameter values are tuned for qualitative behavior, not physical accuracy.
- **Single tumor state**: Real GBM has OPC-like, NPC-like, AC-like, and MES-like programs. CAUSANTA uses one generic tumor cell with ecDNA-driven modifiers.
- **Single ecDNA axis**: Real GBM ecDNA can carry multiple oncogenes with distinct dynamics. CAUSANTA models one EGFR-like axis.
- **Simplified tissue physics**: Endpoint collision detection rather than force-based mechanics. 2D section simulation without 3D-to-2D sectioning logic.
- **No calibration pipeline**: Parameters are hand-set, not calibrated against experimental spatial data. No uncertainty quantification.

**Future directions** (see `instructions.md` for detailed roadmap):
- GBM state programs and niche-aware logic
- Multi-oncogene ecDNA modeling
- Calibration and ensemble benchmarking
- Validation against measured GBM spatial patterns

---

## Citation

If you use CAUSANTA in your research, please cite:

```bibtex
@software{causanta2024,
  title = {CAUSANTA: Causal Analysis Using Somatic And Neighborhood Tissue Architecture},
  year = {2024},
  url = {https://github.com/yourusername/causanta}
}
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
