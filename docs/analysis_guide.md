# CAUSANTA Analysis Guide

## Statistical Methods for Causal Inference

This guide documents the analysis methods in CAUSANTA: causal effect estimation, statistical diagnostics, and interpretation guidelines.

---

## Table of Contents

1. [Overview of Causal Inference](#1-overview-of-causal-inference)
2. [The ecDNA Instrument](#2-the-ecdna-instrument)
3. [Effect Estimation Methods](#3-effect-estimation-methods)
4. [Statistical Diagnostics](#4-statistical-diagnostics)
5. [Working with Lineage Data](#5-working-with-lineage-data)
6. [Spatial Analysis](#6-spatial-analysis)
7. [Time-Series Analysis](#7-time-series-analysis)
8. [Common Pitfalls](#8-common-pitfalls)
9. [Code Examples](#9-code-examples)

---

## 1. Overview of Causal Inference

### The Fundamental Problem

In observational data, **correlation does not imply causation**. When we observe cells with high ecDNA dividing faster, several explanations exist:

| Explanation | Mechanism | Implication |
|-------------|-----------|-------------|
| **Causal** | ecDNA directly causes faster division | Effect is real |
| **Reverse causation** | Fast-dividing cells accumulate more ecDNA | Direction is backwards |
| **Confounding** | A third variable (e.g., O2) causes both | Association is spurious |

Standard regression cannot distinguish these scenarios when confounders are unobserved.

### Why CAUSANTA Enables Causal Inference

CAUSANTA exploits the fact that **ecDNA segregation acts as a natural randomizer**.

When a parent cell divides:

- Parent has 20 ecDNA copies, specific location, specific O2 level
- Daughter A receives `Binomial(20, 0.5)` = 12 copies
- Daughter B receives the remaining 8 copies

**Critically, the siblings share:**
- Same parent (genetics identical)
- Same location (environment identical at division)
- Same time (temporal confounders identical)

**The only difference is ecDNA count, which was assigned randomly.**

This mimics a randomized controlled trial within each cell division.

---

## 2. The ecDNA Instrument

### Instrument Validity Conditions

For ecDNA to be a valid instrumental variable (IV), three conditions must hold:

**1. Relevance**: The instrument must affect the exposure.

- ecDNA copy number directly affects oncogene expression (gene dosage effect)
- More copies = more mRNA = more protein

**2. Independence** (Randomization): The instrument must be independent of confounders.

- Segregation follows `Binomial(N, 0.5)` regardless of cell state
- The number of copies a daughter receives is random

**3. Exclusion Restriction**: The instrument affects the outcome *only* through the exposure.

- ecDNA affects phenotype solely via gene expression
- There is no direct ecDNA → behavior pathway

### The Causal Structure

```
ecDNA (Instrument) → Gene Expression (Exposure) → Phenotype (Outcome)
                                ↑                        ↑
                                └─────── Confounder ─────┘
                                        (e.g., O2)
```

The key insight: **confounders can affect both expression and phenotype, but they cannot affect ecDNA** because segregation is random.

### Segregation Statistics

The binomial segregation model predicts specific properties:

| Statistic | Formula | Example (N=20) |
|-----------|---------|----------------|
| Mean daughter fraction | 0.5 | 10 copies |
| Variance of fraction | p(1-p)/N | 0.0125 |
| Standard deviation | sqrt(0.25/N) | 0.112 |
| 95% CI for fraction | 0.5 +/- 1.96*SD | [0.28, 0.72] |

### Testing Segregation

```python
from scipy import stats
import numpy as np

# From lineage data: compute daughter fraction for each division
daughter_fractions = [
    d['daughter_ecDNA'] / d['parent_ecDNA_before']
    for d in lineage if d['parent_ecDNA_before'] > 0
]

# Test if mean equals expected 0.5
t_stat, p_value = stats.ttest_1samp(daughter_fractions, 0.5)

print(f"Mean fraction: {np.mean(daughter_fractions):.3f}")
print(f"P-value (H0: mean=0.5): {p_value:.4f}")

# If p < 0.05, segregation may be biased (violates independence)
```

---

## 3. Effect Estimation Methods

### Method 1: Ordinary Least Squares (OLS)

The simplest approach: regress outcome directly on ecDNA.

**Model:**
```
Y = beta_0 + beta_1 * ecDNA + epsilon
```

**Example:**
```
division_time = beta_0 + beta_1 * log2(ecDNA + 1) + epsilon
```

**Interpretation:** beta_1 = change in outcome per unit change in log2(ecDNA).

**Problem:** OLS is biased if confounders exist. If oxygen affects both ecDNA accumulation (cells in high O2 survive longer) and division rate (more resources), then beta_1 captures both the causal effect and the confounding.

### Method 2: Two-Stage Least Squares (2SLS)

Uses the instrument to isolate causal variation.

**Stage 1: Predict exposure from instrument**
```
gene_expr = gamma_0 + gamma_1 * ecDNA + eta
gene_expr_hat = gamma_0 + gamma_1 * ecDNA  # Predicted values
```

**Stage 2: Regress outcome on predicted exposure**
```
division_time = beta_0 + beta_1 * gene_expr_hat + epsilon
```

**Interpretation:** beta_1 is the **causal effect** of gene expression on division time. Confounding is removed because we only use variation in expression that comes from the random instrument.

**In CAUSANTA:** Since gene expression is modeled implicitly (ecDNA directly affects phenotype), we can use ecDNA as both instrument and exposure in a "reduced form" analysis.

### Method 3: Sibling Comparison

Compares cells that share a parent, eliminating all shared confounders.

For each division event:
- Parent: ecDNA = N
- Daughter A: ecDNA = N_A
- Daughter B: ecDNA = N_B (where N_A + N_B = N)

**The sibling difference:**
```
delta_ecDNA = N_A - N_B     # Random, mean = 0
delta_Outcome = Y_A - Y_B

# Regression on differences
delta_Outcome = beta * delta_ecDNA + epsilon
```

This design removes **all** confounders shared between siblings (genetics, location at birth, time of origin). The only remaining variation is the random segregation difference.

### Method 4: Fixed Effects Panel Regression

Controls for cell-level and time-level confounders using panel data.

**Model:**
```
Y_it = beta * ecDNA_it + alpha_i + tau_t + epsilon_it
```

Where:
- `Y_it` = outcome for cell i at time t
- `alpha_i` = cell fixed effect (controls lineage confounders)
- `tau_t` = time fixed effect (controls temporal confounders)

This approach tracks individual cells over time and uses within-cell variation in ecDNA (after division events) to estimate causal effects.

---

## 4. Statistical Diagnostics

### Checking Instrument Strength

A weak instrument leads to biased 2SLS estimates.

**F-Statistic Test:**

First-stage regression: `gene_expr ~ ecDNA + controls`

The F-statistic tests whether the coefficient on ecDNA is significantly different from zero.

| F-Statistic | Interpretation |
|-------------|----------------|
| F > 10 | Strong instrument, 2SLS is reliable |
| F in [5, 10] | Moderate, some bias possible |
| F < 5 | Weak instrument, do NOT trust 2SLS |

In CAUSANTA, the instrument is strong by design (ecDNA directly determines phenotype), but always verify.

### Checking for Confounding

Compare estimates across methods:

| Method | Estimate | If Confounded... |
|--------|----------|------------------|
| OLS | beta_OLS | Biased toward confounder effect |
| 2SLS/IV | beta_IV | Unbiased (if instrument valid) |
| Sibling | beta_sib | Unbiased (controls shared confounders) |

**Key diagnostic:** If beta_OLS differs significantly from beta_IV, confounding is present.

### Heterogeneity Analysis

Check if effects vary by subgroup:

```python
# Effect by location (near vs far from vessel)
near_vessel = cells[cells['distance_to_vessel'] < 50]
far_vessel = cells[cells['distance_to_vessel'] >= 50]

effect_near = estimate_effect(near_vessel)
effect_far = estimate_effect(far_vessel)

# If effects differ substantially, there may be:
# - Effect modification (true heterogeneity)
# - Residual confounding (invalid instrument in subgroup)
```

---

## 5. Working with Lineage Data

### Data Structure

The `lineage.tsv` file records every division event:

| Column | Description |
|--------|-------------|
| time_hr | Simulation time of division |
| parent_id | ID of parent cell |
| parent_ecDNA_before | Parent's ecDNA before segregation |
| parent_ecDNA_after | Parent's ecDNA after segregation |
| daughter_id | ID of new daughter cell |
| daughter_ecDNA | Daughter's ecDNA count |
| x, y | Location of division |
| generation | Generation number in lineage |

### Reconstructing Lineage Trees

```python
import networkx as nx

def build_lineage_tree(lineage_records):
    """Build directed graph of cell lineages."""
    G = nx.DiGraph()
    
    for record in lineage_records:
        parent = record['parent_id']
        daughter = record['daughter_id']
        
        # Add nodes with ecDNA info
        G.add_node(parent, ecDNA=record['parent_ecDNA_after'])
        G.add_node(daughter, ecDNA=record['daughter_ecDNA'])
        
        # Add edge (division event)
        G.add_edge(parent, daughter, 
                   time=record['time_hr'],
                   x=record['x'],
                   y=record['y'])
    
    return G

# Find all descendants of original tumor
tree = build_lineage_tree(lineage)
original_tumor_id = 3030  # From config
descendants = nx.descendants(tree, original_tumor_id)
```

### Tracking ecDNA Dynamics

```python
def compute_ecdna_trajectory(lineage, cell_id):
    """Track ecDNA changes through a cell's history."""
    trajectory = []
    current = cell_id
    
    while True:
        # Find this cell's birth record
        records = [r for r in lineage if r['daughter_id'] == current]
        if not records:
            break
        
        record = records[0]
        trajectory.append({
            'time': record['time_hr'],
            'ecDNA': record['daughter_ecDNA'],
            'event': 'birth'
        })
        
        current = record['parent_id']
    
    return list(reversed(trajectory))
```

---

## 6. Spatial Analysis

### Nutrient Gradients

Oxygen and glucose form gradients from vessels into the tissue:

| Distance from Vessel | O2 (mmHg) | State |
|---------------------|-----------|-------|
| 0-50 um | 40-60 | Normoxic |
| 50-100 um | 20-40 | Mildly hypoxic |
| 100-200 um | 10-20 | Hypoxic threshold |
| >200 um | <10 | Severely hypoxic |

The hypoxia threshold (~10 mmHg) is critical: below this, cells switch to the "Go" phenotype (increased migration) and cannot proliferate.

### Spatial Confounding

Location can confound ecDNA-phenotype relationships if selection operates spatially.

```python
def check_spatial_confounding(cells):
    """Test if ecDNA correlates with location."""
    from scipy.stats import spearmanr
    
    distances = compute_vessel_distances(cells)
    ecdna = [c['ecDNA_count'] for c in cells]
    
    corr, p_value = spearmanr(distances, ecdna)
    
    if p_value < 0.05:
        print(f"WARNING: ecDNA correlates with vessel distance")
        print(f"  Correlation: {corr:.3f}, p={p_value:.4f}")
        print(f"  This may indicate spatial confounding")
```

If high-ecDNA cells are more likely to survive in certain locations, this can bias estimates.

### Spatial Regression

Control for location in regressions:

```python
import statsmodels.api as sm

# Add spatial controls
cells_df['x_centered'] = cells_df['x'] - cells_df['x'].mean()
cells_df['y_centered'] = cells_df['y'] - cells_df['y'].mean()
cells_df['dist_to_center'] = np.sqrt(
    cells_df['x_centered']**2 + cells_df['y_centered']**2
)

# Model with spatial controls
model = sm.OLS(
    cells_df['division_time'],
    sm.add_constant(cells_df[[
        'log_ecdna', 'O2_local', 
        'dist_to_center', 'x_centered', 'y_centered'
    ]])
)
results = model.fit()
```

---

## 7. Time-Series Analysis

### Growth Curve Fitting

Fit standard growth models to tumor dynamics:

| Model | Equation | Use Case |
|-------|----------|----------|
| **Exponential** | N(t) = N0 * exp(r*t) | Unlimited growth, constant doubling time |
| **Logistic** | N(t) = K / (1 + ((K-N0)/N0) * exp(-r*t)) | Carrying capacity K, S-shaped curve |
| **Gompertz** | N(t) = K * exp(-exp(a - b*t)) | Asymmetric S-curve, common for tumors |

```python
from scipy.optimize import curve_fit
import numpy as np

def exponential(t, N0, r):
    return N0 * np.exp(r * t)

def logistic(t, K, r, N0):
    return K / (1 + ((K - N0) / N0) * np.exp(-r * t))

# Fit to data
time = summary_df['step'].values
tumor = summary_df['tumor_cells'].values

# Exponential fit
popt_exp, _ = curve_fit(exponential, time, tumor, p0=[1, 0.01])
doubling_time = np.log(2) / popt_exp[1]
print(f"Exponential: r={popt_exp[1]:.4f}, T_d={doubling_time:.1f}h")

# Logistic fit
popt_log, _ = curve_fit(logistic, time, tumor, p0=[1000, 0.01, 1])
print(f"Logistic: K={popt_log[0]:.0f}, r={popt_log[1]:.4f}")
```

### Phase Detection

Identify growth phases automatically:

```python
def detect_growth_phases(time, tumor_count):
    """Detect lag, exponential, and plateau phases."""
    log_count = np.log(np.maximum(tumor_count, 1))
    
    # Local growth rate (derivative of log)
    growth_rate = np.gradient(log_count, time)
    
    threshold = 0.01  # 1% per hour minimum for exponential
    
    phases = []
    for i, rate in enumerate(growth_rate):
        if rate < threshold:
            phases.append('lag' if i < len(growth_rate) // 3 else 'plateau')
        else:
            phases.append('exponential')
    
    return phases
```

---

## 8. Common Pitfalls

### Pitfall 1: Selection Bias

**Problem:** We only observe surviving cells.

If high-ecDNA cells die more frequently, the survivors are "special" (lucky), and their phenotypes may not represent typical high-ecDNA cells.

**Solutions:**
- Use lineage data to include all cells (dead and alive)
- Model death as a competing risk
- Use inverse probability weighting

### Pitfall 2: Time-Varying Confounding

**Problem:** The environment changes over time.

- At t=0: High O2 everywhere
- At t=500: Low O2 near tumor core

Cells born at different times face different environments, even with the same ecDNA count.

**Solutions:**
- Include time fixed effects
- Use within-time-window comparisons
- Control for local environment at time of measurement

### Pitfall 3: Measurement Error

**Problem:** In real data, ecDNA count is estimated with error.

Classical measurement error in X attenuates (biases toward zero) OLS estimates:
```
beta_OLS -> beta_true * reliability
```

**Good news:** IV methods are robust to measurement error in the exposure. This is another reason to use ecDNA as an instrument.

### Pitfall 4: Multiple Testing

When testing multiple causal effects (alpha, beta, delta, gamma), adjust for multiplicity:

```python
from statsmodels.stats.multitest import multipletests

p_values = [p_alpha, p_beta, p_delta, p_gamma]

# Bonferroni correction
alpha_adjusted = 0.05 / len(p_values)

# Or FDR correction (Benjamini-Hochberg)
reject, p_corrected, _, _ = multipletests(p_values, method='fdr_bh')
```

---

## 9. Code Examples

### Complete Analysis Pipeline

```python
"""Complete CAUSANTA analysis pipeline."""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

from causanta.analyze.loader import (
    load_cells_tsv,
    load_lineage_tsv, 
    load_summary_log
)
from causanta.analyze.effects import run_causal_analysis


def analyze_simulation(output_dir: Path):
    """Run complete analysis on simulation output."""
    
    # 1. Load data
    print("Loading data...")
    summary = load_summary_log(output_dir / "summary.log")
    lineage = load_lineage_tsv(output_dir / "lineage.tsv")
    
    cell_files = sorted(output_dir.glob("cells_t*.tsv"))
    final_cells = load_cells_tsv(cell_files[-1])
    
    # 2. Basic statistics
    print("\n=== Simulation Summary ===")
    print(f"Duration: {summary[-1]['step']} hours")
    print(f"Final tumor cells: {summary[-1]['tumor_cells']}")
    print(f"Total divisions: {len(lineage)}")
    
    # 3. Segregation analysis
    print("\n=== Segregation Analysis ===")
    daughter_fractions = []
    for div in lineage:
        if div['parent_ecDNA_before'] > 0:
            frac = div['daughter_ecDNA'] / div['parent_ecDNA_before']
            daughter_fractions.append(frac)
    
    if daughter_fractions:
        mean_frac = np.mean(daughter_fractions)
        t_stat, p_val = stats.ttest_1samp(daughter_fractions, 0.5)
        print(f"Mean daughter fraction: {mean_frac:.3f} (expected: 0.500)")
        print(f"T-test p-value (H0: mean=0.5): {p_val:.4f}")
    
    # 4. Growth analysis
    print("\n=== Growth Analysis ===")
    time = np.array([s['step'] for s in summary])
    tumor = np.array([s['tumor_cells'] for s in summary])
    
    if np.any(tumor > 0):
        log_tumor = np.log(tumor[tumor > 0])
        time_pos = time[tumor > 0]
        slope, intercept = np.polyfit(time_pos, log_tumor, 1)
        doubling_time = np.log(2) / slope if slope > 0 else float('inf')
        print(f"Growth rate: {slope:.4f} per hour")
        print(f"Doubling time: {doubling_time:.1f} hours")
    
    # 5. ecDNA distribution
    print("\n=== ecDNA Distribution (Final) ===")
    tumor_cells = [c for c in final_cells if c['cell_type'] == 6]
    if tumor_cells:
        ecdna_counts = [c['ecDNA_count'] for c in tumor_cells]
        print(f"N tumor cells: {len(tumor_cells)}")
        print(f"ecDNA mean: {np.mean(ecdna_counts):.1f}")
        print(f"ecDNA range: [{min(ecdna_counts)}, {max(ecdna_counts)}]")
    
    # 6. Causal effect estimation
    print("\n=== Causal Effect Estimation ===")
    analysis = run_causal_analysis(output_dir)
    
    return analysis


if __name__ == '__main__':
    import sys
    output_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("output/latest")
    results = analyze_simulation(output_dir)
```

### Visualization Functions

```python
"""Visualization functions for CAUSANTA analysis."""

import matplotlib.pyplot as plt
import numpy as np


def plot_growth_curve(summary, ax=None):
    """Plot tumor growth with phases annotated."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    time = [s['step'] for s in summary]
    tumor = [s['tumor_cells'] for s in summary]
    
    ax.semilogy(time, tumor, 'b-', linewidth=2, label='Tumor cells')
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Tumor cell count (log scale)')
    ax.set_title('Tumor Growth Curve')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    return ax


def plot_ecdna_evolution(lineage, ax=None):
    """Plot ecDNA distribution over time."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))
    
    times = [d['time_hr'] for d in lineage]
    ecdna = [d['daughter_ecDNA'] for d in lineage]
    
    if not times:
        return ax
    
    ax.scatter(times, ecdna, alpha=0.3, s=10)
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('ecDNA count at birth')
    ax.set_title('ecDNA Distribution Over Time')
    
    return ax


def plot_spatial_snapshot(cells, domain_size=1000, ax=None):
    """Plot spatial distribution of cells."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    
    colors = {
        0: '#4477AA',  # Neuron
        1: '#66CCEE',  # Astrocyte
        2: '#228833',  # Oligodendrocyte
        3: '#EE6677',  # Microglia
        4: '#AA3377',  # Endothelial
        5: '#BBBBBB',  # Pericyte
        6: '#CCBB44',  # Tumor
        7: '#EE8866',  # RecruitedImmune
        8: '#555555',  # Necrotic
    }
    
    for cell_type in sorted(colors.keys()):
        type_cells = [c for c in cells if c['cell_type'] == cell_type]
        if type_cells:
            x = [c['x'] for c in type_cells]
            y = [c['y'] for c in type_cells]
            size = 30 if cell_type == 6 else 10
            ax.scatter(x, y, c=colors[cell_type], s=size, alpha=0.7)
    
    ax.set_xlim(0, domain_size)
    ax.set_ylim(0, domain_size)
    ax.set_aspect('equal')
    ax.set_xlabel('X (um)')
    ax.set_ylabel('Y (um)')
    
    return ax
```

---

## Summary

CAUSANTA provides a rigorous framework for studying causal relationships in tumor biology:

| Feature | Benefit |
|---------|---------|
| Known ground truth | Validate causal inference methods |
| Natural instrument (ecDNA) | Enables IV estimation without RCT |
| Rich lineage data | Sibling comparisons, trajectories |
| Spatial structure | Study microenvironment effects |
| Temporal dynamics | Growth curve analysis |

The key to successful analysis:

1. **Verify instrument validity** (check segregation statistics)
2. **Use appropriate methods** (IV, sibling comparison)
3. **Check for confounding** (compare OLS vs IV estimates)
4. **Report uncertainty** (standard errors, confidence intervals)

---

## Further Reading

- **[Tutorial](tutorial.html)**: Step-by-step introduction to causal inference
- **[Immune System](immune_system.html)**: Tumor-immune interaction modeling
- **[Parameters Reference](../causanta/simulate/params/default.json)**: Full configuration options

### Key Papers

- Angrist & Pischke (2009). *Mostly Harmless Econometrics*
- Pearl (2009). *Causality: Models, Reasoning, and Inference*
- Hernan & Robins (2020). *Causal Inference: What If*

---

*CAUSANTA Analysis Guide v0.1.0*
