# Causal Inference with CAUSANTA

## A Practical Tutorial on Identifying Causal Effects in Spatial Biology

This tutorial teaches **causal inference** concepts using CAUSANTA as a training ground. You'll learn to:

1. Distinguish correlation from causation in biological data
2. Use instrumental variables (IV) to identify causal effects
3. Apply Two-Stage Least Squares (2SLS) regression
4. Validate causal methods against known ground truth

---

## The Problem: Correlation is Not Causation

### A Motivating Example

You observe tumor cells with more EGFR expression migrate faster. Is this causal?

**The naive interpretation**: *EGFR causes faster migration*.

**But wait** - what if both are caused by hypoxia?

- Hypoxia → Increased EGFR (HIF-1α stabilization)
- Hypoxia → Increased migration (Go-or-Grow)

This is **confounding**. The EGFR-migration correlation might be entirely spurious.

### The Confounding Problem Visualized

```
       ┌─────────────────────────────────────────────┐
       │           CONFOUNDED RELATIONSHIP           │
       └─────────────────────────────────────────────┘

                         Hypoxia
                        (CONFOUNDER)
                       ╱           ╲
                      ╱             ╲
                     ▼               ▼
               ┌──────────┐    ┌──────────┐
               │   EGFR   │ ── │Migration │
               │Expression│ ?  │  Speed   │
               └──────────┘    └──────────┘

       Is the EGFR→Migration arrow real, or is it
       entirely explained by the hypoxia→both paths?
```

**Standard regression can't tell you.** Even if you control for observed confounders, there might be unobserved ones.

---

## The Solution: Instrumental Variables

### What is an Instrument?

An **instrumental variable (IV)** is something that:

1. **Affects the exposure** (EGFR expression)
2. **Is as-if random** (not caused by confounders)
3. **Only affects the outcome through the exposure** (exclusion restriction)

If we find such a variable, we can use it to isolate the causal effect.

### ecDNA as a Natural Experiment

In CAUSANTA, **extrachromosomal DNA (ecDNA) copy number** is a perfect instrument:

| Property | Why ecDNA Satisfies It |
|----------|----------------------|
| **Relevance** | ecDNA carries EGFR oncogene → more copies = more EGFR expression |
| **Independence** | ecDNA segregates via `Binomial(N, 0.5)` at mitosis - essentially random |
| **Exclusion** | ecDNA affects cells *only through* gene expression |

The binomial segregation is key: when a cell with 20 ecDNA copies divides, each daughter gets `Binomial(20, 0.5)` copies - anywhere from 0 to 20. This randomness is like nature running a controlled experiment!

### The Causal DAG

```
       ┌─────────────────────────────────────────────┐
       │           CAUSANTA CAUSAL STRUCTURE          │
       └─────────────────────────────────────────────┘

         ┌────────────┐
         │ ecDNA_count│ ←── Segregation is RANDOM
         │(INSTRUMENT)│
         └─────┬──────┘
               │ Gene dosage effect
               ▼
         ┌────────────┐
         │   EGFR     │
         │ Expression │ ←── This is the EXPOSURE
         │ (EXPOSURE) │
         └─────┬──────┘
               │
       ┌───────┼───────┐
       ▼       ▼       ▼
   ┌──────┐ ┌──────┐ ┌──────┐
   │Prolif│ │Migra-│ │ VEGF │ ←── These are OUTCOMES
   │ rate │ │tion  │ │secre-│
   └──────┘ └──────┘ │tion  │
       ▲       ▲     └──────┘
       │       │
       └───┬───┘
           │
     ┌─────────┐
     │ Hypoxia │ ←── This is a CONFOUNDER
     │(O2_local)│    (affects outcomes, but NOT ecDNA)
     └─────────┘
```

Hypoxia confounds the EGFR→Outcome relationships, but it **cannot affect ecDNA counts** because segregation is random. This breaks the confounding!

---

## Two-Stage Least Squares (2SLS)

### The Method

2SLS uses the instrument to isolate causal variation:

**Stage 1: Predict exposure from instrument**
```
EGFR_expression = γ₀ + γ₁·ecDNA_count + ε₁
```

This captures only the variation in EGFR that's *caused by* ecDNA (the random part).

**Stage 2: Regress outcome on predicted exposure**
```
Migration_speed = β₀ + β₁·EGFR_predicted + ε₂
```

The coefficient β₁ is the **causal effect** of EGFR on migration.

### Why Does This Work?

The key insight: `EGFR_predicted` contains only variation from ecDNA, which is random. Any correlation between `EGFR_predicted` and migration must be causal, because the randomness ensures no confounding.

### Mathematical Intuition

If we denote:
- Z = ecDNA (instrument)
- X = EGFR expression (exposure)
- Y = Migration (outcome)
- U = Hypoxia (confounder)

The problem with OLS: `Cov(X, U) ≠ 0` (EGFR is correlated with hypoxia)

2SLS solution: Use `X̂ = E[X|Z]` instead. Since `Cov(Z, U) = 0` (ecDNA is random), we have `Cov(X̂, U) = 0`.

The 2SLS estimator:
```
β̂_IV = Cov(Y, Z) / Cov(X, Z)
```

This is the ratio of "how much does the outcome change with the instrument" to "how much does the exposure change with the instrument."

---

## Hands-On: Running the Analysis

### Step 1: Run a Simulation

```bash
python -m causanta.simulate.core --hours 300 --seed 42
```

This generates a tumor with known causal parameters (set in `default.json`):
- `ecDNA_effect_on_division` = 0.30 (α)
- `ecDNA_effect_on_migration` = 0.05 (δ)
- `ecDNA_effect_on_VEGF` = 0.10 (β)

### Step 2: Load the Data

```python
import pandas as pd
from pathlib import Path

# Load final cell states
data_dir = Path("output/run_20260413_093442/data")
cells = pd.read_csv(data_dir / "cells_t000300.tsv", sep="\t")

# Filter to tumor cells only
tumor = cells[cells["cell_type"] == 6].copy()
print(f"Tumor cells: {len(tumor)}")
```

### Step 3: Naive OLS (Wrong!)

```python
import numpy as np
from scipy import stats

# Naive regression: EGFR → Migration
X = tumor["egfr_expression"].values
Y = tumor["migration_rate"].values

slope, intercept, r, p, se = stats.linregress(X, Y)
print(f"OLS estimate: {slope:.4f} (biased by confounding!)")
```

This gives a biased estimate because hypoxia confounds the relationship.

### Step 4: 2SLS (Correct!)

```python
# Stage 1: ecDNA → EGFR
Z = tumor["ecDNA_count"].values
X = tumor["egfr_expression"].values

gamma1, gamma0, _, _, _ = stats.linregress(Z, X)
X_hat = gamma0 + gamma1 * Z  # Predicted EGFR

# Stage 2: EGFR_predicted → Migration
Y = tumor["migration_rate"].values
beta1, beta0, _, _, _ = stats.linregress(X_hat, Y)

print(f"2SLS estimate: {beta1:.4f}")
print(f"True effect: 0.05")  # From simulation parameters
```

### Step 5: Compare to Ground Truth

| Method | Estimate | True Value | Error |
|--------|----------|------------|-------|
| OLS | ~0.12 | 0.05 | 140% |
| 2SLS | ~0.05 | 0.05 | <5% |

The 2SLS estimate should be much closer to the true value!

---

## Checking Instrument Validity

### 1. Relevance: F-Statistic

The first-stage F-statistic should be > 10 (rule of thumb):

```python
# First stage regression
from scipy import stats

slope, intercept, r, p, se = stats.linregress(Z, X)
n = len(Z)
F = (r**2 / 1) / ((1 - r**2) / (n - 2))

print(f"First-stage F = {F:.1f}")
# Should be >> 10 for a strong instrument
```

### 2. Independence: Visual Check

Plot ecDNA against potential confounders:

```python
# ecDNA should NOT correlate with O2
r_o2, p_o2 = stats.pearsonr(
    tumor["ecDNA_count"], 
    tumor["O2_local"]
)
print(f"ecDNA-O2 correlation: r={r_o2:.3f}, p={p_o2:.3f}")
# Should be near zero and non-significant
```

### 3. Exclusion: Domain Knowledge

This must be argued, not tested. In CAUSANTA, we *know* ecDNA only affects phenotypes through gene expression because we coded it that way. In real data, this requires biological reasoning.

---

## The Go-or-Grow Example

### Biological Context

The **Go-or-Grow hypothesis** states that tumor cells trade off between proliferation and migration:
- Well-oxygenated cells: Proliferate (Grow)
- Hypoxic cells: Migrate (Go)

This creates confounding because hypoxia affects both behaviors.

### CAUSANTA Implementation

In the simulation:
- Hypoxic cells get a 2x migration boost (`hypoxia_invasion_boost`)
- Hypoxic cells can't proliferate (below O2 threshold)
- EGFR expression independently modulates both

### Disentangling Effects

Using 2SLS, we can separate:
1. **Direct EGFR effect on migration**: The causal parameter δ
2. **Hypoxia effect**: The Go-or-Grow response

```python
# Multiple regression with IV
# Y = migration
# X = EGFR (instrumented by ecDNA)
# W = is_hypoxic (control variable)

# This requires 2SLS with controls - see analysis_guide.html
```

---

## Summary: The Causal Inference Recipe

1. **Identify the question**: What causal effect do you want to estimate?

2. **Draw the DAG**: Map out all variables and their relationships

3. **Find an instrument**: Something that:
   - Affects the exposure
   - Is as-if random
   - Only works through the exposure

4. **Run 2SLS**:
   - Stage 1: Regress exposure on instrument
   - Stage 2: Regress outcome on predicted exposure

5. **Check validity**:
   - F-stat > 10 (strong instrument)
   - No correlation with confounders
   - Exclusion restriction holds

6. **Validate against ground truth** (in simulation) or **sensitivity analysis** (in real data)

---

## Further Reading

- **[Analysis Guide](analysis_guide.html)**: Detailed statistical methods
- **[Immune System](immune_system.html)**: Tumor-immune dynamics
- **Parameters Reference**: `causanta/simulate/params/default.json`

### Key Papers

- Angrist & Pischke (2009). *Mostly Harmless Econometrics*
- Pearl (2009). *Causality: Models, Reasoning, and Inference*
- Hernán & Robins (2020). *Causal Inference: What If*

---

## Exercises

1. **Vary the effect size**: Change `ecDNA_effect_on_migration` in `default.json` and verify 2SLS recovers the new value.

2. **Weak instrument**: Set `ecDNA_segregation_p` to 0.99 (almost no variation). What happens to the 2SLS estimate?

3. **Multiple outcomes**: Estimate causal effects on proliferation rate and VEGF secretion using the same instrument.

4. **Sensitivity analysis**: How much would unmeasured confounding need to bias the instrument for your conclusions to change?
