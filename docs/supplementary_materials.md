# CAUSANTA Supplementary Materials

## Statistical Methods for Nature Methods Publication

This document provides detailed methodological specifications for the CAUSANTA framework, suitable for peer review and reproducibility.

---

## Table of Contents

1. [Mathematical Framework](#1-mathematical-framework)
2. [Instrumental Variable Estimation](#2-instrumental-variable-estimation)
3. [Bootstrap Confidence Intervals](#3-bootstrap-confidence-intervals)
4. [Power Analysis](#4-power-analysis)
5. [Heterogeneity Analysis](#5-heterogeneity-analysis)
6. [Sensitivity Analysis](#6-sensitivity-analysis)
7. [Causal Discovery](#7-causal-discovery)
8. [Simulation Parameters](#8-simulation-parameters)
9. [Success Criteria](#9-success-criteria)
10. [Software Implementation](#10-software-implementation)

---

## 1. Mathematical Framework

### 1.1 Structural Causal Model

The CAUSANTA simulator implements the following structural equations:

**ecDNA Segregation (Instrument Generation)**
```
N_daughter ~ Binomial(N_parent, 0.5)
```

**Gene Expression (First Stage)**
```
EGFR_i = f(ecDNA_i) = baseline + slope * log2(1 + ecDNA_i) + noise
```

**Phenotypic Effects (Structural Equations)**
```
Division_rate_i = rate_base * (1 + α * log2(1 + ecDNA_i)) * hypoxia_effect(O2_i)
VEGF_i = VEGF_base * (1 + β * sqrt(EGFR_i)) * hypoxia_effect(O2_i)
Migration_i = v_base * (1 + δ * EGFR_i) * hypoxia_effect(O2_i)
Survival_i = a_base / (1 + γ * ecDNA_i) * hypoxia_effect(O2_i)
```

**Confounding Structure**
```
O2_i (unobserved common cause) → EGFR_i
O2_i (unobserved common cause) → Y_i (outcomes)
```

### 1.2 Ground Truth Parameters

| Symbol | Parameter | Default Value | Biological Interpretation |
|--------|-----------|---------------|---------------------------|
| α | ecDNA_effect_on_division | 0.30 | 30% faster division per doubling of ecDNA |
| β | ecDNA_effect_on_VEGF | 0.10 | 10% more VEGF per sqrt(EGFR) unit |
| δ | ecDNA_effect_on_migration | 0.05 | 5% faster migration per EGFR unit |
| γ | ecDNA_effect_on_survival | 0.50 | 50% survival boost per ecDNA copy |

### 1.3 Causal DAG

```
                    O2 (Confounder)
                    ↓         ↓
ecDNA (Z) → EGFR (D) → Outcomes (Y)
   ↑
   |
Binomial(N, 0.5)
```

**Key Properties:**
- Z ⊥ U : ecDNA segregation is independent of confounders (random mechanism)
- Z → D : ecDNA causally affects EGFR expression (gene dosage)
- Z ⫫ Y | D : ecDNA affects outcomes only through EGFR (exclusion restriction)

---

## 2. Instrumental Variable Estimation

### 2.1 Two-Stage Least Squares (2SLS)

**Stage 1 (First Stage):**
```
D_i = π_0 + π_1 * Z_i + π_2 * X_i + η_i
D̂_i = π̂_0 + π̂_1 * Z_i + π̂_2 * X_i
```

**Stage 2 (Second Stage):**
```
Y_i = β_0 + β_1 * D̂_i + β_2 * X_i + ε_i
```

Where:
- Z_i = ecDNA copy number (instrument)
- D_i = EGFR expression (endogenous treatment)
- Y_i = Phenotypic outcome
- X_i = Observed covariates (optional)

### 2.2 Standard Error Correction

The naive second-stage standard errors are incorrect. We use the corrected formula:

```
SE(β̂_IV) = sqrt(σ̂² * (X̂'X̂)^(-1))
```

Where σ̂² is estimated from residuals using original D (not D̂):
```
σ̂² = Σ(Y_i - β̂_0 - β̂_1 * D_i - β̂_2 * X_i)² / (n - k)
```

### 2.3 Instrument Strength Diagnostics

**F-Statistic:**
```
F = (R²_first / k) / ((1 - R²_first) / (n - k - 1))
```

| F-Statistic | Interpretation |
|-------------|----------------|
| F > 10 | Strong instrument (Staiger & Stock rule) |
| F > 16.38 | <5% maximal IV size distortion |
| F > 8.96 | <10% maximal IV size distortion |

**Partial R²:**
```
R²_partial = (R²_full - R²_restricted) / (1 - R²_restricted)
```

### 2.4 Wu-Hausman Endogeneity Test

Tests whether OLS and IV estimates differ significantly.

**Augmented Regression:**
```
Y_i = β_0 + β_1 * D_i + β_2 * X_i + θ * v̂_i + ε_i
```

Where v̂_i = D_i - D̂_i (first-stage residuals).

**Test:** H0: θ = 0 (no endogeneity)
- If p < 0.05, OLS is biased → use IV

### 2.5 Anderson-Rubin Confidence Intervals

Robust to weak instruments:

```
AR(β_0) = n * (SSR_restricted(β_0) - SSR_unrestricted) / SSR_unrestricted
```

Under H0: β = β_0, AR ~ χ²(1)

The AR confidence set inverts this test:
```
CI_AR = {β_0 : AR(β_0) ≤ χ²_{1,α}}
```

---

## 3. Bootstrap Confidence Intervals

### 3.1 BCa (Bias-Corrected and Accelerated) Bootstrap

The BCa method adjusts for bias and skewness in the bootstrap distribution.

**Algorithm:**
1. Compute θ̂ from original sample
2. Draw B bootstrap samples, compute θ̂*_b for each
3. Compute bias correction:
   ```
   ẑ_0 = Φ^(-1)(#{θ̂*_b < θ̂} / B)
   ```
4. Compute acceleration (jackknife influence):
   ```
   â = Σ(θ̂_{(-i)} - θ̂_{(·)})³ / (6 * [Σ(θ̂_{(-i)} - θ̂_{(·)})²]^(3/2))
   ```
5. Adjusted percentiles:
   ```
   α_1 = Φ(ẑ_0 + (ẑ_0 + z_{α/2}) / (1 - â(ẑ_0 + z_{α/2})))
   α_2 = Φ(ẑ_0 + (ẑ_0 + z_{1-α/2}) / (1 - â(ẑ_0 + z_{1-α/2})))
   ```
6. CI = [θ̂*_(α_1), θ̂*_(α_2)]

### 3.2 Bootstrap for IV Estimates

Each bootstrap replicate:
1. Resample (Z_i, D_i, Y_i) with replacement
2. Run full 2SLS procedure
3. Store IV coefficient

### 3.3 Testing Segregation Randomness

**Chi-Square Test:**
```
χ² = Σ (O_k - E_k)² / E_k
```

Where:
- O_k = Observed count in bin k
- E_k = n * P(daughter_fraction ∈ bin_k | Binomial)

**Kolmogorov-Smirnov Test:**
```
D = max|F_n(x) - F_0(x)|
```

Where F_0 is the theoretical CDF under Binomial(N, 0.5).

---

## 4. Power Analysis

### 4.1 Analytical Power for 2SLS

**Asymptotic Distribution:**
```
√n(β̂_IV - β) →d N(0, σ² / (π_1² * Var(Z)))
```

**Standard Error of IV:**
```
SE_IV ≈ sqrt(σ² / (n * R²_first * Var(Z)))
```

**Power:**
```
Power = 1 - Φ(z_{1-α/2} - |β| / SE_IV) + Φ(-z_{1-α/2} - |β| / SE_IV)
```

### 4.2 Required Sample Size

For target power (1 - β) at significance level α:

```
n ≈ (z_{1-α/2} + z_{1-β})² * σ² / (β² * R²_first * Var(Z))
```

### 4.3 Minimum Detectable Effect (MDE)

Given sample size n and target power:

```
MDE = (z_{1-α/2} + z_{1-β}) * SE_IV
```

### 4.4 Power Curve Parameters

| Effect Size | Required n (80% power) | Required n (90% power) |
|-------------|------------------------|------------------------|
| 0.01 | ~10,000 | ~15,000 |
| 0.02 | ~2,500 | ~3,750 |
| 0.05 | ~400 | ~600 |
| 0.10 | ~100 | ~150 |
| 0.20 | ~25 | ~40 |

*(Assuming R²_first = 0.884, α = 0.05)*

---

## 5. Heterogeneity Analysis

### 5.1 Stratified IV Estimation

Run separate IV analyses by stratum:

1. **By Region:**
   - Core: cells within 100μm of tumor center
   - Margin: cells 100-300μm from center
   - Infiltrating: cells >300μm from center

2. **By Hypoxia:**
   - Normoxic: O2 > 20 mmHg
   - Mildly hypoxic: 10-20 mmHg
   - Severely hypoxic: <10 mmHg

3. **By ecDNA Burden:**
   - Low: ecDNA < 10
   - Medium: 10-30
   - High: >30

### 5.2 Cochran's Q Test for Heterogeneity

```
Q = Σ w_k (β̂_k - β̂_pooled)²
```

Where w_k = 1 / SE²_k

Under H0 (homogeneity): Q ~ χ²(K-1)

### 5.3 I² Statistic

```
I² = max(0, (Q - (K-1)) / Q)
```

| I² | Heterogeneity Level |
|----|---------------------|
| 0-25% | Low |
| 25-50% | Moderate |
| 50-75% | Substantial |
| >75% | Considerable |

### 5.4 Radial Effect Profiles

Estimate effects as smooth function of distance from tumor center:

```
β(r) = f(r) estimated via local regression
```

Kernel bandwidth selected by cross-validation.

---

## 6. Sensitivity Analysis

### 6.1 Rosenbaum Bounds

Assess sensitivity to unobserved confounding.

For treatment effect Γ (odds ratio of hidden bias):
```
P(p_upper > α | Γ) = robustness to confounding at level Γ
```

**Interpretation:**
- Γ = 1: No hidden bias (randomization)
- Γ = 2: Unobserved confounder doubles treatment odds
- Critical Γ*: Smallest Γ where p_upper > α

### 6.2 E-Value

Minimum confounding strength to explain away observed effect:

```
E-value = RR + sqrt(RR * (RR - 1))
```

Where RR is the risk ratio estimate.

**Interpretation:**
- E-value = 3: Confounder must have RR ≥ 3 with both treatment and outcome to explain away effect
- Higher E-value → more robust to confounding

### 6.3 Instrument Exclusion Tests

**Falsification Tests:**
1. Test Z → pre-treatment outcomes (should be null)
2. Test Z → known non-pathway outcomes (should be null)
3. Overidentification test (if multiple instruments)

---

## 7. Causal Discovery

### 7.1 PC Algorithm

Constraint-based causal discovery:

1. Start with complete undirected graph
2. Remove edges based on conditional independence tests
3. Orient edges using v-structures and propagation rules

**Independence Test:**
Partial correlation test with Fisher's z-transformation:
```
z = 0.5 * log((1+r)/(1-r)) * sqrt(n-k-3)
```

### 7.2 GES (Greedy Equivalence Search)

Score-based causal discovery:

1. Forward phase: Add edges that improve BIC score
2. Backward phase: Remove edges that improve BIC score
3. Return equivalence class

**BIC Score:**
```
BIC = -2 * log L + k * log(n)
```

### 7.3 Bootstrap Stability

Edge stability across bootstrap samples:
```
Stability(edge) = #{bootstrap samples with edge} / B
```

Threshold: edges with stability > 0.5 are retained.

---

## 8. Simulation Parameters

Default values below are taken from `causanta/simulate/params/default.json`. The 6mm publication run uses `paper_config_6mm.json`, which differs in domain size (6000 μm), total_hours (300), `O2_blood_mmHg` (40 vs 60), `q_O2_transfer_per_hr` (2 vs 5), `hypoxia_threshold_mmHg` (18 vs 10), and tumor `division_time_mean_hr` (24 vs 36).

### 8.1 Domain Configuration

| Parameter | Default Value | Description |
|-----------|---------------|-------------|
| width_um | 1000 | Domain width (μm) |
| height_um | 1000 | Domain height (μm) |
| env_grid_um | 10 | Environment grid spacing (μm) |

### 8.2 Time Configuration

| Parameter | Default Value | Description |
|-----------|---------------|-------------|
| total_hours | 168 | Simulation duration (hours) |
| dt_diffusion_hr | 0.01 | Diffusion sub-step (hours) |
| output_interval_hr | 1 | Output frequency (hours) |
| burnin_hours | 20 | Burn-in equilibration (hours) |

### 8.3 Cell Type Parameters

| Cell Type | Can Divide | Division Time (hr) | Apoptosis Rate (/hr) | Migration Speed (μm/hr) |
|-----------|------------|--------------------|----------------------|-------------------------|
| Neuron | no | — | 1e-4 | 0 |
| Astrocyte | yes | 168 ± 24 | 1e-4 | 3 |
| Oligodendrocyte | no | — | 1e-4 | 1 |
| Microglia | no | — | 1e-4 | 30 |
| Endothelial | yes | 60 ± 12 | 1e-4 | 10 |
| Pericyte | no | — | 1e-4 | 5 |
| Tumor | yes | 36 ± 8 | 5e-5 | 10 |
| RecruitedImmune | yes | 36 ± 8 | 1e-3 | 35 |
| Necrotic | no | — | 0 | 0 |

### 8.4 ecDNA Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Initial copies | 20 | Starting ecDNA in seed tumor |
| Segregation p | 0.5 | Binomial probability |
| α (effect on division) | 0.30 | T_eff = T_base / (1 + α·log₂(1 + ecDNA)) |
| β (effect on VEGF) | 0.10 | VEGF_eff = VEGF_base × (1 + β·ecDNA) |
| δ (effect on migration) | 0.05 | v_eff = v_base × (1 + δ·ecDNA) |
| γ (effect on survival) | 0.50 | a_eff = a_base / (1 + γ·ecDNA) |

---

## 9. Success Criteria

### 9.1 Publication-Ready Metrics

| Criterion | Target | Measurement |
|-----------|--------|-------------|
| IV Accuracy | <15% bias | |β̂_IV - β_true| / β_true |
| CI Coverage | ≥90% | Fraction containing true value |
| Power | ≥80% | At δ = 0.05, n = 1000 |
| Segregation Test | p > 0.05 | KS test vs Binomial |
| First-stage F | >10 | Standard threshold |
| E-value | >2.0 | Robustness to confounding |

### 9.2 Validation Benchmarks

1. **Segregation Validation:**
   - Mean daughter fraction: 0.50 ± 0.02
   - Variance matches Binomial prediction
   - No trend with parent ecDNA count

2. **First Stage Strength:**
   - F > 100 (very strong by design)
   - R² > 0.80 for ecDNA → EGFR

3. **Effect Recovery:**
   - All four effects (α, β, δ, γ) recovered within 2 SE
   - IV estimates unbiased, OLS estimates show confounding bias

4. **Heterogeneity Detection:**
   - Detect planted heterogeneity when effect varies by region
   - Q-test sensitive to true heterogeneity

---

## 10. Software Implementation

### 10.1 Module Structure

```
causanta/
├── simulate/
│   ├── core.py          # Main simulation loop
│   ├── ecdna.py         # ecDNA segregation model
│   ├── sweep.py         # Parameter sweep system
│   └── params/          # Configuration files
├── analyze/
│   ├── iv.py            # 2SLS estimation + diagnostics
│   ├── bootstrap.py     # BCa confidence intervals
│   ├── power.py         # Power analysis
│   ├── heterogeneity.py # Stratified IV analysis
│   ├── effects.py       # Effect estimation
│   └── loader.py        # Data loading utilities
└── visualize/
    ├── nature_style.py  # Publication styling
    └── causal_figures.py # Figure generation
```

### 10.2 Dependencies

- Python ≥ 3.9
- NumPy ≥ 1.21
- SciPy ≥ 1.7
- Matplotlib ≥ 3.5
- (Optional) causal-learn for PC/GES algorithms

### 10.3 Reproducibility

All simulations use:
- Fixed random seed (default: 42)
- Deterministic NumPy operations where possible
- Version-pinned dependencies

### 10.4 Computational Requirements

| Study Type | Estimated Runtime | Memory |
|------------|------------------|--------|
| Single simulation (200hr) | ~5 min | ~500 MB |
| Baseline validation (20 reps) | ~2 hr | ~1 GB |
| Full parameter sweep | ~8-24 hr | ~2 GB |
| Bootstrap (1000 samples) | ~30 min | ~500 MB |

---

## References

1. Angrist JD, Pischke JS. *Mostly Harmless Econometrics*. Princeton, 2009.
2. Efron B, Tibshirani RJ. *An Introduction to the Bootstrap*. Chapman & Hall, 1993.
3. Rosenbaum PR. *Observational Studies*. Springer, 2002.
4. VanderWeele TJ, Ding P. Sensitivity analysis in observational research. *Ann Intern Med* 2017.
5. Spirtes P, Glymour C, Scheines R. *Causation, Prediction, and Search*. MIT Press, 2000.
6. Staiger D, Stock JH. Instrumental variables regression with weak instruments. *Econometrica* 1997.

---

*CAUSANTA Supplementary Materials v1.0*
*For Nature Methods submission*
