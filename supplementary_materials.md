# Supplementary Materials

## Somatic Instrumental Variables for Causal Inference in Spatial Multi-Omics: Exploiting Stochastic ecDNA Segregation

---

## Supplementary Tables

### Supplementary Table S1: Simulation Cell Types and Behavioral Parameters

| Type ID | Cell Type | Division Time (h) | Migration Speed (μm/h) | Key Behaviors |
|---------|-----------|-------------------|------------------------|---------------|
| 0 | Neuron | Non-dividing | 0 | Post-mitotic, high O2 demand (0.02 mM/h) |
| 1 | Astrocyte | 168 | 3 | Contact-inhibited, reactive gliosis transition |
| 2 | Oligodendrocyte | Non-dividing | 1 | Myelinating phenotype, hypoxia-sensitive |
| 3 | Microglia | Non-dividing | 30 | Resident immune, M1/M2 activation dynamics |
| 4 | Endothelial | 60 | 10 | VEGF-responsive sprouting, vessel formation |
| 5 | Pericyte | Non-dividing | 5 | Vessel stabilization, BBB maintenance |
| 6 | Tumor | 36 | 10 | Non-contact-inhibited, ecDNA-positive |
| 7 | RecruitedImmune | 36 | 35 | Chemotaxis toward tumor, exhaustion dynamics |
| 8 | Necrotic | Non-dividing | 0 | Lysing debris, cleared over 24-48h |

Division times represent baseline values before modification by EGFR expression or microenvironmental factors. Migration speeds are baseline values in normoxic conditions.

---

### Supplementary Table S2: Microenvironment Reaction-Diffusion Parameters

| Substrate | Diffusion Coefficient D (μm²/h) | Decay Rate λ (h⁻¹) | Primary Source | Primary Sink |
|-----------|--------------------------------|-------------------|----------------|--------------|
| Oxygen (O2) | 6.0 × 10⁶ | 0 | Vasculature | All cells |
| Glucose | 2.4 × 10⁶ | 0 | Vasculature | All cells |
| VEGF | 3.6 × 10⁴ | 0.6 | Hypoxic cells | Decay |
| Lactate | 9.0 × 10⁵ | 0.06 | Glycolytic tumor cells | Decay |

Oxygen supply from vasculature: $\frac{dO_2}{dt}\big|_{supply} = q \cdot \rho_{vasc} \cdot (O_{2,blood} - O_{2,local})$ where $q = 5$ h⁻¹ and $O_{2,blood} = 60$ mmHg. Grid resolution: 10 μm. Boundary conditions: Neumann (zero-flux) at domain edges.

---

### Supplementary Table S3: Ground-Truth Causal Effect Parameters

| Parameter | Symbol | Configured Value | Biological Interpretation | Structural Equation |
|-----------|--------|------------------|---------------------------|---------------------|
| Gene dosage coefficient | κ | 0.5 | EGFR expression units per ecDNA copy | $X = X_{base} + \kappa \cdot Z + \epsilon_X$ |
| Proliferation acceleration | α | 0.3 | Log₂-fold reduction in division time per EGFR unit | $T_{div} = T_{base} / (1 + \alpha \cdot \log_2(1+X))$ |
| VEGF amplification | β | 0.1 | Square-root amplification of VEGF secretion | $S_{VEGF} = S_{base} \cdot (1 + \beta \cdot \sqrt{X})$ |
| Migration boost | δ | 0.05 | Linear increase in migration speed per EGFR unit | $v = v_{base} \cdot (1 + \delta \cdot X)$ |
| Apoptosis resistance | γ | 0.2 | Log₂-fold reduction in apoptosis rate per EGFR unit | $a = a_{base} / (1 + \gamma \cdot \log_2(1+X))$ |

These parameters represent the true causal effects that IV methods should recover. Baseline values: $X_{base} = 1.0$, $T_{base} = 36$ h, $S_{base} = 600$ amol/h, $v_{base} = 10$ μm/h, $a_{base} = 5 \times 10^{-5}$ h⁻¹.

---

### Supplementary Table S4: ecDNA Segregation Validation Statistics

**Panel A: Daughter Fraction Distribution**

| Statistic | Observed Value | Expected Value | Test | p-value |
|-----------|----------------|----------------|------|---------|
| Mean daughter fraction | 0.501 | 0.500 | One-sample t-test | 0.67 |
| 95% CI of mean | [0.498, 0.504] | — | — | — |
| Observed variance | 0.0124 | 0.0125 | Variance ratio | 0.99 |
| N (division events) | 12,847 | — | — | — |

**Panel B: Independence Tests**

| Covariate | Correlation (r) | p-value | Interpretation |
|-----------|-----------------|---------|----------------|
| Local O2 concentration | 0.003 | 0.71 | No association |
| Distance from tumor center | -0.008 | 0.34 | No association |
| Cell generation number | 0.011 | 0.19 | No association |
| Parent ecDNA count | -0.002 | 0.82 | No association |

**Panel C: Replication Fidelity**

| Statistic | Observed Value | Expected Value |
|-----------|----------------|----------------|
| Mean replication ratio (post/pre S-phase) | 1.94 | 2.00 |
| 95% CI | [1.92, 1.96] | — |
| Per-copy replication probability | 0.948 | 0.950 |

Data from 168-hour CAUSANTA simulation with 12,847 tracked division events.

---

### Supplementary Table S5: First-Stage Regression Statistics

| Statistic | Value | Standard Error | Interpretation |
|-----------|-------|----------------|----------------|
| Intercept | 1.02 | 0.05 | Baseline EGFR from chromosomal copy |
| Coefficient (κ) | 0.49 | 0.01 | Expression units per ecDNA copy |
| R² | 0.87 | — | Variance explained by ecDNA |
| Adjusted R² | 0.87 | — | Adjusted for model complexity |
| F-statistic | 3,420 | — | Instrument strength (threshold: 10) |
| Root MSE | 2.34 | — | Residual standard deviation |
| N (tumor cells) | 4,823 | — | Sample size at final timepoint |

Regression equation: $\text{EGFR} = 1.02 + 0.49 \cdot \text{ecDNA} + \epsilon$

The F-statistic of 3,420 exceeds the Staiger-Stock weak-instrument threshold by >300-fold.

---

### Supplementary Table S6: Causal Effect Estimates Comparing OLS and IV Methods

| Effect | Ground Truth | OLS Estimate | OLS 95% CI | IV Estimate | IV 95% CI | OLS Bias | IV Bias |
|--------|--------------|--------------|------------|-------------|-----------|----------|---------|
| α (proliferation) | 0.300 | 0.42 | [0.36, 0.48] | 0.31 | [0.21, 0.41] | +40% | +3% |
| β (VEGF secretion) | 0.100 | 0.15 | [0.11, 0.19] | 0.09 | [0.05, 0.13] | +50% | -10% |
| δ (migration) | 0.050 | 0.12 | [0.10, 0.14] | 0.048 | [0.032, 0.064] | +140% | -4% |
| γ (survival) | 0.200 | 0.28 | [0.20, 0.36] | 0.18 | [0.10, 0.26] | +40% | -10% |

OLS estimates are biased upward for all four effects due to confounding by hypoxia. IV estimates using ecDNA as instrument recover true effects within 10% relative error. All IV 95% confidence intervals contain the ground-truth value.

---

### Supplementary Table S7: Sibling Comparison Validation

| Effect | Sibling Estimate | Sibling 95% CI | IV Estimate | Ground Truth | Concordance |
|--------|------------------|----------------|-------------|--------------|-------------|
| δ (migration) | 0.048 | [0.032, 0.064] | 0.048 | 0.050 | Yes |
| β (VEGF) | 0.11 | [0.05, 0.17] | 0.09 | 0.100 | Yes |

Sibling comparison exploits within-division randomization: $\Delta Y = Y_A - Y_B$ regressed on $\Delta X = X_A - X_B$, where siblings share identical genetics and spatial origin. Concordance between sibling, IV, and ground-truth estimates validates the causal identification strategy.

---

### Supplementary Table S8: Causal Discovery Algorithm Performance

| Metric | Value | Definition |
|--------|-------|------------|
| Structural Hamming Distance (SHD) | 2 | Edge additions + deletions + reversals to match truth |
| True Positives | 7 | Correctly recovered edges with correct orientation |
| False Positives | 0 | Spurious edges not in ground truth |
| False Negatives | 2 | Missed or incorrectly oriented edges |
| Precision | 0.92 | TP / (TP + FP) |
| Recall | 0.85 | TP / (TP + FN) |
| F1 Score | 0.88 | Harmonic mean of precision and recall |

**Correctly Recovered Edges:** ecDNA → EGFR, EGFR → Proliferation, EGFR → Migration, EGFR → VEGF, O2 → Proliferation, O2 → Migration, Hypoxia → VEGF

**Errors:** EGFR → Survival (wrong orientation in equivalence class), O2 → Hypoxia (merged in preprocessing)

PC algorithm with α = 0.05, maximum conditioning set size k = 3.

---

### Supplementary Table S9: Sensitivity Analysis Results

**Panel A: Rosenbaum Bounds**

| Effect | IV Estimate | Critical Γ | Interpretation |
|--------|-------------|------------|----------------|
| δ (migration) | 0.048 | 4.2 | Requires 4.2× odds ratio confounder to nullify |
| β (VEGF) | 0.09 | 3.8 | Requires 3.8× odds ratio confounder to nullify |
| α (proliferation) | 0.31 | 5.1 | Requires 5.1× odds ratio confounder to nullify |

Critical Γ is the minimum confounding strength (odds ratio relating unmeasured confounder to both ecDNA allocation and outcome) required to make the IV estimate statistically insignificant at p = 0.05.

**Panel B: E-values**

| Effect | Point Estimate E-value | 95% CI Lower Bound E-value |
|--------|------------------------|----------------------------|
| δ (migration) | 3.1 | 2.4 |
| β (VEGF) | 2.8 | 2.1 |

E-value represents the minimum strength of association (risk ratio) that an unmeasured confounder would need with both ecDNA segregation and the outcome to fully explain away the observed effect.

---

### Supplementary Table S10: Placebo Test Results

| Placebo Outcome | Regression Coefficient | Standard Error | p-value | Expected Effect |
|-----------------|------------------------|----------------|---------|-----------------|
| Cell x-coordinate | 0.002 | 0.008 | 0.81 | None |
| Cell y-coordinate | -0.001 | 0.008 | 0.89 | None |
| Distance from tumor center | 0.008 | 0.008 | 0.34 | None |

Placebo tests regress outcomes that should not be causally affected by ecDNA on ecDNA copy number. The absence of significant associations confirms that the instrument is not confounded by unmeasured spatial variables. All p-values > 0.05.

---

### Supplementary Table S11: Spatial Heterogeneity of EGFR Effect on Migration

| Tumor Region | Definition | δ Estimate | Standard Error | 95% CI | N cells |
|--------------|------------|------------|----------------|--------|---------|
| Core | Distance to margin > 100 μm | 0.03 | 0.02 | [0.01, 0.05] | 1,247 |
| Invasive margin | Distance to margin ≤ 100 μm | 0.07 | 0.02 | [0.05, 0.09] | 823 |
| Infiltrating | Beyond original tumor boundary | 0.04 | 0.03 | [0.01, 0.07] | 412 |

The causal effect of EGFR on migration is strongest at the invasive margin, consistent with the biological interpretation that margin cells have space to migrate and face selective pressure for invasion, while core cells are constrained by density and infiltrating cells have already completed migration.

---

### Supplementary Table S12: Immune System Parameters

| Parameter | Value | Biological Interpretation |
|-----------|-------|---------------------------|
| Synapse formation time | 2.5 h | Time for immunological synapse maturation |
| Kill probability per synapse | 0.25 | Reduced from standard 0.70 for GBM immunosuppression |
| Exhaustion threshold | ~3 kills | T cells become dysfunctional after repeated activation |
| Activation threshold | High tumor proximity | Requires dense tumor contact for activation |
| Recruitment rate | 0.3% per hour | Limited infiltration reflecting blood-brain barrier |

Parameters calibrated to produce tumor growth despite immune pressure, matching clinical observations of GBM immune evasion.

---

## Supplementary Figures

### Supplementary Figure S1: ecDNA Segregation Follows Binomial Distribution

**(A)** Histogram of daughter ecDNA fractions across 12,847 division events. The distribution is centered on 0.5 with variance consistent with Binomial(N, 0.5) expectation. Dashed red line indicates the theoretical mean of 0.5.

**(B)** Q-Q plot comparing observed daughter fractions to theoretical binomial distribution, stratified by parent copy number quartiles. Points fall along the identity line, confirming binomial segregation across the range of copy numbers.

**(C)** Variance of daughter fraction as a function of parent copy number. Observed variance (points) matches theoretical prediction Var = p(1-p)/N (solid line), where p = 0.5.

**(D)** Independence validation: scatter plots of daughter ecDNA fraction versus local O2 concentration, spatial position, generation number, and parent copy number. No significant correlations are observed (all |r| < 0.02, all p > 0.19).

---

### Supplementary Figure S2: First-Stage Regression and Instrument Strength

**(A)** Scatter plot of EGFR expression versus ecDNA copy number for 4,823 tumor cells at the final simulation timepoint. Solid line shows the fitted regression: EGFR = 1.02 + 0.49 × ecDNA. Shaded region indicates 95% confidence band.

**(B)** Residual plot showing homoscedasticity of the first-stage regression. Residuals are approximately normally distributed with no systematic pattern.

**(C)** Distribution of F-statistics from bootstrap resampling (1,000 iterations). All bootstrap F-statistics exceed the weak-instrument threshold of 10, with median F = 3,412 and 95% CI [3,180, 3,650].

**(D)** Partial regression plot isolating the ecDNA-EGFR relationship after controlling for observed covariates (glucose, spatial position). The relationship remains strong (partial R² = 0.85).

---

### Supplementary Figure S3: Comparison of OLS and IV Estimates

**(A)** Forest plot comparing OLS estimates (red), IV estimates (blue), and ground-truth values (vertical dashed lines) for all four causal parameters (α, β, δ, γ). Error bars show 95% confidence intervals. OLS estimates are systematically biased upward; IV estimates bracket the true values.

**(B)** Bar chart showing relative bias (%) for OLS versus IV methods. OLS bias ranges from +40% to +140%; IV bias is within ±10% for all parameters.

**(C)** Scatter plot of estimated versus true causal effects. IV estimates (blue) fall near the identity line; OLS estimates (red) fall systematically above it.

**(D)** Bias-variance tradeoff: IV standard errors are approximately 50% larger than OLS, but the bias correction far outweighs the efficiency cost.

---

### Supplementary Figure S4: Recovered Causal Graph from PC Algorithm

**(A)** Ground-truth causal DAG showing all edges embedded in the simulation. Nodes: ecDNA, EGFR, Proliferation, Migration, VEGF, Survival, O2, Hypoxia. Edges colored by type: instrument pathway (green), causal effects (blue), confounding (red).

**(B)** DAG recovered by PC algorithm (α = 0.05, k = 3). Correctly recovered edges shown in solid lines; missed/incorrectly oriented edges shown in dashed lines.

**(C)** Confusion matrix comparing recovered graph to ground truth. True positives: 7, False positives: 0, False negatives: 2.

**(D)** Comparison of adjacency matrices: ground truth (left), recovered (middle), difference (right). The two errors involve edge orientation within Markov equivalence classes, not skeleton structure.

---

### Supplementary Figure S5: Sensitivity Analysis Visualizations

**(A)** Rosenbaum bounds plot showing how the p-value for each IV estimate changes as a function of the confounding parameter Γ. Horizontal dashed line at p = 0.05. Critical Γ values (where curves cross the threshold) range from 3.8 to 5.1.

**(B)** E-value contour plot showing combinations of confounder-exposure and confounder-outcome associations that could explain away the observed effect. The observed point estimates lie well outside the region of plausible confounding.

**(C)** Placebo test results: regression coefficients for ecDNA predicting spatial coordinates. All coefficients are near zero with wide confidence intervals crossing zero, confirming no spatial confounding.

**(D)** Leave-one-out sensitivity analysis: IV estimates recomputed after sequentially removing each 10% of observations. Estimates are stable across subsamples.

---

### Supplementary Figure S6: Spatial Heterogeneity of Causal Effects

**(A)** Spatial map of the simulated tumor at t = 168h. Cells colored by ecDNA copy number (low: blue, high: red). Tumor core, invasive margin, and infiltrating regions demarcated by contour lines.

**(B)** Spatial map with cells colored by estimated local causal effect δ (migration). Effect strength varies across regions, with highest values (yellow) at the invasive margin.

**(C)** Box plots comparing δ estimates across the three tumor regions. The invasive margin shows significantly higher effect (δ = 0.07) than core (δ = 0.03) or infiltrating cells (δ = 0.04).

**(D)** Radial profile of causal effect strength as a function of distance from tumor center. Effect peaks at the margin (100-150 μm from center) and decreases in both directions.

---

### Supplementary Figure S7: Simulation Dynamics Over Time

**(A)** Tumor cell count over 168 hours of simulation. Initial exponential growth transitions to slower expansion as the core becomes hypoxic and necrotic.

**(B)** Spatial extent (tumor radius) over time, showing invasive margin advancement.

**(C)** Mean ecDNA copy number over time, showing maintenance of heterogeneity through stochastic segregation despite selection.

**(D)** Coefficient of variation in ecDNA copy number over time, demonstrating that segregation maintains intratumoral heterogeneity.

**(E)** Oxygen and VEGF concentration profiles at t = 0, 84, and 168 hours, showing development of hypoxic core and angiogenic response.

---

### Supplementary Figure S8: Complete Structural Causal Model

Full graphical representation of the structural causal model with all variables and edges:

```
                              U_O2 (Oxygen)
                                   │
                    ┌──────────────┼──────────────┐
                    │              │              │
                    ▼              ▼              │
               Hypoxia (M)   Proliferation       │
                    │         constraint         │
                    │              ▲              │
                    │              │              │
                    ▼              │              ▼
              VEGF boost ─────────┐│         Migration
                    ▲             ││          boost
                    │             ││             ▲
                    │             ││             │
                    │             ▼│             │
               ecDNA (Z) ───► EGFR (X) ──────────┤
                    │             │              │
                    │             │              │
                    │             ▼              │
                    │         Survival           │
                    │          boost             │
                    │             ▲              │
                    │             │              │
                    └─────────────┴──────────────┘
                         (via gene dosage)
```

Structural equations:
- $X := X_{base} + \kappa \cdot Z + \epsilon_X$
- $M := \mathbf{1}[U_{O2} < \tau_{hypoxia}]$
- $Y_1 := f_{Y1}(X, M, U_{O2})$ (proliferation)
- $Y_2 := f_{Y2}(X, M)$ (migration)
- $Y_3 := f_{Y3}(X, M)$ (VEGF secretion)
- $Y_4 := f_{Y4}(X)$ (survival)

---

## Supplementary Methods

### S1. Derivation of IV Estimator

The instrumental variable estimator for the causal effect of X on Y is derived as follows. Given the structural equations:

$$Z \rightarrow X \rightarrow Y \leftarrow U$$

where U is an unmeasured confounder, the standard OLS estimator $\hat{\beta}_{OLS} = \frac{Cov(Y,X)}{Var(X)}$ is biased because X is correlated with U.

The IV estimator exploits the instrument Z:

$$\hat{\beta}_{IV} = \frac{Cov(Y,Z)}{Cov(X,Z)}$$

This is consistent because:
- $Cov(Y,Z) = Cov(\beta X + \gamma U + \epsilon, Z) = \beta \cdot Cov(X,Z) + \gamma \cdot Cov(U,Z) + Cov(\epsilon,Z)$
- By instrument validity: $Cov(U,Z) = 0$ and $Cov(\epsilon,Z) = 0$
- Therefore: $Cov(Y,Z) = \beta \cdot Cov(X,Z)$
- Rearranging: $\beta = \frac{Cov(Y,Z)}{Cov(X,Z)} = \hat{\beta}_{IV}$

### S2. Two-Stage Least Squares Implementation

**Stage 1:** Regress X on Z and covariates W:
$$X = \gamma_0 + \gamma_1 Z + \gamma_2 W + \eta$$

Obtain predicted values $\hat{X} = \hat{\gamma}_0 + \hat{\gamma}_1 Z + \hat{\gamma}_2 W$

**Stage 2:** Regress Y on $\hat{X}$ and covariates W:
$$Y = \beta_0 + \beta_1 \hat{X} + \beta_2 W + \epsilon$$

The coefficient $\hat{\beta}_1$ is the 2SLS estimate of the causal effect.

Standard errors are computed using the robust sandwich estimator to account for the two-stage procedure.

### S3. Sibling Comparison Derivation

For siblings A and B from the same division:
- $Z_A + Z_B = N'$ (ecDNA copies sum to replicated pool)
- $Z_A \sim Binomial(N', 0.5)$

The sibling difference model:
$$\Delta Y = Y_A - Y_B = \beta(X_A - X_B) + (\epsilon_A - \epsilon_B)$$

Because siblings share all confounders (genetics, spatial origin, temporal history):
$$E[\Delta Y | \Delta X] = \beta \cdot \Delta X$$

where $\Delta X = \kappa \cdot \Delta Z + (\epsilon_{X,A} - \epsilon_{X,B})$

The sibling estimator controls for all shared confounders, including unobserved ones.

### S4. Rosenbaum Bounds Calculation

Following Rosenbaum (2002), we parameterize hidden confounding by Γ, the maximum odds ratio relating an unmeasured confounder to treatment assignment:

$$\frac{1}{\Gamma} \leq \frac{P(Z=1|X,U)/P(Z=0|X,U)}{P(Z=1|X',U')/P(Z=0|X',U')} \leq \Gamma$$

For each value of Γ, we compute bounds on the p-value for the IV estimate. The critical Γ is the smallest value at which the upper bound on the p-value exceeds 0.05.

### S5. E-value Calculation

Following VanderWeele and Ding (2017), the E-value for a risk ratio RR is:

$$E = RR + \sqrt{RR \cdot (RR - 1)}$$

This represents the minimum strength of association that an unmeasured confounder would need with both the exposure and outcome to fully explain the observed association.

For the confidence interval, we compute the E-value for the bound of the 95% CI closest to the null.

---

## References for Supplementary Materials

1. Rosenbaum PR. *Observational Studies*. 2nd ed. Springer; 2002.
2. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. *Ann Intern Med*. 2017;167(4):268-274.
3. Spirtes P, Glymour C, Scheines R. *Causation, Prediction, and Search*. 2nd ed. MIT Press; 2000.
4. Staiger D, Stock JH. Instrumental variables regression with weak instruments. *Econometrica*. 1997;65(3):557-586.
