# Figure Legends

## Main Figures

### Figure 1. Validation of ecDNA as an Instrumental Variable

**(A) ecDNA segregation follows a binomial distribution.** Histogram of daughter ecDNA fractions across 51,658 division events from a 300-hour CAUSANTA simulation on a 6mm × 6mm domain. The distribution is centered on the expected mean of 0.5 (dashed red line), with the observed mean of 0.501 (solid blue line) not significantly different from 0.5 (bootstrap test p = 0.48). The variance of 0.0071 matches the theoretical prediction for Binomial(N, 0.5) segregation. This confirms that ecDNA partitioning during mitosis is random with respect to cell state, satisfying the independence assumption for instrumental variable validity.

**(B) First-stage regression demonstrates strong instrument relevance.** Scatter plot of EGFR expression versus ecDNA copy number for 44,400 tumor cells at the final simulation timepoint (t = 300h). The fitted regression line (EGFR = 2.89 + 1.21 × ecDNA for normoxic cells; 2.5× higher for hypoxic cells) explains 85.1% of the variance in EGFR expression (R² = 0.851). The first-stage F-statistic of 253,699 exceeds the Staiger-Stock weak-instrument threshold of 10 by more than 25,000-fold, confirming that ecDNA copy number is an exceptionally strong instrument for EGFR expression. This strong relationship arises from the direct gene dosage mechanism: each ecDNA circle contains the EGFR locus, and transcription scales linearly with template number. Note: The remaining variance reflects HIF-1α-mediated EGFR upregulation under hypoxia, which adds variance not explained by ecDNA alone.

---

### Figure 2. Instrumental Variable Estimation Corrects Confounding Bias

Forest plot comparing ordinary least squares (OLS, red circles) and instrumental variable (IV, blue squares) estimates for the VEGF secretion (β) and migration (δ) causal parameters embedded in the simulation. Error bars indicate 95% confidence intervals. Vertical dashed green lines indicate ground-truth values configured in the simulation.

OLS estimates are systematically biased upward due to confounding by hypoxia, which independently affects both EGFR expression (via HIF-1α) and cellular phenotypes. The bias is most severe for the VEGF effect β, where OLS overestimates the true effect by 117% (OLS: 0.22 vs. ground truth: 0.10). This reflects the strong hypoxia-VEGF coupling: HIF-1α independently upregulates both EGFR transcription and VEGF secretion, creating a large spurious association. For migration (δ), OLS overestimates by 27% (OLS: 0.063 vs. ground truth: 0.05).

In contrast, IV estimates using ecDNA as an instrument recover the true causal effects within 10% relative error. Both IV estimates fall within 9% of ground-truth values, and the 95% confidence intervals contain the true parameters. The IV estimator successfully purges confounding bias by isolating variation in EGFR expression attributable solely to the random segregation of ecDNA at cell division.

Parameters: β = VEGF secretion amplification, δ = migration speed boost. Additional parameters (α = proliferation, γ = survival) require division time and apoptosis tracking for estimation.

---

### Figure 3. Causal Directed Acyclic Graph (DAG) for ecDNA-Driven Tumor Biology

Graphical representation of the structural causal model underlying the somatic instrumental variable framework. Nodes represent measured and unmeasured variables; directed edges represent causal relationships.

**Instrument pathway (green):** ecDNA copy number (Z) directly determines EGFR expression (X) through gene dosage. This relationship is the basis for instrument relevance.

**Causal effects (blue):** EGFR expression causally affects four downstream phenotypes: proliferation rate, migration speed, VEGF secretion, and cell survival. These are the causal parameters (α, β, δ, γ) that IV methods aim to estimate.

**Confounding paths (red):** Local oxygen concentration (O2) determines hypoxia status, which independently affects proliferation (via cell cycle arrest), migration (via Go-or-Grow switch), and VEGF secretion (via HIF-1α stabilization). These pathways create confounding that biases naive regression estimates. Critically, oxygen cannot influence ecDNA segregation because partitioning occurs via physical mechanisms during cytokinesis that are independent of the microenvironment.

The key identification insight: ecDNA segregation is random (Binomial(N, 0.5)) conditional on parent copy number, breaking the backdoor path through confounders and enabling unbiased causal effect estimation.

---

### Figure 4. Sensitivity Analysis for Unmeasured Confounding

**(A) Rosenbaum bounds.** Plot showing how the upper-bound p-value for each IV estimate changes as a function of the confounding parameter Γ. The critical Γ (Γ*) is the smallest value at which the p-value upper bound exceeds 0.05 (horizontal red line), representing the minimum confounding strength required to explain away the observed effect.

Strikingly, for both the migration effect δ (blue) and VEGF effect β (orange), the critical Γ* exceeds 30 (off-scale on the plot). This means an unmeasured confounder would need to increase both ecDNA allocation and the outcome by more than 30-fold to explain away the observed effects. This extreme robustness reflects the fundamental randomness of ecDNA segregation—because partitioning occurs via physical mechanisms during cytokinesis that are completely independent of the microenvironment, no biological confounder could plausibly achieve such magnitudes.

**(B) E-values.** Bar chart showing E-values for point estimates and confidence interval bounds. The E-value represents the minimum strength of association (risk ratio) that an unmeasured confounder would need with both ecDNA and the outcome to fully explain away the observed effect.

For the migration effect δ, the point estimate E-value is 4.6 and the CI bound E-value is 4.5. For the VEGF effect β, values are 13.4 and 11.7 respectively—these larger values reflect the stronger IV estimate for VEGF secretion. No known biological mechanism links ecDNA segregation to microenvironmental variables at these magnitudes, supporting the validity of the causal estimates.

---

### Figure 5. Spatial Heterogeneity of Causal Effects

**(A) Spatial distribution of ecDNA copy number.** Scatter plot showing the spatial positions of 8,221 tumor cells at the final simulation timepoint (t = 300h) on a 6mm × 6mm domain. Cells are colored by ecDNA copy number (color scale: blue = low, red = high). The tumor originated from three seed cells at positions (3000, 3000), (2700, 2700), and (3300, 3300) and expanded outward. ecDNA copy number varies substantially across the tumor population due to stochastic segregation, creating the within-tumor variation that enables causal inference.

**(B) Regional variation in EGFR effect on migration.** Bar chart comparing the estimated causal effect of EGFR on migration (δ) across three tumor regions. The invasive margin shows the strongest effect, approximately 2-fold higher than the tumor core. Infiltrating cells beyond the original tumor boundary show an intermediate effect.

This spatial heterogeneity has biological plausibility: core cells are constrained by high cell density and contact inhibition, limiting the phenotypic impact of EGFR-driven motility programs; margin cells have space to migrate and face selective pressure for invasion, amplifying the EGFR effect; infiltrating cells have already completed migration, and EGFR may be less rate-limiting once invasion has occurred.

The dashed red line indicates the overall δ estimate of 0.05. This pattern would be masked in bulk analyses but is revealed by spatially resolved causal inference, suggesting that EGFR-targeted therapy might have the greatest impact on invasion when delivered to the tumor margin.

---

## Supplementary Figures

### Supplementary Figure S1. ecDNA Segregation Validation

**(A)** Histogram of daughter ecDNA fractions across all division events. **(B)** Q-Q plot comparing observed daughter fractions to theoretical binomial distribution. **(C)** Variance of daughter fraction as a function of parent copy number, compared to theoretical prediction. **(D)** Scatter plots of daughter ecDNA fraction versus potential confounders (O2, spatial position, generation number), demonstrating independence.

### Supplementary Figure S2. First-Stage Regression Diagnostics

**(A)** Scatter plot with regression line and 95% confidence band. **(B)** Residual plot showing homoscedasticity. **(C)** Bootstrap distribution of F-statistics. **(D)** Partial regression plot controlling for observed covariates.

### Supplementary Figure S3. Complete OLS vs IV Comparison

**(A)** Forest plot with all four parameters. **(B)** Relative bias bar chart. **(C)** Estimated vs true scatter plot. **(D)** Bias-variance tradeoff visualization.

### Supplementary Figure S4. Causal Discovery Results

**(A)** Ground-truth DAG with edge coloring. **(B)** DAG recovered by PC algorithm. **(C)** Confusion matrix. **(D)** Adjacency matrix comparison.

### Supplementary Figure S5. Sensitivity Analysis Details

**(A)** Full Rosenbaum bounds curves. **(B)** E-value contour plot. **(C)** Placebo test coefficients. **(D)** Leave-one-out stability analysis.

### Supplementary Figure S6. Spatial Analysis

**(A)** Spatial map colored by ecDNA. **(B)** Spatial map colored by local effect estimate. **(C)** Regional box plots. **(D)** Radial profile of effect strength.

### Supplementary Figure S7. Simulation Dynamics

**(A)** Tumor cell count over time. **(B)** Tumor radius over time. **(C)** Mean ecDNA copy number over time. **(D)** Coefficient of variation in ecDNA over time. **(E)** Oxygen and VEGF concentration profiles at three timepoints.

### Supplementary Figure S8. Structural Causal Model

Complete graphical representation of the structural causal model with all variables, edges, and structural equations annotated.

### Supplementary Figure S9. Confounding Structure

Visualization of the confounded relationship between EGFR expression and migration, illustrating why naive regression cannot distinguish causal effects from spurious associations.

### Supplementary Figure S10. Immunological Synapse Timeline

Schematic of the temporal progression from T cell receptor recognition through synapse formation, cytotoxic activity, and exhaustion, with parameter values annotated.

---

## Figure Generation

All figures were generated from CAUSANTA simulation output using the script `scripts/generate_figures.py`. Simulation parameters are specified in `causanta/simulate/params/default.json`. Raw data are available in the `output/` directory.

**Software:** Python 3.10+, matplotlib 3.7+, scipy 1.10+, pandas 2.0+

**Simulation:** CAUSANTA v0.2.0, 300-hour simulation, seed 42, 6mm × 6mm domain, paper_config_6mm.json parameters
