# Figure Legends

## Main Figures

### Figure 1. Validation of ecDNA as an Instrumental Variable

**(A) ecDNA segregation follows a binomial distribution.** Histogram of daughter ecDNA fractions (daughter / replicated total) across 19,534 division events from a 300-hour CAUSANTA simulation on a 6 mm × 6 mm domain (6 mm baseline, seed 43). The distribution is centered on the expected mean of 0.5 (dashed red line), with the observed mean of 0.4997 not significantly different from 0.5 (one-sample t-test: t = −0.52, p = 0.60). The observed variance of 0.0069 is consistent with per-cell binomial variation around the prescribed Binomial(N, 0.5) segregation. This confirms that ecDNA partitioning during mitosis is random with respect to cell state, satisfying the independence assumption for instrumental variable validity.

**(B) First-stage regression demonstrates strong instrument relevance.** Scatter plot of EGFR expression versus ecDNA copy number for tumor cells at the final simulation timepoint (t = 300 h, ≈ 7,200–14,800 tumor cells per seed in the 6 mm baseline). Fitting the structural first-stage model EGFR = X_base + κ·Z + κ_hyp·M·(X_base + κ·Z) + ε recovers κ̂ = 1.23 ± 0.040 across the five seeds, matching the configured κ = 1.21 and the 2.5-fold hypoxic upregulation (1 + κ_hyp = 1 + 1.5). In the `removed` scenario (κ_hyp = 0), where no hypoxia confounding exists, the simpler univariate OLS EGFR ~ const + Z recovers intercept 2.91 and slope 1.21 within sampling error of the configured (X_base = 2.89, κ = 1.21). First-stage F-statistics range from 29,057 to 149,960 across the 30 multi-seed runs — more than three orders of magnitude above the Staiger-Stock weak-instrument threshold of 10.

---

### Figure 2. Instrumental Variable Estimation Corrects Confounding Bias

Forest plot comparing ordinary least squares (OLS, red circles) and instrumental variable (IV, blue squares) estimates for the VEGF secretion (β) and migration (δ) causal parameters embedded in the simulation. Error bars indicate 95% confidence intervals. Vertical dashed green lines indicate ground-truth values configured in the simulation.

OLS estimates of the *linear* slope of VEGF on EGFR are systematically biased upward by hypoxia, which independently affects both EGFR expression (via HIF-2α translational upregulation) and VEGF secretion (via HIF-1α transcriptional upregulation). On the 6 mm baseline cells (κ_hyp = 1.5), the multi-seed pipeline reports OLS slope 4.15 ± 0.10 versus 2SLS slope 3.70 ± 0.22 (OLS bias = +12.5% ± 4.1%, n = 5 seeds). When the structural form is inverted (LHS = VEGF/(VEGF_base·gate(M)) − 1, RHS = β·√EGFR), both OLS and 2SLS recover β = 0.100 ± 0.000 exactly, matching the configured ground truth — because once the hypoxia gating is undone, no confounding remains on the structural scale.

The IV estimator successfully purges confounding bias on the linear scale by isolating variation in EGFR expression attributable solely to the random segregation of ecDNA at cell division. The same pattern, at smaller magnitude, holds for migration (Figure 7 multi-seed reliability table).

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

E-values follow VanderWeele & Ding (*Ann Intern Med* 167:268–274, 2017), E = RR + √(RR(RR − 1)), implemented in `causanta.analyze.iv.compute_e_value` (with the standard RR < 1 → 1/RR flip and the rule that the CI E-value uses the bound closest to the null). Evaluated at a 10-ecDNA-copy contrast (≈ 12 EGFR units for κ = 1.21), the point-estimate E-value is 3.07 for the migration effect δ and 1.34 for the VEGF effect β; the corresponding CI E-values are 3.05 and 1.34. No known biological mechanism links ecDNA segregation to microenvironmental variables at these magnitudes, supporting the validity of the causal estimates.

---

### Figure 5. Spatial causal inference: recovery of a planted gradient and a flat negative control

**(A) Planted spatial zones (2 mm gradient run).** Scatter plot of tumor-cell positions colored by the concentric zone they occupy, defined by distance from the tumor-seed centroid: core (< 780 µm; planted δ = 0.03), margin (780–970 µm; planted δ = 0.07), and infiltrating (> 970 µm; planted δ = 0.05). Dashed circles mark the zone boundaries.

**(B) Per-zone 2SLS estimate of δ for two runs.** Grey bars: the constant-δ control run, where the per-zone structural 2SLS estimate is flat at δ̂ = 0.050 in every zone (the framework does not fabricate a gradient). Red bars: the planted-gradient run, where the identical per-zone 2SLS recovers δ̂ = 0.030 / 0.070 / 0.050 for core / margin / infiltrating, each matching the planted ground truth (black markers) to the third decimal. Error bars are 95% CIs (narrow; n ≈ 2,600–4,400 tumor cells per zone).

Together the two runs show the spatially resolved SIV estimator is both **specific** (returns a uniform map when the effect is spatially constant — a negative control against confounder-driven artifacts) and **sensitive** (recovers genuine effect modification when it is present). The planted gradient is embedded via the simulator's optional `ecDNA_migration_spatial` mode (`modulate_migration_speed`); the constant-δ runs leave it off. This matters for real-tissue application, where an apparent core-versus-margin difference in a naive (non-IV) analysis could reflect the spatial structure of the confounder (hypoxia, density) rather than true effect modification, which the per-zone IV estimator separates.

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

### Simulation Overview Figure (`figure_simulation_overview.png`)

Data-driven snapshot of the v19 simulation (canonical 2 mm baseline run), complementing the schematic architecture figure. **(A)** Spatial tissue at the final timestep, every cell coloured by type — the populated parenchyma (neurons, astrocytes, oligodendrocytes, microglia), the endothelial/pericyte vascular network, the central tumor mass, infiltrating recruited-immune cells, and the necrotic core (final composition: Tumor 11,259; Endothelial 5,297; Pericyte 1,270; Neuron 565; Necrotic 439; RecruitedImmune 393; Astrocyte 311; Oligodendrocyte 152; Microglia 89). **(B)** Oxygen field with the 18 mmHg hypoxia boundary (dashed) and tumor/necrotic overlay, showing the hypoxic tumor core. **(C)** Population dynamics over 240 h (tumor, necrotic, recruited-immune; symlog scale). **(D)** The somatic instrument: ecDNA copy number → EGFR expression in tumor cells, coloured by hypoxia — the two-cohort first-stage signature (normoxic slope ≈ κ = 1.21, hypoxic ≈ κ × 2.5). Generated by `scripts/make_figure_simulation_overview.py`.

### Interactive Causal DAG (`docs/causal_dag_interactive.html`)

Self-contained interactive page (no external dependencies) showing the structural causal model the simulator implements and the analysis is benchmarked against: ecDNA (instrument Z) → EGFR (exposure X) → {VEGF, migration, proliferation, survival}, with hypoxia (confounder M) acting on both EGFR (HIF-2α) and the outcomes. Each edge is annotated with its **configured** parameter and the **recovered** estimate from the v19 runs (κ → κ̂ = 1.23; β, δ → IV β̂ = 0.100, δ̂ = 0.050; α → 0.31 on the κ_hyp = 0 runs; γ not identified; OLS-vs-IV bias +8.7%/+12.5%; first-stage F 29,057–149,960; E-values 3.07/1.34; PC F1 1.000/0.954). Nodes are draggable; clicking a node or edge opens a detail panel; a toggle highlights the confounded path that biases OLS versus the instrument path that 2SLS exploits. Generated by `scripts/make_dag_html.py`.

### Figure 9. Structural-parameter recovery and identifiability (`figure_recovery_identifiability.png`)

Summary of which ecDNA→phenotype effects are recovered from the canonical multi-seed runs and which are not (5 seeds per cell). **Left:** recovered / configured ratio (truth = 1.0, dashed line) with 95% CI for the point-identified coefficients — the first-stage slope κ̂ = 1.23, the structural effects β = 0.100 (VEGF) and δ = 0.050 (migration) recovered by 2SLS in every scenario, and α = 0.315 (division) recovered from inter-division intervals on the κ_hyp = 0 runs (0.323 at 2 mm, 0.307 at 6 mm). γ (survival) is a hatched "not identified" row. **Right:** identifiability class and rationale for each parameter — β and δ are deterministic per-cell phenotypes read off the snapshot; α is a dynamic rate identified only where EGFR is exactly reconstructable from ecDNA (κ_hyp = 0, T_base known); γ leaves no footprint because baseline apoptosis (5×10⁻⁵/hr) produces essentially no deaths, so it is a configured design parameter rather than a fit target. Generated by `scripts/make_figure_recovery_identifiability.py`.

### Multi-Scale Simulation Explainer (`docs/simulation_multiscale_explorer.html` + static `figure_simulation_multiscale.{png,pdf}`)

A self-contained interactive HTML (no external dependencies) that walks the reader through how the CAUSANTA simulation works at three scales, with a static three-panel PNG/PDF companion for copy-paste into slides or print:

- **Section** (≈ 6 mm × 6 mm) — the whole tumor section at the chosen time. The radial O₂ field (hypoxic red core → oxygenated blue rim), the diffusion-limited necrotic core, the native + tumor-induced vessel network, and infiltrating recruited-immune cells are all rendered. A draggable dashed "neighbourhood" rectangle selects the patch shown in the next panel.
- **Tissue** (≈ 400 µm patch) — the contents of the selected neighbourhood, with all nine cell types drawn at biologically-weighted abundances per zone (necrotic core / hypoxic interior / proliferating margin / normal parenchyma). A vessel passes through oxygenated zones; recruited-immune cells form dashed "synapses" with the tumor cells they kill (red ring); each tumor cell carries visible ecDNA specks. Click a tumor cell to inspect it.
- **Cell** (single tumor cell) — the ecDNA → EGFR → {α division, β VEGF, δ migration, γ apoptosis} causal chain instantiated for the selected cell, with the simulator's structural equations evaluated live and shown as four gauges. The hypoxic-vs-normoxic state is annotated and the (1 + κ_hyp · M) confounder term is exposed in the EGFR expression.

Time animates (▶ Play): the tumor grows, the necrotic core forms past t ≈ 140 h once it outstrips diffusion, vessels sprout, and immune infiltration ramps from t ≈ 60 h. Toggleable O₂ field, VEGF field, vessels, immune-killing overlays, and cell-type labels. The static PNG (Figure: Multi-Scale Simulation Explainer) captures a representative t = 240 h frame at the proliferating margin (ecDNA = 24, hypoxic), with the same three-panel layout. Generated by `scripts/make_figure_simulation_multiscale.py` and the hand-authored HTML.

### Interactive Parameter Explorer (`docs/simulation_parameter_explorer.html`)

Self-contained interactive page (no external dependencies). Drag sliders for any structural coefficient (α, β, δ, γ, κ, κ_hyp) or the hypoxic fraction, or click a scenario preset (baseline/reduced/removed), and the structural response curves redraw live: the first-stage EGFR vs ecDNA line, the √-form VEGF curve (β), the log₂-form division-time (α) and apoptosis (γ) curves, and the linear migration curve (δ). A panel runs a fresh in-browser Monte-Carlo population on every drag, computing the naive OLS slope of VEGF on EGFR versus the ecDNA-instrumented 2SLS slope and the resulting bias — demonstrating live how the hypoxia→EGFR confounder inflates OLS at κ_hyp > 0 and how the bias collapses to ≈0 at κ_hyp = 0. Each curve is badged with its identifiability class (β, δ, α identified; γ not identified).

### Multi-scale Simulation Explorer (`docs/simulation_multiscale_explorer.html`; static export `figure_simulation_multiscale.png`)

Self-contained reactive diagram explaining how the agent-based simulation works across three linked spatial scales. **Section (~6 mm):** the whole tumour section grown on the reaction–diffusion microenvironment — a proliferating tumour rim, a diffusion-limited necrotic core, the sprouting vascular network, infiltrating immune cells, and the surrounding normal parenchyma; a draggable magnifier selects a neighbourhood. **Tissue (~400 µm):** the selected neighbourhood, whose composition changes with the magnifier's radial zone (necrotic core / hypoxic interior / proliferating margin / normal parenchyma), showing all nine cell types — neurons, astrocytes, oligodendrocytes, microglia (resident immune), endothelial and pericyte vessel cells, ecDNA-bearing tumour cells, infiltrating recruited-immune cells forming lethal synapses (red dashed lines), and necrotic debris — over the O₂/VEGF fields. **Cell (~20 µm):** a clicked tumour cell, showing its ecDNA copies → EGFR (gene dosage, with the 2.5× hypoxic boost) → the four phenotype rates (α division, β VEGF, δ migration, γ survival). A time slider/Play button grows the section over 0–300 h, and a strip highlights the six per-step simulation phases (environment diffusion → cell behaviours → angiogenesis → immune recruitment/killing → cleanup/necrosis → output). Toggles overlay the O₂ and VEGF fields, vessels, immune cells, and labels. Hand-authored; colours match the manuscript cell-type palette.

---

## Figure Generation

The data figures were generated from CAUSANTA simulation output using `scripts/generate_figures.py` and the per-figure scripts `make_figure1_3panel_seed43.py`, `make_figure2_perseed.py`, `make_figure5_spatial.py`, `make_figure6_robustness.py`, `make_figure7_timeseries.py`, `make_figure_simulation_overview.py`, `make_figure_recovery_identifiability.py`, and `generate_simulation_figure.py`. Raw data are in the `output/` directory; the copy-paste-ready figure set (PNG + PDF) is collected in `figures/manuscript_v21/`.

**Software:** Python 3.10+, matplotlib 3.7+, scipy 1.10+, pandas 2.0+

**Simulation:** CAUSANTA v0.2.1, 30-run multi-seed benchmark (2 mm and 6 mm domains × baseline/reduced/removed × seeds 42–46); canonical first-stage / segregation figure from the 6 mm baseline seed-43 run (300 h). Parameters in `causanta/simulate/params/{large,xlarge,multiseed_*}.json`.
