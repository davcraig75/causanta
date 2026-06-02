# Manuscript: Extrachromosomal DNA as a Causal Instrument for Spatial Multi-Omics

This directory contains the manuscript, supporting documents, and figures for the CAUSANTA framework publication.

## Manuscript and Supporting Documents

| File | Description |
|------|-------------|
| [`manuscript.md`](manuscript.md) | Main manuscript in Markdown (GitHub-rendered) |
| [`manuscript.docx`](manuscript.docx) | Main manuscript in Word format (with embedded figures) |
| [`supplementary_materials.md`](supplementary_materials.md) | Supplementary methods in Markdown |
| [`supplementary_materials.docx`](supplementary_materials.docx) | Supplementary methods in Word format |
| [`theory.md`](theory.md) | Theoretical framework + simulation-engine specification |
| [`tutorial.md`](tutorial.md) | Step-by-step causal-inference tutorial |
| [`analysis_guide.md`](analysis_guide.md) | Statistical methods reference |
| [`immune_system.md`](immune_system.md) | Tumor–immune interaction modeling |
| [`figure_legends.docx`](figure_legends.docx) | Editable figure captions |

## Figures

| Figure | Description |
|--------|-------------|
| [`figures/table1_siv_comparison.png`](figures/table1_siv_comparison.png) | Comparison of somatic instrumental variables |
| [`figures/figure1_segregation_firststage.png`](figures/figure1_segregation_firststage.png) | ecDNA segregation validation and first-stage regression |
| [`figures/figure2_ols_vs_iv.png`](figures/figure2_ols_vs_iv.png) | Forest plot comparing OLS vs IV estimates |
| [`figures/figure3_causal_dag.png`](figures/figure3_causal_dag.png) | Structural causal model (DAG) |
| [`figures/figure4_sensitivity.png`](figures/figure4_sensitivity.png) | Sensitivity analysis (Rosenbaum bounds, E-values) |
| [`figures/figure5_spatial_heterogeneity.png`](figures/figure5_spatial_heterogeneity.png) | Spatial heterogeneity of causal effects |
| [`figures/comprehensive_robustness_figure.png`](figures/comprehensive_robustness_figure.png) | Fig 6 — multi-seed robustness (2 scales × 3 scenarios) |
| [`figures/timeseries_figure.png`](figures/timeseries_figure.png) | Fig 7 — time-resolved discovery F1 and OLS bias |
| [`figures/figure8_simulation_architecture.png`](figures/figure8_simulation_architecture.png) | Fig 8 — simulation architecture schematic |
| [`figures/figure_recovery_identifiability.png`](figures/figure_recovery_identifiability.png) | Fig 9 — structural-parameter recovery & identifiability (κ, β, δ, α; γ not identified) |
| [`figures/figure_simulation_overview.png`](figures/figure_simulation_overview.png) | Data-driven simulation overview (2 mm baseline) |
| [`figures/figure_simulation_multiscale.png`](figures/figure_simulation_multiscale.png) | Static export of the multi-scale explainer |
| [`figures/figure_legends.md`](figures/figure_legends.md) | All figure captions (Markdown) |

PDF versions of every figure live alongside the PNGs.

## Interactive Explorers

Open directly in a browser; no server required.

| File | Description |
|------|-------------|
| [`causal_dag_interactive.html`](causal_dag_interactive.html) | Interactive causal DAG (draggable nodes) |
| [`simulation_parameter_explorer.html`](simulation_parameter_explorer.html) | Drag α/β/δ/γ/κ/κ_hyp; live in-browser OLS-vs-2SLS Monte-Carlo |
| [`simulation_multiscale_explorer.html`](simulation_multiscale_explorer.html) | Section→Tissue→Cell reactive simulation diagram (all 9 cell types, O₂/VEGF fields, vessels, immune killing, necrosis, ecDNA→EGFR→effects, time animation) |

## Viewing the Manuscript

The Markdown manuscript renders directly on GitHub. For local PDF generation:

```bash
pandoc manuscript.md -o manuscript.pdf --resource-path=figures/
```

## Citation

```bibtex
@article{craig2024ecdna,
  title = {Extrachromosomal DNA as a Causal Instrument for Spatial Multi-Omics},
  author = {Craig, David W. and Rodin, Andrei S.},
  year = {2024},
  institution = {City of Hope}
}
```

## Authors

- **David W. Craig** — Department of Integrative Translational Science, City of Hope
- **Andrei S. Rodin** — Department of Computational Medicine, City of Hope

*Correspondence: dacraig@coh.org*
