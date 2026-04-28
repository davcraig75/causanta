# Manuscript: Extrachromosomal DNA as a Causal Instrument for Spatial Multi-Omics

This directory contains the manuscript and supplementary materials for the CAUSANTA framework publication.

## Contents

| File | Description |
|------|-------------|
| [`manuscript.md`](manuscript.md) | Main manuscript in Markdown format (GitHub-rendered) |
| [`manuscript.docx`](manuscript.docx) | Main manuscript in Word format (with embedded figures) |
| [`supplementary_materials.docx`](supplementary_materials.docx) | Supplementary materials in Word format |

## Figures

| Figure | Description |
|--------|-------------|
| [`figures/figure1_segregation.png`](figures/figure1_segregation.png) | ecDNA segregation validation (binomial distribution) |
| [`figures/figure2_forest_ols_iv.png`](figures/figure2_forest_ols_iv.png) | Forest plot comparing OLS vs IV estimates |
| [`figures/figure4_sensitivity.png`](figures/figure4_sensitivity.png) | Sensitivity analysis (Rosenbaum bounds, E-values) |
| [`figures/figure5_spatial.png`](figures/figure5_spatial.png) | Spatial heterogeneity of causal effects |
| [`figures/table1_siv_comparison.png`](figures/table1_siv_comparison.png) | Comparison of somatic instrumental variables |

## Viewing the Manuscript

The markdown manuscript is designed to render properly on GitHub. To view with full styling:

1. **On GitHub:** Simply click on [`manuscript.md`](manuscript.md) - figures and tables will render automatically
2. **Locally with VS Code:** Use a Markdown preview extension with CSS support
3. **As PDF:** Use pandoc: `pandoc manuscript.md -o manuscript.pdf --css=assets/style.css`

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

## Authors

- **David W. Craig** - Department of Integrative Translational Science, City of Hope
- **Andrei S. Rodin** - Department of Computational Medicine, City of Hope

*Correspondence: dacraig@coh.org*
