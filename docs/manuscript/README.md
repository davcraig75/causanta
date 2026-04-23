# Manuscript: Extrachromosomal DNA as a Causal Instrument for Spatial Multi-Omics

This directory contains the manuscript and supplementary materials for the CAUSANTA framework publication.

## Contents

| File | Description |
|------|-------------|
| [`manuscript.md`](manuscript.md) | Main manuscript in Markdown format (GitHub-rendered) |
| [`Causal.v042226.v2.docx`](Causal.v042226.v2.docx) | Original Word document |
| [`supplementary_materials.docx`](supplementary_materials.docx) | Supplementary materials (Word) |

## Figures

| Figure | Description |
|--------|-------------|
| [`figures/figure1_dag.png`](figures/figure1_dag.png) | Causal DAG showing ecDNA as instrumental variable |
| [`figures/figure2_segregation.png`](figures/figure2_segregation.png) | Validation of ecDNA segregation as random process |
| [`figures/figure3_iv_estimates.png`](figures/figure3_iv_estimates.png) | IV vs OLS estimates compared to ground truth |
| [`figures/figure4_spatial.png`](figures/figure4_spatial.png) | Spatial heterogeneity of causal effects |

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
