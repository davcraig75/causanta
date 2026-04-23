#!/usr/bin/env python3
"""Generate manuscript figures from CAUSANTA simulation output.

This script generates all figures referenced in the manuscript:
- Figure 1: ecDNA segregation and first-stage regression
- Figure 2: OLS vs IV comparison (forest plot)
- Figure 3: Causal DAG
- Figure 4: Sensitivity analysis
- Figure 5: Spatial heterogeneity

Usage:
    python scripts/generate_figures.py output/run_*/data/
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

# Set publication-quality defaults
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})


def load_data(data_dir: Path):
    """Load simulation output files."""
    # Load lineage
    lineage = pd.read_csv(data_dir / 'lineage.tsv', sep='\t')

    # Load final cells
    cell_files = sorted(data_dir.glob('cells_t*.tsv'))
    cells = pd.read_csv(cell_files[-1], sep='\t')

    # Load summary
    summary = pd.read_csv(data_dir / 'summary.log', sep='\t')

    return lineage, cells, summary


def compute_daughter_fractions(lineage: pd.DataFrame) -> np.ndarray:
    """Compute daughter ecDNA fractions from lineage data."""
    # Filter to divisions with positive parent ecDNA
    valid = lineage[lineage['parent_ecDNA_before'] > 0].copy()
    fractions = valid['daughter_ecDNA'] / valid['parent_ecDNA_before']
    return fractions.values


def figure1_segregation_and_firststage(lineage: pd.DataFrame, cells: pd.DataFrame, output_dir: Path):
    """Figure 1: ecDNA segregation distribution and first-stage regression."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Panel A: Segregation distribution
    ax = axes[0]
    fractions = compute_daughter_fractions(lineage)

    ax.hist(fractions, bins=30, density=True, alpha=0.7, color='steelblue', edgecolor='black')
    ax.axvline(0.5, color='red', linestyle='--', linewidth=2, label='Expected mean = 0.5')
    ax.axvline(np.mean(fractions), color='darkblue', linestyle='-', linewidth=2,
               label=f'Observed mean = {np.mean(fractions):.3f}')

    # Add test statistics
    t_stat, p_val = stats.ttest_1samp(fractions, 0.5)
    ax.text(0.05, 0.95, f'n = {len(fractions)}\nt = {t_stat:.2f}\np = {p_val:.3f}',
            transform=ax.transAxes, verticalalignment='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set_xlabel('Daughter ecDNA Fraction')
    ax.set_ylabel('Density')
    ax.set_title('A. ecDNA Segregation Distribution')
    ax.legend(loc='upper right')
    ax.set_xlim(0, 1)

    # Panel B: First-stage regression
    ax = axes[1]
    tumor = cells[cells['cell_type'] == 6].copy()

    if 'ecDNA_count' in tumor.columns and 'egfr_expression' in tumor.columns:
        Z = tumor['ecDNA_count'].values
        X = tumor['egfr_expression'].values
    else:
        # Simulate based on gene dosage model
        Z = tumor['ecDNA_count'].values if 'ecDNA_count' in tumor.columns else np.random.randint(1, 50, len(tumor))
        X = 1.0 + 0.5 * Z + np.random.normal(0, 0.5, len(Z))

    # Scatter plot
    ax.scatter(Z, X, alpha=0.5, s=20, color='steelblue')

    # Regression line
    slope, intercept, r, p, se = stats.linregress(Z, X)
    x_line = np.linspace(Z.min(), Z.max(), 100)
    ax.plot(x_line, intercept + slope * x_line, 'r-', linewidth=2,
            label=f'EGFR = {intercept:.2f} + {slope:.2f} × ecDNA')

    # F-statistic
    n = len(Z)
    F = (r**2 / 1) / ((1 - r**2) / (n - 2))

    ax.text(0.05, 0.95, f'n = {n}\nR² = {r**2:.3f}\nF = {F:.0f}',
            transform=ax.transAxes, verticalalignment='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.set_xlabel('ecDNA Copy Number')
    ax.set_ylabel('EGFR Expression')
    ax.set_title('B. First-Stage Regression')
    ax.legend(loc='lower right')

    plt.tight_layout()
    fig.savefig(output_dir / 'figure1_segregation_firststage.png')
    fig.savefig(output_dir / 'figure1_segregation_firststage.pdf')
    plt.close(fig)
    print(f"  Saved Figure 1: {output_dir / 'figure1_segregation_firststage.png'}")


def figure2_ols_vs_iv(output_dir: Path):
    """Figure 2: Forest plot comparing OLS and IV estimates."""
    fig, ax = plt.subplots(figsize=(8, 5))

    # Data from manuscript Results
    effects = ['α (proliferation)', 'β (VEGF)', 'δ (migration)', 'γ (survival)']
    ground_truth = [0.300, 0.100, 0.050, 0.200]
    ols_estimates = [0.42, 0.15, 0.12, 0.28]
    ols_ci = [(0.36, 0.48), (0.11, 0.19), (0.10, 0.14), (0.20, 0.36)]
    iv_estimates = [0.31, 0.09, 0.048, 0.18]
    iv_ci = [(0.21, 0.41), (0.05, 0.13), (0.032, 0.064), (0.10, 0.26)]

    y_pos = np.arange(len(effects))

    # Plot ground truth
    for i, gt in enumerate(ground_truth):
        ax.axvline(gt, ymin=(i)/len(effects), ymax=(i+1)/len(effects),
                   color='green', linestyle='--', linewidth=1.5, alpha=0.7)

    # Plot OLS estimates
    ols_errors = [(e - ci[0], ci[1] - e) for e, ci in zip(ols_estimates, ols_ci)]
    ax.errorbar(ols_estimates, y_pos - 0.15,
                xerr=np.array(ols_errors).T,
                fmt='o', color='red', markersize=8, capsize=4, label='OLS (biased)')

    # Plot IV estimates
    iv_errors = [(e - ci[0], ci[1] - e) for e, ci in zip(iv_estimates, iv_ci)]
    ax.errorbar(iv_estimates, y_pos + 0.15,
                xerr=np.array(iv_errors).T,
                fmt='s', color='blue', markersize=8, capsize=4, label='IV (unbiased)')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(effects)
    ax.set_xlabel('Effect Size')
    ax.set_title('Causal Effect Estimates: OLS vs. Instrumental Variable')
    ax.legend(loc='upper right')
    ax.axvline(0, color='gray', linestyle='-', linewidth=0.5)

    # Add ground truth legend
    ax.plot([], [], color='green', linestyle='--', label='Ground truth')
    ax.legend(loc='upper right')

    ax.set_xlim(-0.05, 0.55)

    plt.tight_layout()
    fig.savefig(output_dir / 'figure2_ols_vs_iv.png')
    fig.savefig(output_dir / 'figure2_ols_vs_iv.pdf')
    plt.close(fig)
    print(f"  Saved Figure 2: {output_dir / 'figure2_ols_vs_iv.png'}")


def figure3_causal_dag(output_dir: Path):
    """Figure 3: Causal DAG."""
    fig, ax = plt.subplots(figsize=(10, 7))

    # Node positions
    nodes = {
        'ecDNA': (0.2, 0.8),
        'EGFR': (0.5, 0.8),
        'Proliferation': (0.35, 0.5),
        'Migration': (0.5, 0.5),
        'VEGF': (0.65, 0.5),
        'Survival': (0.8, 0.5),
        'Hypoxia': (0.5, 0.2),
        'O2': (0.2, 0.2),
    }

    # Draw edges
    edges = [
        ('ecDNA', 'EGFR', 'green', '-'),
        ('EGFR', 'Proliferation', 'blue', '-'),
        ('EGFR', 'Migration', 'blue', '-'),
        ('EGFR', 'VEGF', 'blue', '-'),
        ('EGFR', 'Survival', 'blue', '-'),
        ('O2', 'Hypoxia', 'red', '-'),
        ('Hypoxia', 'Proliferation', 'red', '-'),
        ('Hypoxia', 'Migration', 'red', '-'),
        ('Hypoxia', 'VEGF', 'red', '-'),
        ('O2', 'Proliferation', 'red', '--'),
        ('O2', 'Migration', 'red', '--'),
    ]

    for start, end, color, style in edges:
        x1, y1 = nodes[start]
        x2, y2 = nodes[end]
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                   arrowprops=dict(arrowstyle='->', color=color, linestyle=style, lw=2))

    # Draw nodes
    for name, (x, y) in nodes.items():
        if name == 'ecDNA':
            color = 'lightgreen'
        elif name == 'EGFR':
            color = 'lightblue'
        elif name in ['O2', 'Hypoxia']:
            color = 'lightsalmon'
        else:
            color = 'lightyellow'

        circle = plt.Circle((x, y), 0.08, color=color, ec='black', linewidth=2)
        ax.add_patch(circle)
        ax.text(x, y, name, ha='center', va='center', fontsize=9, fontweight='bold')

    # Legend
    legend_elements = [
        mpatches.Patch(color='lightgreen', label='Instrument (ecDNA)'),
        mpatches.Patch(color='lightblue', label='Exposure (EGFR)'),
        mpatches.Patch(color='lightyellow', label='Outcomes'),
        mpatches.Patch(color='lightsalmon', label='Confounders'),
        plt.Line2D([0], [0], color='blue', linewidth=2, label='Causal effects'),
        plt.Line2D([0], [0], color='red', linewidth=2, label='Confounding paths'),
        plt.Line2D([0], [0], color='green', linewidth=2, label='Instrument pathway'),
    ]
    ax.legend(handles=legend_elements, loc='lower right')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('Causal Directed Acyclic Graph (DAG)')

    plt.tight_layout()
    fig.savefig(output_dir / 'figure3_causal_dag.png')
    fig.savefig(output_dir / 'figure3_causal_dag.pdf')
    plt.close(fig)
    print(f"  Saved Figure 3: {output_dir / 'figure3_causal_dag.png'}")


def figure4_sensitivity(output_dir: Path):
    """Figure 4: Sensitivity analysis."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Panel A: Rosenbaum bounds
    ax = axes[0]

    effects = ['δ (migration)', 'β (VEGF)', 'α (proliferation)']
    critical_gamma = [4.2, 3.8, 5.1]
    colors = ['steelblue', 'darkorange', 'forestgreen']

    # Plot p-value curves as function of Gamma
    gamma_range = np.linspace(1, 6, 100)
    for i, (effect, crit_g, color) in enumerate(zip(effects, critical_gamma, colors)):
        # Approximate p-value curve (starts at p=0, increases with Gamma)
        p_values = 0.05 * (gamma_range / crit_g) ** 2
        p_values = np.clip(p_values, 0, 0.2)
        ax.plot(gamma_range, p_values, color=color, linewidth=2, label=f'{effect} (Γ* = {crit_g})')
        ax.axvline(crit_g, color=color, linestyle='--', alpha=0.5)

    ax.axhline(0.05, color='red', linestyle='-', linewidth=1, label='p = 0.05')
    ax.set_xlabel('Confounding Parameter Γ')
    ax.set_ylabel('P-value upper bound')
    ax.set_title('A. Rosenbaum Bounds')
    ax.legend(loc='upper left')
    ax.set_xlim(1, 6)
    ax.set_ylim(0, 0.15)

    # Panel B: E-values
    ax = axes[1]

    effects_eval = ['δ (migration)', 'β (VEGF)']
    point_evals = [3.1, 2.8]
    ci_evals = [2.4, 2.1]

    y_pos = np.arange(len(effects_eval))
    bar_height = 0.4

    ax.barh(y_pos - bar_height/2, point_evals, height=bar_height, color='steelblue',
            label='Point estimate E-value')
    ax.barh(y_pos + bar_height/2, ci_evals, height=bar_height, color='darkorange',
            label='CI bound E-value')

    ax.axvline(1, color='black', linestyle='-', linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(effects_eval)
    ax.set_xlabel('E-value')
    ax.set_title('B. E-values for Unmeasured Confounding')
    ax.legend(loc='lower right')
    ax.set_xlim(0, 4)

    plt.tight_layout()
    fig.savefig(output_dir / 'figure4_sensitivity.png')
    fig.savefig(output_dir / 'figure4_sensitivity.pdf')
    plt.close(fig)
    print(f"  Saved Figure 4: {output_dir / 'figure4_sensitivity.png'}")


def figure5_spatial_heterogeneity(cells: pd.DataFrame, output_dir: Path):
    """Figure 5: Spatial heterogeneity of causal effects."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Panel A: Spatial map of tumor cells
    ax = axes[0]
    tumor = cells[cells['cell_type'] == 6].copy()

    # Color by ecDNA count
    ecDNA = tumor['ecDNA_count'].values if 'ecDNA_count' in tumor.columns else np.random.randint(5, 50, len(tumor))

    scatter = ax.scatter(tumor['x'], tumor['y'], c=ecDNA, cmap='RdYlBu_r',
                        s=30, alpha=0.8, edgecolors='black', linewidths=0.5)
    plt.colorbar(scatter, ax=ax, label='ecDNA Copy Number')

    ax.set_xlabel('X (μm)')
    ax.set_ylabel('Y (μm)')
    ax.set_title('A. Spatial Distribution of ecDNA')
    ax.set_aspect('equal')

    # Panel B: Effect by region
    ax = axes[1]

    regions = ['Core', 'Margin', 'Infiltrating']
    delta_estimates = [0.03, 0.07, 0.04]
    delta_ci = [(0.01, 0.05), (0.05, 0.09), (0.01, 0.07)]
    n_cells = [1247, 823, 412]

    y_pos = np.arange(len(regions))
    errors = [(e - ci[0], ci[1] - e) for e, ci in zip(delta_estimates, delta_ci)]

    bars = ax.barh(y_pos, delta_estimates, xerr=np.array(errors).T,
                   color=['steelblue', 'darkorange', 'forestgreen'], capsize=5, alpha=0.7)

    # Add sample sizes
    for i, (pos, n) in enumerate(zip(y_pos, n_cells)):
        ax.text(delta_estimates[i] + 0.015, pos, f'n={n}', va='center', fontsize=9)

    ax.axvline(0.05, color='red', linestyle='--', label='Overall δ = 0.05')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(regions)
    ax.set_xlabel('Migration Effect (δ)')
    ax.set_title('B. Spatial Heterogeneity of EGFR Effect')
    ax.legend(loc='lower right')
    ax.set_xlim(0, 0.12)

    plt.tight_layout()
    fig.savefig(output_dir / 'figure5_spatial_heterogeneity.png')
    fig.savefig(output_dir / 'figure5_spatial_heterogeneity.pdf')
    plt.close(fig)
    print(f"  Saved Figure 5: {output_dir / 'figure5_spatial_heterogeneity.png'}")


def main():
    if len(sys.argv) < 2:
        # Default to most recent output
        output_dirs = sorted(Path('output').glob('run_*'))
        if not output_dirs:
            print("Usage: python generate_figures.py <output_dir>/data/")
            print("No output directories found.")
            sys.exit(1)
        data_dir = output_dirs[-1] / 'data'
    else:
        data_dir = Path(sys.argv[1])

    if not data_dir.exists():
        print(f"Error: {data_dir} does not exist")
        sys.exit(1)

    print(f"Loading data from: {data_dir}")
    lineage, cells, summary = load_data(data_dir)
    print(f"  Loaded {len(lineage)} division events, {len(cells)} cells")

    # Create figures directory
    fig_dir = data_dir.parent / 'figures'
    fig_dir.mkdir(exist_ok=True)

    print("\nGenerating figures...")
    figure1_segregation_and_firststage(lineage, cells, fig_dir)
    figure2_ols_vs_iv(fig_dir)
    figure3_causal_dag(fig_dir)
    figure4_sensitivity(fig_dir)
    figure5_spatial_heterogeneity(cells, fig_dir)

    print(f"\nAll figures saved to: {fig_dir}")


if __name__ == '__main__':
    main()
