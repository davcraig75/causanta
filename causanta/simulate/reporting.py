"""Comprehensive report generation for CAUSANTA simulations.

Generates self-contained HTML reports with:
- Simulation metadata and parameters
- Causal DAG visualization
- Population dynamics charts
- Spatial snapshots
- Environment heatmaps
- Causal effect estimation results
"""

from __future__ import annotations

import base64
import io
import json
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from ..analyze.effects import (
    CausalEffectAnalysis,
    run_causal_analysis,
)
from ..analyze.loader import (
    load_cells_tsv,
    load_lineage_tsv,
    load_summary_log,
)
from ..graph import CausalDAG, build_causanta_dag
from .config import SimulationConfig

# Try to import matplotlib, but make it optional
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def _fig_to_base64(fig) -> str:
    """Convert matplotlib figure to base64 PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return img_base64


def _create_placeholder_image(text: str) -> str:
    """Create a placeholder image when matplotlib is not available."""
    # Simple SVG placeholder
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="400" height="200">
        <rect width="100%" height="100%" fill="#f0f0f0"/>
        <text x="50%" y="50%" text-anchor="middle" fill="#666" font-size="14">
            {text}
        </text>
    </svg>'''
    return base64.b64encode(svg.encode()).decode('utf-8')


class SimulationReport:
    """Generator for comprehensive simulation reports."""

    # Color scheme matching the simulation
    CELL_COLORS = {
        0: "#4477AA",  # Neuron
        1: "#66CCEE",  # Astrocyte
        2: "#228833",  # Oligodendrocyte
        3: "#EE6677",  # Microglia
        4: "#AA3377",  # Endothelial
        5: "#BBBBBB",  # Pericyte
        6: "#CCBB44",  # Tumor
        7: "#EE8866",  # RecruitedImmune
        8: "#555555",  # Necrotic
    }

    CELL_NAMES = {
        0: "Neuron",
        1: "Astrocyte",
        2: "Oligodendrocyte",
        3: "Microglia",
        4: "Endothelial",
        5: "Pericyte",
        6: "Tumor",
        7: "RecruitedImmune",
        8: "Necrotic",
    }

    def __init__(
        self,
        output_dir: Path,
        config: SimulationConfig,
        runtime_seconds: float = 0.0,
        data_dir: Path | None = None,
    ):
        self.output_dir = Path(output_dir)
        self.config = config
        self.runtime_seconds = runtime_seconds

        # Data directory (separate from report output directory)
        self.data_dir = Path(data_dir) if data_dir else self.output_dir

        # Extract parameters for DAG
        tumor_config = config.cell_types.get(6)
        dag_params = {}
        if tumor_config:
            dag_params = {
                "alpha": tumor_config.ecDNA_effect_on_division,
                "beta": tumor_config.ecDNA_effect_on_VEGF,
                "delta": tumor_config.ecDNA_effect_on_migration,
                "gamma": tumor_config.ecDNA_effect_on_survival,
                "O2_prolif_threshold_mmHg": tumor_config.O2_prolif_threshold_mmHg,
                "glucose_prolif_threshold_mM": tumor_config.glucose_prolif_threshold_mM,
            }
        dag_params["angiogenesis_threshold_nM"] = config.angiogenesis.angiogenesis_threshold_nM
        dag_params["q_O2_transfer_per_hr"] = config.environment.q_O2_transfer_per_hr

        # Build causal DAG
        self.dag = build_causanta_dag(dag_params)

        # Load data from data directory
        self.summary_data = load_summary_log(self.data_dir / "summary.log")
        self.lineage_data = load_lineage_tsv(self.data_dir / "lineage.tsv")

        # Find cell snapshot files in data directory
        self.cells_files = sorted(self.data_dir.glob("cells_t*.tsv"))
        self.env_files = sorted(self.data_dir.glob("environment_t*.tsv"))

    def generate(self) -> Path:
        """Generate complete HTML report.

        Returns:
            Path to the generated report file.
        """
        # Run causal analysis
        tumor_config = self.config.cell_types.get(6)
        config_params = {
            "alpha": tumor_config.ecDNA_effect_on_division if tumor_config else 0.0,
            "beta": tumor_config.ecDNA_effect_on_VEGF if tumor_config else 0.0,
            "delta": tumor_config.ecDNA_effect_on_migration if tumor_config else 0.0,
            "gamma": tumor_config.ecDNA_effect_on_survival if tumor_config else 0.0,
            "base_division_time": tumor_config.division_time_mean_hr if tumor_config else 36.0,
            "base_vegf": tumor_config.VEGF_secretion_amol_hr if tumor_config else 600.0,
            "base_migration": tumor_config.migration_speed_um_hr if tumor_config else 10.0,
        }
        causal_analysis = run_causal_analysis(self.data_dir, config_params)

        # Generate all chart images
        charts = self._generate_charts()

        # Build HTML
        html = self._build_html(charts, causal_analysis)

        # Write report
        report_path = self.output_dir / "report.html"
        with open(report_path, "w") as f:
            f.write(html)

        return report_path

    def _generate_charts(self) -> dict[str, str]:
        """Generate all chart images as base64 strings."""
        charts = {}

        if not HAS_MATPLOTLIB:
            # Return placeholders if matplotlib not available
            charts["population"] = _create_placeholder_image("Install matplotlib for charts")
            charts["tumor_growth"] = _create_placeholder_image("Install matplotlib for charts")
            charts["ecdna_distribution"] = _create_placeholder_image("Install matplotlib for charts")
            charts["environment"] = _create_placeholder_image("Install matplotlib for charts")
            charts["spatial_t0"] = _create_placeholder_image("Install matplotlib for charts")
            charts["spatial_mid"] = _create_placeholder_image("Install matplotlib for charts")
            charts["spatial_final"] = _create_placeholder_image("Install matplotlib for charts")
            charts["causal_effects"] = _create_placeholder_image("Install matplotlib for charts")
            return charts

        # Population over time
        charts["population"] = self._chart_population()

        # Tumor growth curve
        charts["tumor_growth"] = self._chart_tumor_growth()

        # ecDNA distribution evolution
        charts["ecdna_distribution"] = self._chart_ecdna_distribution()

        # Environment metrics over time
        charts["environment"] = self._chart_environment()

        # Spatial snapshots
        if self.cells_files:
            charts["spatial_t0"] = self._chart_spatial_snapshot(self.cells_files[0], "t=0")
            mid_idx = len(self.cells_files) // 2
            charts["spatial_mid"] = self._chart_spatial_snapshot(
                self.cells_files[mid_idx], f"t={mid_idx}"
            )
            charts["spatial_final"] = self._chart_spatial_snapshot(
                self.cells_files[-1], f"t={len(self.cells_files)-1}"
            )
        else:
            placeholder = _create_placeholder_image("No spatial data")
            charts["spatial_t0"] = placeholder
            charts["spatial_mid"] = placeholder
            charts["spatial_final"] = placeholder

        # Environment heatmaps
        if self.env_files:
            charts["heatmap_o2"] = self._chart_environment_heatmap(
                self.env_files[-1], "O2", "O2 (mmHg)", "YlOrRd_r"
            )
            charts["heatmap_vegf"] = self._chart_environment_heatmap(
                self.env_files[-1], "VEGF", "VEGF (nM)", "Purples"
            )
        else:
            charts["heatmap_o2"] = _create_placeholder_image("No environment data")
            charts["heatmap_vegf"] = _create_placeholder_image("No environment data")

        # Causal effect comparison
        charts["causal_effects"] = self._chart_causal_effects()

        return charts

    def _chart_population(self) -> str:
        """Generate population by cell type over time chart."""
        if not self.summary_data:
            return _create_placeholder_image("No summary data")

        fig, ax = plt.subplots(figsize=(10, 6))

        steps = [d.get("step", i) for i, d in enumerate(self.summary_data)]
        total = [d.get("total_cells", 0) for d in self.summary_data]
        tumor = [d.get("tumor_cells", 0) for d in self.summary_data]
        necrotic = [d.get("necrotic_cells", 0) for d in self.summary_data]
        immune = [d.get("immune_cells", 0) for d in self.summary_data]

        ax.fill_between(steps, total, alpha=0.3, label="Total", color="#888888")
        ax.plot(steps, tumor, label="Tumor", color=self.CELL_COLORS[6], linewidth=2)
        ax.plot(steps, necrotic, label="Necrotic", color=self.CELL_COLORS[8], linewidth=2)
        ax.plot(steps, immune, label="Immune", color=self.CELL_COLORS[7], linewidth=2)

        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("Cell Count")
        ax.set_title("Cell Population Dynamics")
        ax.legend(loc="upper left")
        ax.grid(True, alpha=0.3)

        return _fig_to_base64(fig)

    def _chart_tumor_growth(self) -> str:
        """Generate tumor growth curve with log scale option."""
        if not self.summary_data:
            return _create_placeholder_image("No summary data")

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        steps = [d.get("step", i) for i, d in enumerate(self.summary_data)]
        tumor = [d.get("tumor_cells", 0) for d in self.summary_data]

        # Linear scale
        ax1.plot(steps, tumor, color=self.CELL_COLORS[6], linewidth=2)
        ax1.fill_between(steps, tumor, alpha=0.3, color=self.CELL_COLORS[6])
        ax1.set_xlabel("Time (hours)")
        ax1.set_ylabel("Tumor Cell Count")
        ax1.set_title("Tumor Growth (Linear)")
        ax1.grid(True, alpha=0.3)

        # Log scale (add 1 to avoid log(0))
        tumor_log = [max(t, 1) for t in tumor]
        ax2.semilogy(steps, tumor_log, color=self.CELL_COLORS[6], linewidth=2)
        ax2.fill_between(steps, tumor_log, alpha=0.3, color=self.CELL_COLORS[6])
        ax2.set_xlabel("Time (hours)")
        ax2.set_ylabel("Tumor Cell Count (log)")
        ax2.set_title("Tumor Growth (Logarithmic)")
        ax2.grid(True, alpha=0.3, which='both')

        plt.tight_layout()
        return _fig_to_base64(fig)

    def _chart_ecdna_distribution(self) -> str:
        """Generate ecDNA distribution evolution histogram."""
        if not self.cells_files:
            return _create_placeholder_image("No cell data")

        # Sample a few timepoints
        n_files = len(self.cells_files)
        indices = [0, n_files // 4, n_files // 2, 3 * n_files // 4, n_files - 1]
        indices = sorted(set(max(0, min(i, n_files - 1)) for i in indices))

        fig, axes = plt.subplots(1, len(indices), figsize=(15, 4), sharey=True)
        if len(indices) == 1:
            axes = [axes]

        for i, idx in enumerate(indices):
            cells = load_cells_tsv(self.cells_files[idx])
            tumor_ecdna = [c["ecDNA_count"] for c in cells if c["cell_type"] == 6]

            if tumor_ecdna:
                axes[i].hist(tumor_ecdna, bins=20, color=self.CELL_COLORS[6],
                           alpha=0.7, edgecolor='black')
                axes[i].axvline(np.mean(tumor_ecdna), color='red', linestyle='--',
                              label=f'Mean: {np.mean(tumor_ecdna):.1f}')
                axes[i].legend(fontsize=8)

            # Extract timestep from filename
            fname = self.cells_files[idx].stem
            t = fname.replace("cells_t", "")
            axes[i].set_title(f"t={int(t)}h")
            axes[i].set_xlabel("ecDNA Count")
            if i == 0:
                axes[i].set_ylabel("Frequency")

        plt.suptitle("ecDNA Distribution Evolution in Tumor Cells", fontsize=12)
        plt.tight_layout()
        return _fig_to_base64(fig)

    def _chart_environment(self) -> str:
        """Generate environment metrics over time chart."""
        if not self.summary_data:
            return _create_placeholder_image("No summary data")

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))

        steps = [d.get("step", i) for i, d in enumerate(self.summary_data)]

        # O2
        o2_mean = [d.get("mean_O2_mmHg", 0) for d in self.summary_data]
        o2_min = [d.get("min_O2_mmHg", 0) for d in self.summary_data]
        axes[0, 0].plot(steps, o2_mean, label="Mean", color="#e74c3c", linewidth=2)
        axes[0, 0].plot(steps, o2_min, label="Min", color="#c0392b", linewidth=1, linestyle="--")
        axes[0, 0].axhline(self.config.environment.hypoxia_threshold_mmHg,
                          color="gray", linestyle=":", label="Hypoxia threshold")
        axes[0, 0].set_ylabel("O2 (mmHg)")
        axes[0, 0].set_title("Oxygen")
        axes[0, 0].legend(fontsize=8)
        axes[0, 0].grid(True, alpha=0.3)

        # Glucose
        glucose = [d.get("mean_glucose_mM", 0) for d in self.summary_data]
        axes[0, 1].plot(steps, glucose, color="#27ae60", linewidth=2)
        axes[0, 1].set_ylabel("Glucose (mM)")
        axes[0, 1].set_title("Glucose")
        axes[0, 1].grid(True, alpha=0.3)

        # VEGF
        vegf = [d.get("max_VEGF_nM", 0) for d in self.summary_data]
        axes[1, 0].plot(steps, vegf, color="#9b59b6", linewidth=2)
        axes[1, 0].axhline(self.config.angiogenesis.angiogenesis_threshold_nM,
                          color="gray", linestyle=":", label="Angio threshold")
        axes[1, 0].set_xlabel("Time (hours)")
        axes[1, 0].set_ylabel("Max VEGF (nM)")
        axes[1, 0].set_title("VEGF")
        axes[1, 0].legend(fontsize=8)
        axes[1, 0].grid(True, alpha=0.3)

        # Lactate
        lactate = [d.get("mean_lactate_mM", 0) for d in self.summary_data]
        axes[1, 1].plot(steps, lactate, color="#f39c12", linewidth=2)
        axes[1, 1].set_xlabel("Time (hours)")
        axes[1, 1].set_ylabel("Lactate (mM)")
        axes[1, 1].set_title("Lactate")
        axes[1, 1].grid(True, alpha=0.3)

        plt.suptitle("Environment Dynamics", fontsize=12)
        plt.tight_layout()
        return _fig_to_base64(fig)

    def _chart_spatial_snapshot(self, cells_path: Path, title: str) -> str:
        """Generate spatial scatter plot of cells."""
        cells = load_cells_tsv(cells_path)
        if not cells:
            return _create_placeholder_image("No cells")

        fig, ax = plt.subplots(figsize=(8, 8))

        # Plot by cell type (non-tumor first, then tumor on top)
        for cell_type in sorted(self.CELL_COLORS.keys()):
            if cell_type == 6:  # Skip tumor for now
                continue
            type_cells = [c for c in cells if c["cell_type"] == cell_type]
            if type_cells:
                x = [c["x"] for c in type_cells]
                y = [c["y"] for c in type_cells]
                ax.scatter(x, y, c=self.CELL_COLORS[cell_type], s=10, alpha=0.5,
                          label=self.CELL_NAMES[cell_type])

        # Plot tumor cells on top
        tumor_cells = [c for c in cells if c["cell_type"] == 6]
        if tumor_cells:
            x = [c["x"] for c in tumor_cells]
            y = [c["y"] for c in tumor_cells]
            # Size by ecDNA
            sizes = [max(20, c["ecDNA_count"] * 2) for c in tumor_cells]
            ax.scatter(x, y, c=self.CELL_COLORS[6], s=sizes, alpha=0.8,
                      label=f"Tumor (n={len(tumor_cells)})", edgecolors='black', linewidths=0.5)

        ax.set_xlim(0, self.config.domain.width_um)
        ax.set_ylim(0, self.config.domain.height_um)
        ax.set_xlabel("X (um)")
        ax.set_ylabel("Y (um)")
        ax.set_title(f"Spatial Distribution - {title}")
        ax.set_aspect('equal')
        ax.legend(loc='upper right', fontsize=8, markerscale=0.8)

        return _fig_to_base64(fig)

    def _chart_environment_heatmap(
        self,
        env_path: Path,
        field: str,
        label: str,
        cmap: str,
    ) -> str:
        """Generate environment field heatmap."""
        # Read environment data
        import csv
        data = []
        with open(env_path, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                data.append({
                    "x_grid": int(row["x_grid"]),
                    "y_grid": int(row["y_grid"]),
                    field: float(row[field]),
                })

        if not data:
            return _create_placeholder_image("No environment data")

        # Build grid
        nx = max(d["x_grid"] for d in data) + 1
        ny = max(d["y_grid"] for d in data) + 1
        grid = np.zeros((ny, nx))

        for d in data:
            grid[d["y_grid"], d["x_grid"]] = d[field]

        fig, ax = plt.subplots(figsize=(8, 8))

        im = ax.imshow(grid, origin='lower', cmap=cmap,
                      extent=[0, self.config.domain.width_um,
                              0, self.config.domain.height_um])
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label(label)

        ax.set_xlabel("X (um)")
        ax.set_ylabel("Y (um)")
        ax.set_title(f"{label} - Final Timepoint")

        return _fig_to_base64(fig)

    def _chart_causal_effects(self) -> str:
        """Generate causal effect comparison bar chart."""
        tumor_config = self.config.cell_types.get(6)
        if not tumor_config:
            return _create_placeholder_image("No tumor config")

        # Get configured values
        effects = {
            "α (division)": tumor_config.ecDNA_effect_on_division,
            "β (VEGF)": tumor_config.ecDNA_effect_on_VEGF,
            "δ (migration)": tumor_config.ecDNA_effect_on_migration,
            "γ (survival)": tumor_config.ecDNA_effect_on_survival,
        }

        fig, ax = plt.subplots(figsize=(10, 6))

        names = list(effects.keys())
        values = list(effects.values())
        x = np.arange(len(names))

        bars = ax.bar(x, values, color=['#3498db', '#9b59b6', '#e74c3c', '#27ae60'],
                     alpha=0.8, edgecolor='black')

        ax.set_xticks(x)
        ax.set_xticklabels(names)
        ax.set_ylabel("Effect Size")
        ax.set_title("Configured Causal Effect Sizes (Ground Truth)")
        ax.grid(True, alpha=0.3, axis='y')

        # Add value labels
        for bar, val in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                   f'{val:.3f}', ha='center', va='bottom', fontsize=10)

        plt.tight_layout()
        return _fig_to_base64(fig)

    def _build_html(
        self,
        charts: dict[str, str],
        causal_analysis: CausalEffectAnalysis,
    ) -> str:
        """Build complete HTML report with GitHub README-style design."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Mermaid DAG
        mermaid_dag = self.dag.to_mermaid()

        # Format runtime
        if self.runtime_seconds > 0:
            mins, secs = divmod(self.runtime_seconds, 60)
            if mins > 0:
                runtime_str = f"{int(mins)}m {secs:.1f}s"
            else:
                runtime_str = f"{secs:.1f}s"
        else:
            runtime_str = "N/A"

        # Get final stats
        final_stats = self.summary_data[-1] if self.summary_data else {}
        initial_stats = self.summary_data[0] if self.summary_data else {}

        # Build HTML with GitHub-style clean design
        html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CAUSANTA Simulation Report</title>
    <script src="https://cdn.jsdelivr.net/npm/mermaid/dist/mermaid.min.js"></script>
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
            background: #ffffff;
            color: #24292f;
            margin: 0;
            padding: 40px 20px;
            line-height: 1.6;
            font-size: 16px;
        }}
        .container {{
            max-width: 980px;
            margin: 0 auto;
        }}
        h1 {{
            font-size: 2em;
            font-weight: 600;
            border-bottom: 1px solid #d0d7de;
            padding-bottom: 0.3em;
            margin-bottom: 16px;
        }}
        h2 {{
            font-size: 1.5em;
            font-weight: 600;
            border-bottom: 1px solid #d0d7de;
            padding-bottom: 0.3em;
            margin-top: 48px;
            margin-bottom: 16px;
        }}
        h3 {{
            font-size: 1.25em;
            font-weight: 600;
            margin-top: 24px;
            margin-bottom: 16px;
        }}
        h4 {{
            font-size: 1em;
            font-weight: 600;
            margin-top: 24px;
            margin-bottom: 16px;
        }}
        p {{
            margin-bottom: 16px;
            color: #57606a;
        }}
        .lead {{
            font-size: 1.1em;
            color: #57606a;
            margin-bottom: 24px;
        }}
        a {{ color: #0969da; text-decoration: none; }}
        a:hover {{ text-decoration: underline; }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 16px;
            margin: 24px 0;
        }}
        .stat-card {{
            background: #f6f8fa;
            border: 1px solid #d0d7de;
            border-radius: 6px;
            padding: 16px;
            text-align: center;
        }}
        .stat-value {{
            font-size: 2em;
            font-weight: 600;
            color: #24292f;
            line-height: 1.2;
        }}
        .stat-label {{
            font-size: 0.85em;
            color: #57606a;
            margin-top: 4px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 16px 0;
            font-size: 14px;
        }}
        th, td {{
            padding: 8px 16px;
            text-align: left;
            border: 1px solid #d0d7de;
        }}
        th {{
            background: #f6f8fa;
            font-weight: 600;
        }}
        tr:hover td {{
            background: #f6f8fa;
        }}
        .img-container {{
            margin: 24px 0;
            text-align: center;
        }}
        .img-container img {{
            max-width: 100%;
            height: auto;
            border: 1px solid #d0d7de;
            border-radius: 6px;
        }}
        .img-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 24px;
            margin: 24px 0;
        }}
        .img-card {{
            border: 1px solid #d0d7de;
            border-radius: 6px;
            overflow: hidden;
        }}
        .img-card img {{
            width: 100%;
            display: block;
        }}
        .img-card-title {{
            padding: 12px 16px;
            background: #f6f8fa;
            font-weight: 600;
            font-size: 14px;
            border-top: 1px solid #d0d7de;
        }}
        .mermaid {{
            background: #f6f8fa;
            padding: 24px;
            border-radius: 6px;
            border: 1px solid #d0d7de;
            margin: 16px 0;
        }}
        code {{
            background: rgba(175, 184, 193, 0.2);
            padding: 0.2em 0.4em;
            border-radius: 6px;
            font-family: ui-monospace, SFMono-Regular, 'SF Mono', Menlo, Consolas, monospace;
            font-size: 85%;
        }}
        .badge {{
            display: inline-block;
            padding: 2px 8px;
            border-radius: 12px;
            font-size: 12px;
            font-weight: 500;
        }}
        .badge-success {{ background: #dafbe1; color: #1a7f37; }}
        .badge-warning {{ background: #fff8c5; color: #9a6700; }}
        .badge-error {{ background: #ffebe9; color: #cf222e; }}
        .effect-good {{ color: #1a7f37; }}
        .effect-warn {{ color: #9a6700; }}
        .effect-bad {{ color: #cf222e; }}
        .legend-list {{
            list-style: none;
            padding: 0;
            margin: 16px 0;
        }}
        .legend-list li {{
            padding: 4px 0;
            display: flex;
            align-items: center;
            gap: 8px;
        }}
        .legend-dot {{
            width: 12px;
            height: 12px;
            border-radius: 50%;
            display: inline-block;
        }}
        .note {{
            background: #ddf4ff;
            border: 1px solid #54aeff;
            border-radius: 6px;
            padding: 16px;
            margin: 16px 0;
        }}
        .note-title {{
            font-weight: 600;
            margin-bottom: 8px;
            color: #0969da;
        }}
        .footer {{
            margin-top: 48px;
            padding-top: 24px;
            border-top: 1px solid #d0d7de;
            text-align: center;
            color: #57606a;
            font-size: 14px;
        }}
        hr {{
            border: none;
            border-top: 1px solid #d0d7de;
            margin: 32px 0;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>CAUSANTA Simulation Report</h1>
        <p class="lead">
            This report summarizes a tumor growth simulation using the CAUSANTA framework,
            which models extrachromosomal DNA (ecDNA) dynamics and their causal effects on
            tumor progression in a spatially-explicit microenvironment.
        </p>

        <!-- Quick Stats -->
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">{self.config.time.total_hours}h</div>
                <div class="stat-label">Simulated Time</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{final_stats.get('tumor_cells', 0)}</div>
                <div class="stat-label">Final Tumor Cells</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{len(self.lineage_data)}</div>
                <div class="stat-label">Division Events</div>
            </div>
            <div class="stat-card">
                <div class="stat-value">{runtime_str}</div>
                <div class="stat-label">Runtime</div>
            </div>
        </div>

        <table>
            <tr><th>Parameter</th><th>Value</th><th>Description</th></tr>
            <tr><td>Domain Size</td><td>{self.config.domain.width_um} x {self.config.domain.height_um} um</td><td>Physical simulation area</td></tr>
            <tr><td>Grid Resolution</td><td>{self.config.domain.env_grid_um} um</td><td>Environment discretization</td></tr>
            <tr><td>Initial Cells</td><td>{initial_stats.get('total_cells', 0)}</td><td>Total cells at t=0</td></tr>
            <tr><td>RNG Seed</td><td><code>{self.config.rng_seed}</code></td><td>For reproducibility</td></tr>
            <tr><td>Generated</td><td>{timestamp}</td><td>Report timestamp</td></tr>
        </table>

        <!-- Population Dynamics -->
        <h2>Population Dynamics</h2>
        <p>
            The simulation tracks multiple cell populations over time. Normal brain tissue cells
            (neurons, astrocytes, oligodendrocytes, microglia) maintain homeostasis while the
            tumor grows. Cell counts change due to proliferation, apoptosis, necrosis, and
            immune-mediated killing.
        </p>

        <div class="img-grid">
            <div class="img-card">
                <img src="data:image/png;base64,{charts.get('population', '')}" alt="Population">
                <div class="img-card-title">Cell Population Over Time</div>
            </div>
            <div class="img-card">
                <img src="data:image/png;base64,{charts.get('tumor_growth', '')}" alt="Tumor Growth">
                <div class="img-card-title">Tumor Growth (Linear & Log Scale)</div>
            </div>
        </div>

        <p>
            The tumor growth curve reveals the proliferation dynamics. Exponential growth (linear
            on log scale) indicates unrestricted proliferation, while subexponential growth
            suggests nutrient limitation or space constraints. Growth may plateau when the tumor
            outpaces its vascular supply.
        </p>

        <!-- ecDNA Evolution -->
        <h3>ecDNA Distribution Evolution</h3>
        <p>
            Extrachromosomal DNA (ecDNA) segregates randomly during cell division, creating
            heterogeneity in copy number across the tumor population. High ecDNA cells may have
            proliferative advantages due to amplified oncogenes (EGFR, MYC). The histograms below
            show how the distribution evolves over time.
        </p>

        <div class="img-container">
            <img src="data:image/png;base64,{charts.get('ecdna_distribution', '')}" alt="ecDNA Distribution">
        </div>

        <!-- Spatial Snapshots -->
        <h2>Spatial Organization</h2>
        <p>
            The spatial distribution of cells reveals tumor morphology and microenvironmental
            interactions. Tumor cells (yellow) expand from the initial seed, displacing or
            invading normal tissue. The vascular network (magenta endothelial cells) provides
            oxygen and glucose. Immune cells may infiltrate the tumor periphery.
        </p>

        <div class="img-grid">
            <div class="img-card">
                <img src="data:image/png;base64,{charts.get('spatial_t0', '')}" alt="Spatial t=0">
                <div class="img-card-title">Initial State (t=0)</div>
            </div>
            <div class="img-card">
                <img src="data:image/png;base64,{charts.get('spatial_mid', '')}" alt="Spatial mid">
                <div class="img-card-title">Mid Simulation</div>
            </div>
            <div class="img-card">
                <img src="data:image/png;base64,{charts.get('spatial_final', '')}" alt="Spatial final">
                <div class="img-card-title">Final State (t={self.config.time.total_hours}h)</div>
            </div>
        </div>

        <!-- Environment -->
        <h2>Microenvironment</h2>
        <p>
            The tumor microenvironment is characterized by gradients in oxygen, glucose, and
            signaling molecules. Rapidly proliferating tumor cells consume oxygen and glucose
            faster than diffusion can replenish them, creating hypoxic regions. Hypoxic cells
            secrete VEGF to stimulate angiogenesis. Lactate accumulates as a byproduct of
            glycolysis.
        </p>

        <div class="img-container">
            <img src="data:image/png;base64,{charts.get('environment', '')}" alt="Environment">
        </div>

        <h3>Spatial Gradients (Final Timepoint)</h3>
        <p>
            The heatmaps show the spatial distribution of key substrates at the end of the
            simulation. Oxygen depletion (red/yellow = lower O2) near the tumor core can trigger
            necrosis and drive angiogenic signaling. VEGF accumulation (purple) indicates where
            new blood vessels may sprout.
        </p>

        <div class="img-grid">
            <div class="img-card">
                <img src="data:image/png;base64,{charts.get('heatmap_o2', '')}" alt="O2 Heatmap">
                <div class="img-card-title">Oxygen Concentration</div>
            </div>
            <div class="img-card">
                <img src="data:image/png;base64,{charts.get('heatmap_vegf', '')}" alt="VEGF Heatmap">
                <div class="img-card-title">VEGF Concentration</div>
            </div>
        </div>

        <!-- Causal Structure -->
        <h2>Causal Structure</h2>
        <p>
            CAUSANTA embeds a known causal directed acyclic graph (DAG) into the simulation.
            This ground truth structure enables validation of causal inference methods. The
            key insight is that ecDNA copy number acts as a natural instrument because it
            segregates randomly during mitosis, mimicking a randomized experiment.
        </p>

        <p>{self.dag.description}</p>

        <div class="mermaid">
{mermaid_dag}
        </div>

        <h4>Node Types</h4>
        <ul class="legend-list">
            <li><span class="legend-dot" style="background:#2e7d32"></span> <strong>Instrument</strong>: ecDNA copy number (randomized by segregation)</li>
            <li><span class="legend-dot" style="background:#1565c0"></span> <strong>Exposure</strong>: Gene expression levels driven by ecDNA</li>
            <li><span class="legend-dot" style="background:#e65100"></span> <strong>Mediator</strong>: Intermediate cellular states</li>
            <li><span class="legend-dot" style="background:#c2185b"></span> <strong>Outcome</strong>: Observable phenotypes (proliferation, migration)</li>
            <li><span class="legend-dot" style="background:#7b1fa2"></span> <strong>Confounder</strong>: Shared causes (e.g., microenvironment)</li>
        </ul>

        <!-- Causal Effect Estimation -->
        <h2>Causal Effect Estimation</h2>
        <p>
            The table below compares the configured (ground truth) causal effect parameters
            with estimates recovered from the simulation data. Effect parameters control how
            ecDNA copy number influences cell behaviors: &alpha; (division rate), &beta; (VEGF
            secretion), &delta; (migration speed), and &gamma; (survival/apoptosis).
        </p>

        <table>
            <tr>
                <th>Effect</th>
                <th>Parameter</th>
                <th>Configured</th>
                <th>Estimated</th>
                <th>Error</th>
                <th>Method</th>
            </tr>
            {self._format_effect_row(causal_analysis.alpha)}
            {self._format_effect_row(causal_analysis.beta)}
            {self._format_effect_row(causal_analysis.delta)}
            {self._format_effect_row(causal_analysis.gamma)}
        </table>

        <div class="img-container">
            <img src="data:image/png;base64,{charts.get('causal_effects', '')}" alt="Causal Effects">
        </div>

        <h3>ecDNA Segregation Statistics</h3>
        <p>
            During cell division, ecDNA molecules segregate stochastically between daughter
            cells. This table verifies that the segregation mechanism produces the expected
            binomial distribution. Deviations may indicate implementation issues or biological
            constraints (e.g., unequal segregation due to selection).
        </p>

        <table>
            <tr><th>Metric</th><th>Observed</th><th>Expected</th></tr>
            <tr>
                <td>Number of divisions</td>
                <td>{causal_analysis.segregation_stats.get('n_divisions', 0)}</td>
                <td>-</td>
            </tr>
            <tr>
                <td>Mean daughter fraction</td>
                <td>{causal_analysis.segregation_stats.get('mean_daughter_fraction', 0):.3f}</td>
                <td>0.500</td>
            </tr>
            <tr>
                <td>Segregation variance</td>
                <td>{causal_analysis.segregation_stats.get('segregation_variance', 0):.4f}</td>
                <td>{causal_analysis.segregation_stats.get('expected_variance', 0):.4f}</td>
            </tr>
        </table>

        <!-- Configuration Details -->
        <h2>Configuration Details</h2>
        <p>
            The following tables document the key simulation parameters. Tumor cell parameters
            define proliferation, migration, and metabolic rates. Environment parameters control
            substrate dynamics and thresholds for biological responses.
        </p>

        <h3>Tumor Cell Parameters</h3>
        {self._format_tumor_config()}

        <h3>Environment Parameters</h3>
        <table>
            <tr><th>Parameter</th><th>Value</th><th>Description</th></tr>
            <tr><td>O2 blood</td><td>{self.config.environment.O2_blood_mmHg} mmHg</td><td>Arterial oxygen partial pressure</td></tr>
            <tr><td>Hypoxia threshold</td><td>{self.config.environment.hypoxia_threshold_mmHg} mmHg</td><td>O2 level triggering hypoxic response</td></tr>
            <tr><td>Angiogenesis threshold</td><td>{self.config.angiogenesis.angiogenesis_threshold_nM} nM</td><td>VEGF level to initiate sprouting</td></tr>
            <tr><td>Vessel maturation</td><td>{self.config.angiogenesis.vessel_maturation_hr} hr</td><td>Time for new vessels to mature</td></tr>
        </table>

        <hr>

        <div class="footer">
            <p><strong>CAUSANTA</strong> v0.1.0</p>
            <p>Causal Analysis Using Somatic And Neighborhood Tissue Architecture</p>
        </div>
    </div>

    <script>
        mermaid.initialize({{ startOnLoad: true, theme: 'neutral' }});
    </script>
</body>
</html>'''

        return html

    def _format_effect_row(self, effect) -> str:
        """Format a table row for causal effect comparison."""
        rel_err = effect.relative_error
        if rel_err < 0.1:
            err_class = "effect-good"
            badge = '<span class="badge badge-success">Good</span>'
        elif rel_err < 0.5:
            err_class = "effect-warn"
            badge = '<span class="badge badge-warning">Fair</span>'
        else:
            err_class = "effect-bad"
            badge = '<span class="badge badge-error">Poor</span>'

        return f'''<tr>
    <td><strong>{effect.parameter_symbol}</strong></td>
    <td><code>{effect.parameter_name}</code></td>
    <td>{effect.configured_value:.4f}</td>
    <td>{effect.estimated_value:.4f} ± {effect.standard_error:.4f}</td>
    <td class="{err_class}">{rel_err*100:.1f}% {badge}</td>
    <td>{effect.method}</td>
</tr>'''

    def _format_tumor_config(self) -> str:
        """Format tumor cell configuration table."""
        tumor = self.config.cell_types.get(6)
        if not tumor:
            return "<p>No tumor configuration found.</p>"

        return f'''<table>
    <tr><th>Parameter</th><th>Value</th><th>Description</th></tr>
    <tr><td>Division time</td><td>{tumor.division_time_mean_hr} ± {tumor.division_time_std_hr} hr</td><td>Cell cycle duration (mean ± std)</td></tr>
    <tr><td>Migration speed</td><td>{tumor.migration_speed_um_hr} um/hr</td><td>Base migration velocity</td></tr>
    <tr><td>O2 consumption</td><td>{tumor.O2_consumption_amol_hr:,.0f} amol/hr</td><td>Oxygen uptake rate per cell</td></tr>
    <tr><td>Glucose consumption</td><td>{tumor.glucose_consumption_amol_hr:,.0f} amol/hr</td><td>Glucose uptake rate per cell</td></tr>
    <tr><td>VEGF secretion</td><td>{tumor.VEGF_secretion_amol_hr} amol/hr</td><td>Angiogenic factor production</td></tr>
    <tr><td>Initial ecDNA</td><td>{tumor.ecDNA_init_count} copies</td><td>Starting ecDNA copy number</td></tr>
    <tr><td>ecDNA cargo</td><td><code>{', '.join(tumor.ecDNA_cargo) if tumor.ecDNA_cargo else 'None'}</code></td><td>Amplified oncogenes</td></tr>
    <tr><td>&alpha; (division)</td><td>{tumor.ecDNA_effect_on_division}</td><td>ecDNA effect on proliferation rate</td></tr>
    <tr><td>&beta; (VEGF)</td><td>{tumor.ecDNA_effect_on_VEGF}</td><td>ecDNA effect on VEGF secretion</td></tr>
    <tr><td>&delta; (migration)</td><td>{tumor.ecDNA_effect_on_migration}</td><td>ecDNA effect on migration speed</td></tr>
    <tr><td>&gamma; (survival)</td><td>{tumor.ecDNA_effect_on_survival}</td><td>ecDNA effect on apoptosis resistance</td></tr>
</table>'''


def generate_report(
    output_dir: Path,
    config: SimulationConfig,
    runtime_seconds: float = 0.0,
    data_dir: Path | None = None,
) -> Path:
    """Generate simulation report.

    Args:
        output_dir: Path to write report.html
        config: Simulation configuration
        runtime_seconds: How long the simulation took
        data_dir: Path to simulation data directory (if different from output_dir)

    Returns:
        Path to generated report.html
    """
    report = SimulationReport(output_dir, config, runtime_seconds, data_dir=data_dir)
    return report.generate()
