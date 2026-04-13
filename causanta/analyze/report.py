"""Analysis report generation for CAUSANTA.

Generates comprehensive HTML reports with:
- Interactive DAG editor
- All causal analysis results
- Visualizations and comparisons to ground truth
"""

from __future__ import annotations

import base64
import io
import json
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from .loader import SimulationData
from .iv import estimate_iv_effects, estimate_segregation_iv
from .discovery import discover_causal_structure, compare_to_ground_truth
from .matching import propensity_score_matching, inverse_propensity_weighting
from .regression import linear_regression_adjustment, doubly_robust_estimation
from .sem import mediation_analysis, full_sem
from .sensitivity import rosenbaum_bounds, omitted_variable_bias, placebo_test

# Try to import matplotlib
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def _fig_to_base64(fig) -> str:
    """Convert matplotlib figure to base64 PNG."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return img_base64


def _create_effect_comparison_chart(
    methods: list[str],
    estimates: list[float],
    errors: list[float],
    ground_truth: float | None = None,
) -> str:
    """Create chart comparing effect estimates across methods."""
    if not HAS_MATPLOTLIB:
        return ""

    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.arange(len(methods))
    colors = ['#3498db', '#e74c3c', '#27ae60', '#9b59b6', '#f39c12', '#1abc9c']

    bars = ax.bar(x, estimates, yerr=errors, capsize=5,
                  color=colors[:len(methods)], alpha=0.8, edgecolor='black')

    if ground_truth is not None:
        ax.axhline(ground_truth, color='red', linestyle='--',
                   linewidth=2, label=f'Ground Truth: {ground_truth:.3f}')
        ax.legend()

    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=45, ha='right')
    ax.set_ylabel('Effect Estimate')
    ax.set_title('Causal Effect Estimates Across Methods')
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    return _fig_to_base64(fig)


def _create_dag_comparison_chart(
    precision: float,
    recall: float,
    f1: float,
) -> str:
    """Create bar chart for DAG discovery metrics."""
    if not HAS_MATPLOTLIB:
        return ""

    fig, ax = plt.subplots(figsize=(8, 5))

    metrics = ['Precision', 'Recall', 'F1 Score']
    values = [precision, recall, f1]
    colors = ['#3498db', '#27ae60', '#9b59b6']

    bars = ax.bar(metrics, values, color=colors, alpha=0.8, edgecolor='black')

    ax.set_ylim(0, 1)
    ax.set_ylabel('Score')
    ax.set_title('Causal Discovery Performance')
    ax.grid(True, alpha=0.3, axis='y')

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{val:.2f}', ha='center', fontsize=12)

    plt.tight_layout()
    return _fig_to_base64(fig)


def generate_analysis_report(
    data: SimulationData,
    output_path: Path,
    methods: list[str] | None = None,
    dag_config: dict | None = None,
    include_json: bool = False,
) -> Path:
    """Generate comprehensive analysis report.

    Args:
        data: Loaded simulation data
        output_path: Where to save the report
        methods: Which analysis methods to run
        dag_config: Optional DAG configuration
        include_json: Whether to also output JSON results

    Returns:
        Path to generated report.
    """
    if methods is None:
        methods = ["iv", "discovery", "matching", "regression", "sem", "sensitivity"]

    results = {}
    charts = {}

    # Run each analysis method
    print("  Running IV estimation...")
    if "iv" in methods:
        iv_results = estimate_iv_effects(data)
        segregation_iv = estimate_segregation_iv(data)
        results["iv"] = {
            "effects": [r.to_dict() for r in iv_results],
            "segregation": segregation_iv.to_dict() if segregation_iv else None,
        }

    print("  Running causal discovery...")
    if "discovery" in methods:
        pc_dag = discover_causal_structure(data, algorithm="pc")
        ges_dag = discover_causal_structure(data, algorithm="ges")

        # Get ground truth edges for comparison
        ground_truth_edges = [
            ("ecDNA_count", "VEGF_secretion"),
            ("ecDNA_count", "migration_rate"),
            ("O2_local", "VEGF_secretion"),
        ]

        pc_comparison = compare_to_ground_truth(pc_dag, ground_truth_edges)
        ges_comparison = compare_to_ground_truth(ges_dag, ground_truth_edges)

        results["discovery"] = {
            "pc": {
                "dag": pc_dag.to_dict(),
                "comparison": pc_comparison,
            },
            "ges": {
                "dag": ges_dag.to_dict(),
                "comparison": ges_comparison,
            },
        }

        # Create comparison chart
        if HAS_MATPLOTLIB:
            charts["discovery"] = _create_dag_comparison_chart(
                pc_comparison["precision"],
                pc_comparison["recall"],
                pc_comparison["f1_score"],
            )

    print("  Running propensity score matching...")
    if "matching" in methods:
        matching_result = propensity_score_matching(data)
        ipw_result = inverse_propensity_weighting(data)
        results["matching"] = {
            "psm": matching_result.to_dict(),
            "ipw": ipw_result,
        }

    print("  Running regression analysis...")
    if "regression" in methods:
        reg_vegf = linear_regression_adjustment(data, outcome_name="VEGF_secretion")
        reg_migration = linear_regression_adjustment(data, outcome_name="migration_rate")
        dr_result = doubly_robust_estimation(data)
        results["regression"] = {
            "vegf": reg_vegf.to_dict(),
            "migration": reg_migration.to_dict(),
            "doubly_robust": dr_result,
        }

    print("  Running SEM analysis...")
    if "sem" in methods:
        mediation = mediation_analysis(data)
        full = full_sem(data)
        results["sem"] = {
            "mediation": mediation.to_dict(),
            "full_model": full.to_dict(),
        }

    print("  Running sensitivity analysis...")
    if "sensitivity" in methods:
        rosenbaum = rosenbaum_bounds(data)
        ovb = omitted_variable_bias(data)
        placebo = placebo_test(data)
        results["sensitivity"] = {
            "rosenbaum": rosenbaum.to_dict(),
            "omitted_variable": ovb.to_dict(),
            "placebo": placebo,
        }

    # Create effect comparison chart
    if HAS_MATPLOTLIB:
        effect_methods = []
        effect_estimates = []
        effect_errors = []

        if "regression" in results:
            effect_methods.append("OLS")
            effect_estimates.append(results["regression"]["vegf"]["treatment_coef"])
            effect_errors.append(results["regression"]["vegf"]["treatment_se"])

        if "matching" in results:
            effect_methods.append("PSM")
            effect_estimates.append(results["matching"]["psm"]["att"])
            effect_errors.append(results["matching"]["psm"]["att_se"])

            effect_methods.append("IPW")
            effect_estimates.append(results["matching"]["ipw"].get("ate", 0))
            effect_errors.append(results["matching"]["ipw"].get("ate_se", 0))

        if "regression" in results and "doubly_robust" in results["regression"]:
            effect_methods.append("DR")
            effect_estimates.append(results["regression"]["doubly_robust"].get("ate", 0))
            effect_errors.append(results["regression"]["doubly_robust"].get("ate_se", 0))

        if effect_methods:
            # Get ground truth from config
            ground_truth = None
            if data.config:
                tumor_cfg = data.config.get("cell_types", {}).get("6", {})
                ground_truth = tumor_cfg.get("ecDNA_effect_on_VEGF", None)

            charts["effects"] = _create_effect_comparison_chart(
                effect_methods, effect_estimates, effect_errors, ground_truth
            )

    # Generate HTML
    html = _build_analysis_html(data, results, charts, dag_config)

    output_path = Path(output_path)
    with open(output_path, "w") as f:
        f.write(html)

    # Optionally write JSON
    if include_json:
        json_path = output_path.with_suffix(".json")
        with open(json_path, "w") as f:
            json.dump(results, f, indent=2, default=str)

    return output_path


def _build_analysis_html(
    data: SimulationData,
    results: dict[str, Any],
    charts: dict[str, str],
    dag_config: dict | None,
) -> str:
    """Build complete HTML report."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Build DAG editor (Mermaid + interactive controls)
    dag_mermaid = _build_dag_mermaid(dag_config)

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CAUSANTA Analysis Report</title>
    <script src="https://cdn.jsdelivr.net/npm/mermaid/dist/mermaid.min.js"></script>
    <style>
        :root {{
            --bg-dark: #1a1a2e;
            --bg-card: #16213e;
            --text-primary: #eee;
            --text-secondary: #aaa;
            --accent: #e94560;
            --accent2: #0f3460;
            --success: #27ae60;
            --warning: #f39c12;
        }}
        * {{ box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: var(--bg-dark);
            color: var(--text-primary);
            margin: 0;
            padding: 20px;
            line-height: 1.6;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        h1 {{ color: var(--accent); border-bottom: 2px solid var(--accent); padding-bottom: 10px; }}
        h2 {{ color: var(--accent); margin-top: 40px; border-left: 4px solid var(--accent); padding-left: 15px; }}
        h3 {{ color: var(--text-primary); margin-top: 25px; }}
        .card {{
            background: var(--bg-card);
            border-radius: 8px;
            padding: 20px;
            margin: 20px 0;
            box-shadow: 0 4px 6px rgba(0,0,0,0.3);
        }}
        .grid {{ display: grid; gap: 20px; }}
        .grid-2 {{ grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
        }}
        th, td {{
            padding: 10px;
            text-align: left;
            border-bottom: 1px solid #333;
        }}
        th {{ background: var(--accent2); }}
        .mermaid {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            margin: 15px 0;
        }}
        .dag-editor {{
            background: var(--bg-card);
            padding: 20px;
            border-radius: 8px;
            margin: 20px 0;
        }}
        .dag-controls {{
            display: flex;
            gap: 10px;
            margin-bottom: 15px;
            flex-wrap: wrap;
        }}
        .dag-controls button {{
            padding: 8px 16px;
            border: none;
            border-radius: 4px;
            background: var(--accent);
            color: white;
            cursor: pointer;
        }}
        .dag-controls button:hover {{ opacity: 0.8; }}
        .dag-controls input, .dag-controls select {{
            padding: 8px;
            border-radius: 4px;
            border: 1px solid #444;
            background: #2d3436;
            color: white;
        }}
        .metric-good {{ color: var(--success); }}
        .metric-bad {{ color: var(--accent); }}
        .img-container {{ text-align: center; margin: 15px 0; }}
        .img-container img {{ max-width: 100%; border-radius: 4px; }}
        code {{ background: #2d3436; padding: 2px 6px; border-radius: 3px; }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #333;
            text-align: center;
            color: var(--text-secondary);
        }}
        .stat-box {{
            background: var(--accent2);
            padding: 15px;
            border-radius: 8px;
            text-align: center;
        }}
        .stat-value {{ font-size: 2em; font-weight: bold; color: var(--accent); }}
        .stat-label {{ color: var(--text-secondary); font-size: 0.9em; }}
        .collapsible {{
            cursor: pointer;
            padding: 10px;
            background: var(--accent2);
            border: none;
            text-align: left;
            width: 100%;
            color: var(--text-primary);
            border-radius: 4px;
            margin-top: 10px;
        }}
        .collapsible:after {{ content: ' +'; float: right; }}
        .collapsible.active:after {{ content: ' -'; }}
        .content {{
            padding: 0 18px;
            display: none;
            overflow: hidden;
            background: var(--bg-card);
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>CAUSANTA Causal Analysis Report</h1>

        <!-- Data Summary -->
        <div class="card">
            <h3>Data Summary</h3>
            <div class="grid grid-2">
                <div class="stat-box">
                    <div class="stat-value">{len(data.tumor_cells)}</div>
                    <div class="stat-label">Tumor Cells</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{data.n_divisions}</div>
                    <div class="stat-label">Division Events</div>
                </div>
            </div>
            <table>
                <tr><th>Parameter</th><th>Value</th></tr>
                <tr><td>Analysis Date</td><td>{timestamp}</td></tr>
                <tr><td>Simulation Directory</td><td><code>{data.output_dir}</code></td></tr>
                <tr><td>Timesteps</td><td>{data.n_timesteps}</td></tr>
            </table>
        </div>

        <!-- DAG Editor Section -->
        <h2>Causal DAG (Interactive Editor)</h2>
        <div class="dag-editor">
            <div class="dag-controls">
                <select id="dag-preset">
                    <option value="default">Default CAUSANTA Model</option>
                    <option value="simple">Simple Model</option>
                    <option value="full">Full Model</option>
                </select>
                <button onclick="updateDAG()">Update DAG</button>
                <button onclick="exportDAG()">Export JSON</button>
            </div>
            <div class="mermaid" id="dag-display">
{dag_mermaid}
            </div>
            <p><em>The DAG above shows the assumed causal structure. Modify the preset or export to customize.</em></p>
        </div>

        {_build_iv_section(results.get("iv", {}))}
        {_build_discovery_section(results.get("discovery", {}), charts.get("discovery", ""))}
        {_build_matching_section(results.get("matching", {}))}
        {_build_regression_section(results.get("regression", {}))}
        {_build_sem_section(results.get("sem", {}))}
        {_build_sensitivity_section(results.get("sensitivity", {}))}

        <!-- Effect Comparison -->
        {_build_effects_comparison(charts.get("effects", ""), results)}

        <div class="footer">
            <p>Generated by CAUSANTA v0.2.0</p>
            <p>Causal Analysis Using Somatic And Neighborhood Tissue Architecture</p>
        </div>
    </div>

    <script>
        mermaid.initialize({{ startOnLoad: true, theme: 'default' }});

        function updateDAG() {{
            const preset = document.getElementById('dag-preset').value;
            let mermaidCode = '';

            if (preset === 'simple') {{
                mermaidCode = `graph TD
                    ecDNA[ecDNA Count] --> VEGF[VEGF Secretion]
                    ecDNA --> Migration[Migration Rate]
                    style ecDNA fill:#e1f5fe,stroke:#01579b
                    style VEGF fill:#e8f5e9,stroke:#1b5e20
                    style Migration fill:#e8f5e9,stroke:#1b5e20`;
            }} else if (preset === 'full') {{
                mermaidCode = `graph TD
                    ecDNA[ecDNA Count] --> EGFR[EGFR Expression]
                    ecDNA --> MYC[MYC Expression]
                    EGFR --> Prolif[Proliferation]
                    MYC --> Prolif
                    ecDNA --> VEGF[VEGF Secretion]
                    VEGF --> Angio[Angiogenesis]
                    Angio --> O2[O2 Local]
                    O2 --> Prolif
                    ecDNA --> Migration[Migration]
                    ecDNA --> Survival[Survival]
                    style ecDNA fill:#e1f5fe,stroke:#01579b
                    style O2 fill:#ffebee,stroke:#b71c1c`;
            }} else {{
                mermaidCode = `{dag_mermaid.strip()}`;
            }}

            document.getElementById('dag-display').innerHTML = mermaidCode;
            mermaid.init(undefined, document.getElementById('dag-display'));
        }}

        function exportDAG() {{
            const dagData = {{
                preset: document.getElementById('dag-preset').value,
                timestamp: new Date().toISOString()
            }};
            const blob = new Blob([JSON.stringify(dagData, null, 2)], {{type: 'application/json'}});
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'dag_config.json';
            a.click();
        }}

        // Collapsible sections
        document.querySelectorAll('.collapsible').forEach(btn => {{
            btn.addEventListener('click', function() {{
                this.classList.toggle('active');
                const content = this.nextElementSibling;
                content.style.display = content.style.display === 'block' ? 'none' : 'block';
            }});
        }});
    </script>
</body>
</html>'''

    return html


def _build_dag_mermaid(dag_config: dict | None) -> str:
    """Build Mermaid DAG diagram."""
    return '''graph TD
    ecDNA[ecDNA Count] -->|β| VEGF[VEGF Secretion]
    ecDNA -->|δ| Migration[Migration Rate]
    ecDNA -->|α| Prolif[Proliferation Rate]
    ecDNA -->|γ| Survival[Survival]
    VEGF --> Angio[Angiogenesis]
    Angio --> O2[O2 Local]
    O2 -->|confounder| Prolif
    O2 -->|confounder| VEGF

    style ecDNA fill:#e1f5fe,stroke:#01579b
    style VEGF fill:#fff3e0,stroke:#e65100
    style Migration fill:#e8f5e9,stroke:#1b5e20
    style Prolif fill:#e8f5e9,stroke:#1b5e20
    style Survival fill:#e8f5e9,stroke:#1b5e20
    style O2 fill:#ffebee,stroke:#b71c1c'''


def _build_iv_section(iv_results: dict) -> str:
    """Build IV estimation section."""
    if not iv_results:
        return ""

    effects = iv_results.get("effects", [])
    seg = iv_results.get("segregation")

    rows = ""
    for e in effects:
        weak = "Yes" if e.get("weak_instrument") else "No"
        rows += f'''<tr>
            <td>{e.get("outcome", "")}</td>
            <td>{e.get("second_stage_coef", 0):.4f}</td>
            <td>{e.get("second_stage_se", 0):.4f}</td>
            <td>{e.get("ols_coef", 0):.4f}</td>
            <td>{e.get("first_stage_f_stat", 0):.1f}</td>
            <td class="{'metric-bad' if weak == 'Yes' else 'metric-good'}">{weak}</td>
        </tr>'''

    seg_html = ""
    if seg:
        seg_html = f'''
        <h4>Segregation-Based IV</h4>
        <p>Using ecDNA segregation randomness as instrument:</p>
        <ul>
            <li>Effect estimate: {seg.get("second_stage_coef", 0):.4f}</li>
            <li>F-statistic: {seg.get("first_stage_f_stat", 0):.1f}</li>
        </ul>'''

    return f'''
    <h2>Instrumental Variable Estimation</h2>
    <div class="card">
        <p>Using ecDNA copy number as an instrument for gene expression effects.</p>
        <table>
            <tr>
                <th>Outcome</th>
                <th>IV Estimate</th>
                <th>SE</th>
                <th>OLS Estimate</th>
                <th>First-Stage F</th>
                <th>Weak Instrument</th>
            </tr>
            {rows}
        </table>
        {seg_html}
    </div>'''


def _build_discovery_section(discovery_results: dict, chart: str) -> str:
    """Build causal discovery section."""
    if not discovery_results:
        return ""

    pc = discovery_results.get("pc", {})
    ges = discovery_results.get("ges", {})

    pc_comp = pc.get("comparison", {})
    ges_comp = ges.get("comparison", {})

    chart_html = f'<div class="img-container"><img src="data:image/png;base64,{chart}" alt="Discovery Metrics"></div>' if chart else ""

    return f'''
    <h2>Causal Discovery</h2>
    <div class="card">
        <div class="grid grid-2">
            <div>
                <h4>PC Algorithm</h4>
                <ul>
                    <li>Edges found: {len(pc.get("dag", {}).get("edges", []))}</li>
                    <li>Precision: {pc_comp.get("precision", 0):.2f}</li>
                    <li>Recall: {pc_comp.get("recall", 0):.2f}</li>
                    <li>F1: {pc_comp.get("f1_score", 0):.2f}</li>
                    <li>SHD: {pc_comp.get("structural_hamming_distance", 0)}</li>
                </ul>
            </div>
            <div>
                <h4>GES Algorithm</h4>
                <ul>
                    <li>Edges found: {len(ges.get("dag", {}).get("edges", []))}</li>
                    <li>Precision: {ges_comp.get("precision", 0):.2f}</li>
                    <li>Recall: {ges_comp.get("recall", 0):.2f}</li>
                    <li>F1: {ges_comp.get("f1_score", 0):.2f}</li>
                    <li>SHD: {ges_comp.get("structural_hamming_distance", 0)}</li>
                </ul>
            </div>
        </div>
        {chart_html}
    </div>'''


def _build_matching_section(matching_results: dict) -> str:
    """Build propensity score matching section."""
    if not matching_results:
        return ""

    psm = matching_results.get("psm", {})
    ipw = matching_results.get("ipw", {})

    return f'''
    <h2>Propensity Score Methods</h2>
    <div class="card">
        <div class="grid grid-2">
            <div>
                <h4>Propensity Score Matching</h4>
                <ul>
                    <li>Matched pairs: {psm.get("n_matched", 0)}</li>
                    <li>ATT: {psm.get("att", 0):.4f} (SE: {psm.get("att_se", 0):.4f})</li>
                    <li>Method: {psm.get("method", "")}</li>
                </ul>
            </div>
            <div>
                <h4>Inverse Propensity Weighting</h4>
                <ul>
                    <li>ATE: {ipw.get("ate", 0):.4f} (SE: {ipw.get("ate_se", 0):.4f})</li>
                    <li>N treated: {ipw.get("n_treated", 0)}</li>
                    <li>N control: {ipw.get("n_control", 0)}</li>
                </ul>
            </div>
        </div>
    </div>'''


def _build_regression_section(regression_results: dict) -> str:
    """Build regression analysis section."""
    if not regression_results:
        return ""

    vegf = regression_results.get("vegf", {})
    dr = regression_results.get("doubly_robust", {})

    return f'''
    <h2>Regression Analysis</h2>
    <div class="card">
        <h4>Linear Regression (VEGF Outcome)</h4>
        <table>
            <tr><th>Metric</th><th>Value</th></tr>
            <tr><td>Treatment Effect</td><td>{vegf.get("treatment_coef", 0):.4f}</td></tr>
            <tr><td>Standard Error</td><td>{vegf.get("treatment_se", 0):.4f}</td></tr>
            <tr><td>P-value</td><td>{vegf.get("treatment_pvalue", 0):.4f}</td></tr>
            <tr><td>R²</td><td>{vegf.get("r_squared", 0):.4f}</td></tr>
        </table>

        <h4>Doubly Robust Estimation</h4>
        <ul>
            <li>ATE: {dr.get("ate", 0):.4f}</li>
            <li>95% CI: [{dr.get("ate_ci_lower", 0):.4f}, {dr.get("ate_ci_upper", 0):.4f}]</li>
        </ul>
    </div>'''


def _build_sem_section(sem_results: dict) -> str:
    """Build SEM analysis section."""
    if not sem_results:
        return ""

    mediation = sem_results.get("mediation", {})

    direct = mediation.get("direct_effects", {})
    indirect = mediation.get("indirect_effects", {})
    total = mediation.get("total_effects", {})

    return f'''
    <h2>Structural Equation Modeling</h2>
    <div class="card">
        <h4>Mediation Analysis</h4>
        <table>
            <tr><th>Effect Type</th><th>Value</th></tr>
            <tr><td>Direct Effect</td><td>{list(direct.values())[0] if direct else 0:.4f}</td></tr>
            <tr><td>Total Indirect</td><td>{indirect.get("total_indirect", 0):.4f}</td></tr>
            <tr><td>Total Effect</td><td>{list(total.values())[0] if total else 0:.4f}</td></tr>
        </table>

        <button class="collapsible">Show Path Coefficients</button>
        <div class="content">
            <table>
                <tr><th>From</th><th>To</th><th>Coefficient</th><th>P-value</th></tr>
                {"".join(f"<tr><td>{p.get('from', '')}</td><td>{p.get('to', '')}</td><td>{p.get('coefficient', 0):.4f}</td><td>{p.get('p', 0):.4f}</td></tr>" for p in mediation.get("paths", []))}
            </table>
        </div>
    </div>'''


def _build_sensitivity_section(sensitivity_results: dict) -> str:
    """Build sensitivity analysis section."""
    if not sensitivity_results:
        return ""

    rosenbaum = sensitivity_results.get("rosenbaum", {})
    ovb = sensitivity_results.get("omitted_variable", {})

    return f'''
    <h2>Sensitivity Analysis</h2>
    <div class="card">
        <h4>Rosenbaum Bounds</h4>
        <p>{rosenbaum.get("interpretation", "")}</p>
        <ul>
            <li>Breakdown point (Gamma): {rosenbaum.get("breakdown_point", 0):.2f}</li>
            <li>Robustness value: {rosenbaum.get("robustness_value", 0):.2f}</li>
        </ul>

        <h4>Omitted Variable Bias</h4>
        <p>{ovb.get("interpretation", "")}</p>
        <ul>
            <li>Original estimate: {ovb.get("original_estimate", 0):.4f}</li>
            <li>Robustness value (RV): {ovb.get("robustness_value", 0):.3f}</li>
        </ul>
    </div>'''


def _build_effects_comparison(chart: str, results: dict) -> str:
    """Build effect comparison section."""
    if not chart:
        return ""

    return f'''
    <h2>Effect Estimates Comparison</h2>
    <div class="card">
        <p>Comparison of causal effect estimates across different methods.</p>
        <div class="img-container">
            <img src="data:image/png;base64,{chart}" alt="Effect Comparison">
        </div>
    </div>'''
