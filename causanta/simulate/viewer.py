"""Self-contained HTML viewer generation for CAUSANTA simulation results.

Generates a standalone HTML file with SVG-based cell visualization
that properly renders cell shapes without overlap.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from .cells import CellPopulation, TYPE_NAMES
from .config import SimulationConfig
from .domain import Domain
from .environment import EnvironmentFields
from .shapes import parse_svg_path_bounds


def cells_to_json(population: CellPopulation, config: SimulationConfig) -> list[dict]:
    """Convert cell population to JSON-serializable list with full shape info."""
    cells = []
    for cell in population.iter_living():
        # Get shape path from config or cell
        shape_path = cell.shape_path
        if not shape_path:
            tp = config.cell_types.get(cell.cell_type)
            if tp:
                shape_path = tp.shape_path

        cells.append({
            "cell_id": cell.cell_id,
            "parent_id": cell.parent_id,
            "cell_type": cell.cell_type,
            "cell_type_name": cell.cell_type_name,
            "x": cell.x,
            "y": cell.y,
            "diameter": round(cell.diameter, 2),
            "angle": round(cell.angle, 4),
            "shape_path": shape_path,
            "ecDNA_count": cell.ecDNA_count,
            "O2_local": round(cell.O2_local, 2),
            "glucose_local": round(cell.glucose_local, 3),
            "is_hypoxic": cell.is_hypoxic,
            "cell_cycle_phase": cell.cell_cycle_phase,
            "generation": cell.generation,
            "size_scale": round(cell.size_scale, 2),
        })
    return cells


def environment_to_json(
    env: EnvironmentFields,
    domain: Domain,
    downsample: int = 1,
) -> list[dict]:
    """Convert environment fields to JSON-serializable list."""
    records = []
    for gj in range(0, domain.ny, downsample):
        for gi in range(0, domain.nx, downsample):
            x_um, y_um = domain.grid_to_um(gi, gj)
            records.append({
                "x": round(x_um, 1),
                "y": round(y_um, 1),
                "O2": round(float(env.O2[gj, gi]), 2),
                "glucose": round(float(env.glucose[gj, gi]), 3),
                "VEGF": round(float(env.VEGF[gj, gi]), 4),
                "lactate": round(float(env.lactate[gj, gi]), 3),
                "vascular": round(float(env.vascular_density[gj, gi]), 2),
            })
    return records


def _generate_cell_svg(
    cells_data: list[dict],
    type_colors: dict[str, str],
    domain_width: float,
    domain_height: float,
    svg_width: int = 800,
    svg_height: int = 800,
) -> str:
    """Generate SVG string for cell visualization with proper shapes.

    Renders each cell's shape path directly (no symbols) for maximum compatibility.
    """
    scale = min(svg_width / domain_width, svg_height / domain_height)

    svg_parts = [
        f'<svg width="{svg_width}" height="{svg_height}" '
        f'viewBox="0 0 {svg_width} {svg_height}" '
        f'xmlns="http://www.w3.org/2000/svg" '
        f'style="background: #1a1a2e; border-radius: 8px; display: block;">'
    ]

    # Add grid lines for reference
    svg_parts.append('<g class="grid-layer" stroke="#333" stroke-width="0.5" opacity="0.3">')
    for i in range(0, int(domain_width) + 1, 100):
        x = i * scale
        svg_parts.append(f'  <line x1="{x}" y1="0" x2="{x}" y2="{svg_height}"/>')
    for i in range(0, int(domain_height) + 1, 100):
        y = i * scale
        svg_parts.append(f'  <line x1="0" y1="{y}" x2="{svg_width}" y2="{y}"/>')
    svg_parts.append("</g>")

    # Render cells sorted by type for consistent layering
    svg_parts.append('<g class="cells-layer">')

    # Sort cells: necrotic first (bottom), then normal cells, tumor last (top)
    def cell_sort_key(c):
        if c["cell_type_name"] == "Necrotic":
            return (0, c["cell_type"])
        elif c["cell_type_name"] == "Tumor":
            return (2, c["cell_type"])
        else:
            return (1, c["cell_type"])

    sorted_cells = sorted(cells_data, key=cell_sort_key)

    for cell in sorted_cells:
        x = cell["x"] * scale
        y = cell["y"] * scale
        # Scale factor: shape is in [-1,1], cell diameter maps to actual size
        cell_scale = cell["diameter"] * scale / 2.0  # radius in pixels
        angle_deg = math.degrees(cell["angle"])
        color = type_colors.get(cell["cell_type_name"], "#888")
        shape_path = cell.get("shape_path", "")

        # Highlight cells with ecDNA
        if cell["ecDNA_count"] > 5:
            svg_parts.append(
                f'  <circle cx="{x}" cy="{y}" r="{cell_scale + 3}" '
                f'fill="none" stroke="#FFD700" stroke-width="2.5"/>'
            )

        if shape_path:
            # Render the path with transform: translate to position, scale, rotate
            # The path is centered at origin in [-1,1] space
            svg_parts.append(
                f'  <path d="{shape_path}" fill="{color}" stroke="#222" stroke-width="{0.1/cell_scale:.4f}" '
                f'transform="translate({x},{y}) rotate({angle_deg}) scale({cell_scale})">'
                f'<title>ID: {cell["cell_id"]}, Type: {cell["cell_type_name"]}, '
                f'ecDNA: {cell["ecDNA_count"]}, O2: {cell["O2_local"]}</title></path>'
            )
        else:
            # Fallback to circle
            svg_parts.append(
                f'  <circle cx="{x}" cy="{y}" r="{cell_scale}" '
                f'fill="{color}" stroke="#222" stroke-width="1">'
                f'<title>ID: {cell["cell_id"]}, Type: {cell["cell_type_name"]}</title></circle>'
            )

    svg_parts.append("</g>")

    # Add scale bar
    scale_bar_um = 100
    scale_bar_px = scale_bar_um * scale
    svg_parts.append(
        f'<g class="scale-bar" transform="translate(20, {svg_height - 30})">'
        f'  <rect width="{scale_bar_px}" height="4" fill="#eee"/>'
        f'  <text x="{scale_bar_px/2}" y="18" text-anchor="middle" fill="#eee" font-size="12">'
        f'{scale_bar_um} um</text>'
        f'</g>'
    )

    svg_parts.append("</svg>")
    return "\n".join(svg_parts)


def generate_standalone_html(
    population: CellPopulation,
    env: EnvironmentFields,
    domain: Domain,
    config: SimulationConfig,
    time_hr: float,
    output_path: Path,
    title: str = "CAUSANTA Simulation Results",
) -> Path:
    """Generate a self-contained HTML file with SVG-based cell visualization.

    Renders cells with their actual shapes from the config, ensuring
    no visual overlap due to proper collision detection.
    """
    # Convert data to JSON with shape info
    cells_data = cells_to_json(population, config)
    env_data = environment_to_json(env, domain, downsample=2)

    # Build color and shape mapping from config
    type_colors = {}
    type_shapes = {}
    for tid, tc in config.cell_types.items():
        type_colors[tc.type_name] = tc.color
        type_shapes[tc.type_name] = tc.shape_path

    # Calculate stats
    cell_counts = {}
    total_ecDNA = 0
    for cell in cells_data:
        tname = cell["cell_type_name"]
        cell_counts[tname] = cell_counts.get(tname, 0) + 1
        total_ecDNA += cell["ecDNA_count"]

    # Generate cell SVG
    cell_svg = _generate_cell_svg(
        cells_data,
        type_colors,
        domain.width_um,
        domain.height_um,
        svg_width=700,
        svg_height=700,
    )

    # Create Vega-Lite specs for charts (keep these for statistics)
    env_spec = {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "title": "Oxygen Distribution",
        "width": 280,
        "height": 280,
        "data": {"values": env_data},
        "mark": "rect",
        "encoding": {
            "x": {"field": "x", "type": "quantitative", "bin": {"maxbins": 40}, "title": "X (um)"},
            "y": {"field": "y", "type": "quantitative", "bin": {"maxbins": 40}, "title": "Y (um)"},
            "color": {"field": "O2", "type": "quantitative", "aggregate": "mean",
                      "scale": {"scheme": "viridis"}, "title": "O2 (mmHg)"},
        },
    }

    vegf_spec = {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "title": "VEGF Distribution",
        "width": 280,
        "height": 280,
        "data": {"values": env_data},
        "mark": "rect",
        "encoding": {
            "x": {"field": "x", "type": "quantitative", "bin": {"maxbins": 40}, "title": "X (um)"},
            "y": {"field": "y", "type": "quantitative", "bin": {"maxbins": 40}, "title": "Y (um)"},
            "color": {"field": "VEGF", "type": "quantitative", "aggregate": "mean",
                      "scale": {"scheme": "reds"}, "title": "VEGF (nM)"},
        },
    }

    hist_spec = {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "title": "Cell Count by Type",
        "width": 350,
        "height": 220,
        "data": {"values": [{"type": k, "count": v} for k, v in sorted(cell_counts.items())]},
        "mark": "bar",
        "encoding": {
            "x": {"field": "type", "type": "nominal", "title": "Cell Type", "sort": "-y"},
            "y": {"field": "count", "type": "quantitative", "title": "Count"},
            "color": {"field": "type", "type": "nominal",
                      "scale": {"domain": list(type_colors.keys()), "range": list(type_colors.values())},
                      "legend": None},
        },
    }

    ecDNA_data = [c for c in cells_data if c["ecDNA_count"] > 0]
    ecDNA_spec = {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "title": "ecDNA Distribution",
        "width": 350,
        "height": 220,
        "data": {"values": ecDNA_data},
        "mark": "boxplot",
        "encoding": {
            "x": {"field": "cell_type_name", "type": "nominal", "title": "Cell Type"},
            "y": {"field": "ecDNA_count", "type": "quantitative", "title": "ecDNA Count"},
            "color": {"field": "cell_type_name", "type": "nominal",
                      "scale": {"domain": list(type_colors.keys()), "range": list(type_colors.values())},
                      "legend": None},
        },
    }

    # Generate legend SVG
    legend_items = []
    for i, (tname, color) in enumerate(type_colors.items()):
        y = 20 + i * 25
        shape_path = type_shapes.get(tname, "")
        if shape_path:
            legend_items.append(
                f'<g transform="translate(15, {y})">'
                f'<svg viewBox="-1 -1 2 2" width="18" height="18" style="overflow:visible">'
                f'<path d="{shape_path}" fill="{color}" stroke="#333" stroke-width="0.1"/>'
                f'</svg>'
                f'<text x="25" y="14" fill="#eee" font-size="13">{tname}</text></g>'
            )
        else:
            legend_items.append(
                f'<g transform="translate(15, {y})">'
                f'<circle cx="9" cy="9" r="8" fill="{color}" stroke="#333"/>'
                f'<text x="25" y="14" fill="#eee" font-size="13">{tname}</text></g>'
            )

    legend_svg = (
        f'<svg width="180" height="{len(type_colors) * 25 + 30}">'
        + "\n".join(legend_items)
        + "</svg>"
    )

    # Generate HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
    <style>
        :root {{
            --bg-dark: #1a1a2e;
            --bg-card: #16213e;
            --text-primary: #eee;
            --text-secondary: #aaa;
            --accent: #e94560;
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
        h2 {{ color: var(--accent); margin-top: 30px; font-size: 1.3em; }}
        .card {{
            background: var(--bg-card);
            border-radius: 8px;
            padding: 20px;
            margin: 20px 0;
            box-shadow: 0 4px 6px rgba(0,0,0,0.3);
        }}
        .grid {{ display: grid; gap: 20px; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); }}
        .stat-box {{
            background: #0f3460;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
        }}
        .stat-value {{ font-size: 2em; font-weight: bold; color: var(--accent); }}
        .stat-label {{ color: var(--text-secondary); font-size: 0.9em; }}
        .vis-container {{ display: flex; justify-content: center; margin: 15px 0; }}
        .spatial-view {{
            display: flex;
            gap: 20px;
            align-items: flex-start;
            flex-wrap: wrap;
            justify-content: center;
        }}
        .cell-svg-container {{
            border: 1px solid #333;
            border-radius: 8px;
            overflow: hidden;
        }}
        .legend-container {{
            background: #0f3460;
            padding: 15px;
            border-radius: 8px;
        }}
        table {{ width: 100%; border-collapse: collapse; margin: 15px 0; }}
        th, td {{ padding: 10px; text-align: left; border-bottom: 1px solid #333; }}
        th {{ background: #0f3460; }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #333;
            text-align: center;
            color: var(--text-secondary);
        }}
        .note {{ font-size: 0.85em; color: var(--text-secondary); margin-top: 10px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>{title}</h1>

        <div class="card">
            <h2>Simulation Summary</h2>
            <div class="grid">
                <div class="stat-box">
                    <div class="stat-value">{len(cells_data)}</div>
                    <div class="stat-label">Total Cells</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{cell_counts.get('Tumor', 0)}</div>
                    <div class="stat-label">Tumor Cells</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{total_ecDNA}</div>
                    <div class="stat-label">Total ecDNA</div>
                </div>
                <div class="stat-box">
                    <div class="stat-value">{time_hr:.0f}h</div>
                    <div class="stat-label">Simulation Time</div>
                </div>
            </div>
            <table>
                <tr><th>Parameter</th><th>Value</th></tr>
                <tr><td>Domain Size</td><td>{domain.width_um} x {domain.height_um} um</td></tr>
                <tr><td>Grid Resolution</td><td>{domain.env_grid_um} um ({domain.nx} x {domain.ny})</td></tr>
                <tr><td>Mean O2</td><td>{float(env.O2.mean()):.1f} mmHg</td></tr>
                <tr><td>Max VEGF</td><td>{float(env.VEGF.max()):.3f} nM</td></tr>
            </table>
        </div>

        <div class="card">
            <h2>Spatial Distribution (Shape-Accurate)</h2>
            <div class="spatial-view">
                <div class="cell-svg-container">
                    {cell_svg}
                </div>
                <div class="legend-container">
                    <strong style="color: #e94560;">Cell Types</strong>
                    {legend_svg}
                    <div class="note">
                        Cells with >5 ecDNA highlighted in gold outline.<br/>
                        Hover over cells for details.
                    </div>
                </div>
            </div>
        </div>

        <div class="card">
            <h2>Environment</h2>
            <div class="grid">
                <div class="vis-container"><div id="o2-vis"></div></div>
                <div class="vis-container"><div id="vegf-vis"></div></div>
            </div>
        </div>

        <div class="card">
            <h2>Cell Statistics</h2>
            <div class="grid">
                <div class="vis-container"><div id="hist-vis"></div></div>
                <div class="vis-container"><div id="ecdna-vis"></div></div>
            </div>
        </div>

        <div class="footer">
            Generated by CAUSANTA | Cells rendered with shape-accurate collision detection
        </div>
    </div>

    <script>
        const envSpec = {json.dumps(env_spec)};
        const vegfSpec = {json.dumps(vegf_spec)};
        const histSpec = {json.dumps(hist_spec)};
        const ecdnaSpec = {json.dumps(ecDNA_spec)};

        const opts = {{
            theme: 'dark',
            actions: {{ export: true, source: false, editor: false }},
            config: {{
                background: '#16213e',
                title: {{ color: '#eee' }},
                axis: {{ labelColor: '#aaa', titleColor: '#eee', gridColor: '#333', domainColor: '#666' }},
                legend: {{ labelColor: '#eee', titleColor: '#eee' }},
                view: {{ stroke: 'transparent' }}
            }}
        }};

        vegaEmbed('#o2-vis', envSpec, opts);
        vegaEmbed('#vegf-vis', vegfSpec, opts);
        vegaEmbed('#hist-vis', histSpec, opts);
        vegaEmbed('#ecdna-vis', ecdnaSpec, opts);
    </script>
</body>
</html>"""

    output_path = Path(output_path)
    with open(output_path, "w") as f:
        f.write(html)

    return output_path


def generate_vega_json(
    population: CellPopulation,
    env: EnvironmentFields,
    domain: Domain,
    config: SimulationConfig,
    time_hr: float,
    output_path: Path,
) -> Path:
    """Generate a standalone Vega JSON spec with embedded data.

    Note: This uses standard Vega-Lite marks. For shape-accurate rendering,
    use generate_standalone_html() which produces SVG directly.
    """
    cells_data = cells_to_json(population, config)
    env_data = environment_to_json(env, domain, downsample=2)

    # Map cell types to Vega shape symbols (approximations)
    vega_shapes = {
        "Neuron": "triangle-up",
        "Astrocyte": "cross",
        "Oligodendrocyte": "circle",
        "Microglia": "diamond",
        "Endothelial": "square",
        "Pericyte": "square",
        "Tumor": "circle",
        "RecruitedImmune": "circle",
        "Necrotic": "square",
    }

    # Add shape field to cells
    for cell in cells_data:
        cell["vega_shape"] = vega_shapes.get(cell["cell_type_name"], "circle")

    # Build color mapping
    type_colors = {}
    for tid, tc in config.cell_types.items():
        type_colors[tc.type_name] = tc.color

    spec = {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "title": f"CAUSANTA Simulation (t = {time_hr:.0f} hours)",
        "width": 700,
        "height": 700,
        "layer": [
            # Environment heatmap layer
            {
                "data": {"values": env_data},
                "mark": {"type": "rect", "opacity": 0.4},
                "encoding": {
                    "x": {"field": "x", "type": "quantitative", "bin": {"maxbins": 50}},
                    "y": {"field": "y", "type": "quantitative", "bin": {"maxbins": 50}},
                    "color": {
                        "field": "O2", "type": "quantitative", "aggregate": "mean",
                        "scale": {"scheme": "viridis"}, "legend": {"title": "O2 (mmHg)"},
                    },
                },
            },
            # Cell point layer with shapes
            {
                "data": {"values": cells_data},
                "mark": {"type": "point", "filled": True, "opacity": 0.9, "stroke": "#333", "strokeWidth": 0.5},
                "encoding": {
                    "x": {"field": "x", "type": "quantitative", "title": "X (um)"},
                    "y": {"field": "y", "type": "quantitative", "title": "Y (um)"},
                    "color": {
                        "field": "cell_type_name", "type": "nominal",
                        "scale": {"domain": list(type_colors.keys()), "range": list(type_colors.values())},
                        "legend": {"title": "Cell Type"},
                    },
                    "shape": {
                        "field": "vega_shape", "type": "nominal",
                        "legend": None,
                    },
                    "size": {
                        "field": "diameter", "type": "quantitative",
                        "scale": {"range": [20, 400]},
                        "legend": {"title": "Diameter"},
                    },
                    "tooltip": [
                        {"field": "cell_id", "title": "ID"},
                        {"field": "cell_type_name", "title": "Type"},
                        {"field": "ecDNA_count", "title": "ecDNA"},
                        {"field": "diameter", "title": "Diameter"},
                        {"field": "O2_local", "title": "O2"},
                    ],
                },
            },
        ],
    }

    output_path = Path(output_path)
    with open(output_path, "w") as f:
        json.dump(spec, f, indent=2)

    return output_path
