"""Vega.js specification generation for CAUSANTA visualization.

Generates a self-contained Vega JSON spec that renders cell agents
as SVG path marks over an environment heatmap, with interactive
time slider, layer toggles, and cell type filtering.
"""

from __future__ import annotations

import json
from pathlib import Path

from .config import SimulationConfig


def generate_vega_spec(
    output_dir: Path,
    config: SimulationConfig,
    total_steps: int,
    output_interval: int = 1,
) -> dict:
    """Generate a Vega.js spec for visualizing simulation output.

    The spec loads cells and environment TSV files and renders them
    as an interactive visualization.
    """
    w = config.domain.width_um
    h = config.domain.height_um

    # Build color scale from cell type configs
    type_colors = {}
    type_names = {}
    for tid, tc in sorted(config.cell_types.items()):
        type_colors[tid] = tc.color
        type_names[tid] = tc.type_name

    color_domain = [type_names[tid] for tid in sorted(type_colors.keys())]
    color_range = [type_colors[tid] for tid in sorted(type_colors.keys())]

    # Time steps available
    steps = list(range(0, total_steps + 1, output_interval))

    spec = {
        "$schema": "https://vega.github.io/schema/vega/v5.json",
        "width": 800,
        "height": 800,
        "padding": 10,
        "autosize": "none",

        "signals": [
            {
                "name": "timeStep",
                "value": 0,
                "bind": {
                    "input": "range",
                    "min": 0,
                    "max": total_steps,
                    "step": output_interval,
                    "name": "Time (hours): ",
                },
            },
            {
                "name": "envChannel",
                "value": "O2",
                "bind": {
                    "input": "select",
                    "options": ["O2", "glucose", "VEGF", "lactate", "pH", "ECM_density", "vascular_density"],
                    "name": "Environment: ",
                },
            },
            {
                "name": "envOpacity",
                "value": 0.4,
                "bind": {
                    "input": "range",
                    "min": 0,
                    "max": 1,
                    "step": 0.05,
                    "name": "Env opacity: ",
                },
            },
            {
                "name": "showCells",
                "value": True,
                "bind": {"input": "checkbox", "name": "Show cells: "},
            },
            {
                "name": "showEnv",
                "value": True,
                "bind": {"input": "checkbox", "name": "Show environment: "},
            },
            {
                "name": "cellScale",
                "value": 1.0,
                "bind": {
                    "input": "range",
                    "min": 0.2,
                    "max": 3.0,
                    "step": 0.1,
                    "name": "Cell scale: ",
                },
            },
        ],

        "data": [
            {
                "name": "cells",
                "url": {"signal": "'cells_t' + pad(timeStep, 6, '0', 'left') + '.tsv'"},
                "format": {"type": "tsv"},
            },
            {
                "name": "environment",
                "url": {"signal": "'environment_t' + pad(timeStep, 6, '0', 'left') + '.tsv'"},
                "format": {"type": "tsv"},
            },
        ],

        "scales": [
            {
                "name": "x",
                "type": "linear",
                "domain": [0, w],
                "range": "width",
            },
            {
                "name": "y",
                "type": "linear",
                "domain": [0, h],
                "range": "height",
            },
            {
                "name": "color",
                "type": "ordinal",
                "domain": color_domain,
                "range": color_range,
            },
            {
                "name": "envColor",
                "type": "linear",
                "domain": [0, 1],
                "range": {"scheme": "viridis"},
            },
        ],

        "marks": [
            # Environment heatmap layer
            {
                "type": "rect",
                "from": {"data": "environment"},
                "encode": {
                    "enter": {
                        "width": {"value": 800 / (w / config.domain.env_grid_um)},
                        "height": {"value": 800 / (h / config.domain.env_grid_um)},
                    },
                    "update": {
                        "x": {"scale": "x", "field": "x_um"},
                        "y": {"scale": "y", "field": "y_um"},
                        "fill": {"scale": "envColor", "signal": "datum[envChannel]"},
                        "opacity": {"signal": "showEnv ? envOpacity : 0"},
                    },
                },
            },
            # Cell marks layer
            {
                "type": "symbol",
                "from": {"data": "cells"},
                "encode": {
                    "update": {
                        "x": {"scale": "x", "field": "x"},
                        "y": {"scale": "y", "field": "y"},
                        "size": {
                            "signal": "pow(datum.size_scale * cellScale * 5, 2)",
                        },
                        "fill": {"scale": "color", "field": "cell_type_name"},
                        "opacity": {
                            "signal": "showCells ? (datum.cell_type == 8 ? 0.3 : (datum.is_hypoxic == 'True' ? 0.7 : 1.0)) : 0",
                        },
                        "shape": {"value": "circle"},
                        "stroke": {"value": "#333"},
                        "strokeWidth": {"value": 0.5},
                        "tooltip": {
                            "signal": "{'Cell ID': datum.cell_id, 'Type': datum.cell_type_name, 'ecDNA': datum.ecDNA_count, 'O2': datum.O2_local, 'Phase': datum.cell_cycle_phase, 'Generation': datum.generation}",
                        },
                    },
                },
            },
        ],

        "legends": [
            {
                "fill": "color",
                "title": "Cell Type",
                "orient": "right",
            },
        ],
    }

    return spec


def write_vega_spec(
    output_dir: Path,
    config: SimulationConfig,
    total_steps: int,
    output_interval: int = 1,
) -> Path:
    """Generate and write Vega spec to output directory."""
    spec = generate_vega_spec(output_dir, config, total_steps, output_interval)
    spec_path = output_dir / "visualization.vg.json"
    with open(spec_path, "w") as f:
        json.dump(spec, f, indent=2)
    return spec_path


def write_html_viewer(output_dir: Path) -> Path:
    """Write a self-contained HTML file that loads the Vega spec."""
    html_content = """<!DOCTYPE html>
<html>
<head>
    <title>CAUSANTA Visualization</title>
    <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
    <style>
        body { font-family: sans-serif; margin: 20px; background: #1a1a2e; color: #eee; }
        h1 { color: #e94560; }
        #vis { background: #16213e; padding: 10px; border-radius: 8px; }
    </style>
</head>
<body>
    <h1>CAUSANTA Tissue Simulation</h1>
    <div id="vis"></div>
    <script>
        vegaEmbed('#vis', 'visualization.vg.json', {
            theme: 'dark',
            actions: { export: true, source: true, editor: true }
        }).catch(console.error);
    </script>
</body>
</html>"""
    html_path = output_dir / "index.html"
    with open(html_path, "w") as f:
        f.write(html_content)
    return html_path
