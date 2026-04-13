"""Generate animated visualization of CAUSANTA simulation.

Creates an HTML file with JavaScript-based animation showing tumor growth over time.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np


def read_cells_tsv(path: Path) -> list[dict]:
    """Read cells from TSV file."""
    cells = []
    with open(path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            values = line.strip().split('\t')
            cell = {}
            for i, col in enumerate(header):
                val = values[i] if i < len(values) else ''
                # Convert numeric fields
                if col in ('cell_id', 'parent_id', 'cell_type', 'x', 'y', 'ecDNA_count', 'generation'):
                    cell[col] = int(val) if val else 0
                elif col in ('diameter', 'angle', 'size_scale', 'O2_local', 'glucose_local'):
                    cell[col] = float(val) if val else 0.0
                else:
                    cell[col] = val
            cells.append(cell)
    return cells


def generate_frame_svg(
    cells: list[dict],
    type_colors: dict[str, str],
    type_shapes: dict[str, str],
    domain_width: float,
    domain_height: float,
    time_hr: int,
    svg_size: int = 500,
) -> str:
    """Generate SVG string for a single frame."""
    scale = svg_size / max(domain_width, domain_height)

    parts = []

    # Render cells sorted: normal cells first, tumor on top
    sorted_cells = sorted(cells, key=lambda c: (c.get('cell_type_name', '') == 'Tumor', c.get('cell_type', 0)))

    for cell in sorted_cells:
        x = cell.get('x', 0) * scale
        y = cell.get('y', 0) * scale
        diameter = cell.get('diameter', 10) * scale
        cell_scale = diameter / 2
        angle_deg = math.degrees(cell.get('angle', 0))
        type_name = cell.get('cell_type_name', 'Unknown')
        color = type_colors.get(type_name, '#888')
        shape_path = type_shapes.get(type_name, '')

        # Highlight tumor cells
        if type_name == 'Tumor':
            parts.append(
                f'<circle cx="{x}" cy="{y}" r="{cell_scale + 2}" '
                f'fill="none" stroke="#FFD700" stroke-width="2"/>'
            )

        if shape_path:
            parts.append(
                f'<path d="{shape_path}" fill="{color}" stroke="#222" stroke-width="{0.08/cell_scale:.4f}" '
                f'transform="translate({x},{y}) rotate({angle_deg}) scale({cell_scale})"/>'
            )
        else:
            parts.append(
                f'<circle cx="{x}" cy="{y}" r="{cell_scale}" fill="{color}" stroke="#222" stroke-width="0.5"/>'
            )

    return '\n'.join(parts)


def generate_movie_html(
    output_dir: Path,
    movie_path: Path,
    domain_width: float = 1000,
    domain_height: float = 1000,
    fps: int = 4,
) -> Path:
    """Generate animated HTML from simulation output.

    Args:
        output_dir: Directory containing cells_t*.tsv files
        movie_path: Output path for the HTML movie
        domain_width: Domain width in micrometers
        domain_height: Domain height in micrometers
        fps: Frames per second for playback

    Returns:
        Path to generated HTML file.
    """
    # Find all cell TSV files
    cell_files = sorted(output_dir.glob('cells_t*.tsv'))
    if not cell_files:
        raise ValueError(f"No cells_t*.tsv files found in {output_dir}")

    print(f"Found {len(cell_files)} frames")

    # Load params for colors and shapes
    params_path = output_dir / 'params.json'
    type_colors = {}
    type_shapes = {}
    if params_path.exists():
        with open(params_path) as f:
            params = json.load(f)
        for tid, tc in params.get('cell_types', {}).items():
            type_colors[tc.get('type_name', '')] = tc.get('color', '#888')
            type_shapes[tc.get('type_name', '')] = tc.get('shape_path', '')

    # Default colors if not found
    if not type_colors:
        type_colors = {
            'Neuron': '#4477AA',
            'Astrocyte': '#66CCEE',
            'Oligodendrocyte': '#228833',
            'Microglia': '#EE6677',
            'Endothelial': '#AA3377',
            'Pericyte': '#BBBBBB',
            'Tumor': '#CCBB44',
            'RecruitedImmune': '#EE8866',
            'Necrotic': '#555555',
        }

    # Generate frames
    svg_size = 600
    frames = []
    tumor_counts = []

    for i, cell_file in enumerate(cell_files):
        # Extract time from filename
        time_str = cell_file.stem.replace('cells_t', '')
        time_hr = int(time_str)

        cells = read_cells_tsv(cell_file)
        tumor_count = sum(1 for c in cells if c.get('cell_type_name') == 'Tumor')
        tumor_counts.append(tumor_count)

        frame_svg = generate_frame_svg(
            cells, type_colors, type_shapes,
            domain_width, domain_height, time_hr, svg_size
        )
        frames.append({
            'time': time_hr,
            'svg': frame_svg,
            'tumor_count': tumor_count,
            'total_cells': len(cells),
        })

        if i % 20 == 0:
            print(f"  Frame {i+1}/{len(cell_files)}: t={time_hr}h, tumor={tumor_count}")

    # Generate HTML
    frame_interval = 1000 // fps  # milliseconds per frame

    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>CAUSANTA Tumor Growth Animation</title>
    <style>
        :root {{
            --bg-dark: #1a1a2e;
            --bg-card: #16213e;
            --text-primary: #eee;
            --accent: #e94560;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: var(--bg-dark);
            color: var(--text-primary);
            margin: 0;
            padding: 20px;
            display: flex;
            flex-direction: column;
            align-items: center;
        }}
        h1 {{ color: var(--accent); margin-bottom: 10px; }}
        .container {{
            background: var(--bg-card);
            border-radius: 12px;
            padding: 20px;
            box-shadow: 0 4px 20px rgba(0,0,0,0.4);
        }}
        .controls {{
            display: flex;
            gap: 15px;
            align-items: center;
            margin-bottom: 15px;
        }}
        button {{
            background: var(--accent);
            border: none;
            color: white;
            padding: 10px 20px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 14px;
        }}
        button:hover {{ opacity: 0.9; }}
        .info {{
            background: #0f3460;
            padding: 10px 20px;
            border-radius: 6px;
            font-family: monospace;
        }}
        .slider-container {{
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        input[type="range"] {{
            width: 200px;
        }}
        #canvas {{
            background: #1a1a2e;
            border-radius: 8px;
            border: 1px solid #333;
        }}
        .legend {{
            display: flex;
            gap: 15px;
            margin-top: 15px;
            flex-wrap: wrap;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 5px;
            font-size: 12px;
        }}
        .legend-color {{
            width: 16px;
            height: 16px;
            border-radius: 3px;
            border: 1px solid #333;
        }}
    </style>
</head>
<body>
    <h1>Tumor Growth Simulation</h1>

    <div class="container">
        <div class="controls">
            <button id="playBtn" onclick="togglePlay()">Play</button>
            <button onclick="stepBack()">Step Back</button>
            <button onclick="stepForward()">Step Forward</button>
            <div class="slider-container">
                <span>Speed:</span>
                <input type="range" id="speedSlider" min="1" max="20" value="{fps}"
                       onchange="updateSpeed(this.value)">
                <span id="speedLabel">{fps} fps</span>
            </div>
        </div>

        <div class="info">
            <span>Time: <strong id="timeLabel">0h</strong></span> |
            <span>Tumor cells: <strong id="tumorLabel">0</strong></span> |
            <span>Total cells: <strong id="totalLabel">0</strong></span> |
            <span>Frame: <strong id="frameLabel">1/{len(frames)}</strong></span>
        </div>

        <svg id="canvas" width="{svg_size}" height="{svg_size}"
             viewBox="0 0 {svg_size} {svg_size}">
            <!-- Grid -->
            <g stroke="#333" stroke-width="0.5" opacity="0.3">
                {"".join(f'<line x1="{i*svg_size/10}" y1="0" x2="{i*svg_size/10}" y2="{svg_size}"/>' for i in range(11))}
                {"".join(f'<line x1="0" y1="{i*svg_size/10}" x2="{svg_size}" y2="{i*svg_size/10}"/>' for i in range(11))}
            </g>
            <g id="cellsLayer"></g>
            <!-- Scale bar -->
            <g transform="translate(20, {svg_size - 25})">
                <rect width="{svg_size/10}" height="4" fill="#eee"/>
                <text x="{svg_size/20}" y="18" text-anchor="middle" fill="#eee" font-size="11">100 um</text>
            </g>
        </svg>

        <div class="legend">
            <div class="legend-item">
                <div class="legend-color" style="background: #CCBB44; box-shadow: 0 0 5px #FFD700;"></div>
                <span>Tumor</span>
            </div>
            <div class="legend-item">
                <div class="legend-color" style="background: #4477AA;"></div>
                <span>Neuron</span>
            </div>
            <div class="legend-item">
                <div class="legend-color" style="background: #66CCEE;"></div>
                <span>Astrocyte</span>
            </div>
            <div class="legend-item">
                <div class="legend-color" style="background: #EE6677;"></div>
                <span>Microglia</span>
            </div>
            <div class="legend-item">
                <div class="legend-color" style="background: #AA3377;"></div>
                <span>Endothelial</span>
            </div>
        </div>
    </div>

    <script>
        const frames = {json.dumps([{'time': f['time'], 'svg': f['svg'], 'tumor': f['tumor_count'], 'total': f['total_cells']} for f in frames])};

        let currentFrame = 0;
        let playing = false;
        let interval = null;
        let fps = {fps};

        const cellsLayer = document.getElementById('cellsLayer');
        const timeLabel = document.getElementById('timeLabel');
        const tumorLabel = document.getElementById('tumorLabel');
        const totalLabel = document.getElementById('totalLabel');
        const frameLabel = document.getElementById('frameLabel');
        const playBtn = document.getElementById('playBtn');
        const speedLabel = document.getElementById('speedLabel');

        function showFrame(idx) {{
            if (idx < 0 || idx >= frames.length) return;
            currentFrame = idx;
            const frame = frames[idx];
            cellsLayer.innerHTML = frame.svg;
            timeLabel.textContent = frame.time + 'h';
            tumorLabel.textContent = frame.tumor;
            totalLabel.textContent = frame.total;
            frameLabel.textContent = (idx + 1) + '/' + frames.length;
        }}

        function togglePlay() {{
            playing = !playing;
            playBtn.textContent = playing ? 'Pause' : 'Play';
            if (playing) {{
                interval = setInterval(() => {{
                    currentFrame = (currentFrame + 1) % frames.length;
                    showFrame(currentFrame);
                }}, 1000 / fps);
            }} else {{
                clearInterval(interval);
            }}
        }}

        function stepForward() {{
            if (playing) togglePlay();
            showFrame((currentFrame + 1) % frames.length);
        }}

        function stepBack() {{
            if (playing) togglePlay();
            showFrame((currentFrame - 1 + frames.length) % frames.length);
        }}

        function updateSpeed(value) {{
            fps = parseInt(value);
            speedLabel.textContent = fps + ' fps';
            if (playing) {{
                clearInterval(interval);
                interval = setInterval(() => {{
                    currentFrame = (currentFrame + 1) % frames.length;
                    showFrame(currentFrame);
                }}, 1000 / fps);
            }}
        }}

        // Keyboard controls
        document.addEventListener('keydown', (e) => {{
            if (e.code === 'Space') {{ togglePlay(); e.preventDefault(); }}
            if (e.code === 'ArrowRight') stepForward();
            if (e.code === 'ArrowLeft') stepBack();
        }});

        // Show first frame
        showFrame(0);
    </script>
</body>
</html>'''

    movie_path = Path(movie_path)
    with open(movie_path, 'w') as f:
        f.write(html)

    print(f"\nGenerated movie: {movie_path}")
    print(f"  Frames: {len(frames)}")
    print(f"  Duration: {len(frames) / fps:.1f} seconds at {fps} fps")
    print(f"  Tumor growth: {tumor_counts[0]} -> {tumor_counts[-1]} cells")

    return movie_path


if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: python -m causanta.simulate.movie <output_dir> [movie.html]")
        sys.exit(1)

    output_dir = Path(sys.argv[1])
    movie_path = Path(sys.argv[2]) if len(sys.argv) > 2 else output_dir / 'movie.html'

    generate_movie_html(output_dir, movie_path)
