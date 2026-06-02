"""Nature Methods styling for publication figures.

Follows Nature Methods guidelines:
- Single column: 88mm (3.46 in)
- Double column: 180mm (7.09 in)
- Maximum height: 230mm (9.06 in)
- Font: Arial/Helvetica, 5-7pt for labels
- Resolution: 300 DPI minimum
- Colors: Accessible, colorblind-friendly
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np


# Nature Methods color palette (colorblind-friendly)
NATURE_COLORS = {
    # Primary colors
    "blue": "#0077BB",
    "orange": "#EE7733",
    "green": "#009988",
    "red": "#CC3311",
    "purple": "#AA3377",
    "cyan": "#33BBEE",
    "gray": "#BBBBBB",

    # Semantic colors
    "ecDNA": "#0077BB",      # Blue - instrument
    "EGFR": "#EE7733",       # Orange - treatment
    "outcome": "#009988",    # Green - outcomes
    "hypoxia": "#CC3311",    # Red - confounder
    "IV": "#AA3377",         # Purple - IV estimates
    "OLS": "#BBBBBB",        # Gray - OLS (biased)
    "true": "#000000",       # Black - ground truth

    # Region colors
    "core": "#CC3311",
    "margin": "#EE7733",
    "infiltrating": "#0077BB",

    # Confidence intervals
    "ci_fill": "#0077BB33",  # Transparent blue
}

# Figure sizes (in inches)
NATURE_FIGSIZE = {
    "single": (3.46, 3.0),      # Single column
    "single_tall": (3.46, 5.0),
    "double": (7.09, 4.0),      # Double column
    "double_tall": (7.09, 7.0),
    "full_page": (7.09, 9.0),
}


def apply_nature_style() -> None:
    """Apply Nature Methods styling to matplotlib."""
    plt.style.use("default")

    mpl.rcParams.update({
        # Font
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,

        # Lines
        "lines.linewidth": 1.0,
        "lines.markersize": 4,

        # Axes
        "axes.linewidth": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelpad": 4,
        "axes.titlepad": 8,

        # Ticks
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.direction": "out",
        "ytick.direction": "out",

        # Legend
        "legend.frameon": False,
        "legend.borderpad": 0.5,

        # Figure
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,

        # Grid
        "axes.grid": False,
    })


def save_figure(
    fig: plt.Figure,
    path: str | Path,
    formats: list[str] | None = None,
    dpi: int = 300,
) -> list[Path]:
    """Save figure in publication-ready formats.

    Args:
        fig: Matplotlib figure
        path: Output path (without extension)
        formats: List of formats to save (default: ["pdf", "png", "svg"])
        dpi: Resolution for raster formats

    Returns:
        List of saved file paths
    """
    if formats is None:
        formats = ["pdf", "png", "svg"]

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    saved = []
    for fmt in formats:
        out_path = path.with_suffix(f".{fmt}")
        fig.savefig(
            out_path,
            format=fmt,
            dpi=dpi if fmt in ["png", "tiff"] else None,
            bbox_inches="tight",
            pad_inches=0.05,
            facecolor="white",
            edgecolor="none",
        )
        saved.append(out_path)

    return saved


def add_panel_label(
    ax: plt.Axes,
    label: str,
    x: float = -0.15,
    y: float = 1.05,
    fontsize: int = 10,
    fontweight: str = "bold",
) -> None:
    """Add panel label (a, b, c, etc.) to axes.

    Args:
        ax: Matplotlib axes
        label: Panel label (e.g., "a", "b")
        x: X position in axes coordinates
        y: Y position in axes coordinates
        fontsize: Font size
        fontweight: Font weight
    """
    ax.text(
        x, y, label,
        transform=ax.transAxes,
        fontsize=fontsize,
        fontweight=fontweight,
        va="top",
        ha="left",
    )


def create_colorbar(
    fig: plt.Figure,
    ax: plt.Axes,
    mappable: Any,
    label: str,
    orientation: str = "vertical",
    shrink: float = 0.8,
) -> plt.colorbar:
    """Create a Nature-styled colorbar.

    Args:
        fig: Matplotlib figure
        ax: Axes to attach colorbar to
        mappable: The mappable (e.g., imshow result)
        label: Colorbar label
        orientation: "vertical" or "horizontal"
        shrink: Shrink factor

    Returns:
        Colorbar object
    """
    cbar = fig.colorbar(
        mappable,
        ax=ax,
        orientation=orientation,
        shrink=shrink,
        pad=0.02,
    )
    cbar.set_label(label, fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    cbar.outline.set_linewidth(0.5)
    return cbar


def significance_stars(p_value: float) -> str:
    """Convert p-value to significance stars.

    Args:
        p_value: Statistical p-value

    Returns:
        Significance stars string
    """
    if p_value < 0.001:
        return "***"
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    else:
        return "ns"


def format_ci(estimate: float, ci_lower: float, ci_upper: float, decimals: int = 3) -> str:
    """Format estimate with confidence interval.

    Args:
        estimate: Point estimate
        ci_lower: Lower CI bound
        ci_upper: Upper CI bound
        decimals: Number of decimal places

    Returns:
        Formatted string like "0.100 [0.080, 0.120]"
    """
    fmt = f".{decimals}f"
    return f"{estimate:{fmt}} [{ci_lower:{fmt}}, {ci_upper:{fmt}}]"


def truncate_colormap(cmap_name: str, minval: float = 0.0, maxval: float = 1.0, n: int = 256):
    """Truncate a colormap to a subset of its range.

    Args:
        cmap_name: Name of colormap
        minval: Minimum value (0-1)
        maxval: Maximum value (0-1)
        n: Number of colors

    Returns:
        New colormap
    """
    cmap = plt.cm.get_cmap(cmap_name)
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        f"trunc({cmap_name},{minval:.2f},{maxval:.2f})",
        cmap(np.linspace(minval, maxval, n))
    )
    return new_cmap
