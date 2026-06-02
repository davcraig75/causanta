#!/usr/bin/env python3
"""Generate publication-quality Figure 8: CAUSANTA simulation architecture.

A single 3-panel figure summarising the simulator:
  A. Multi-scale model state (tissue / cellular / molecular)
  B. Per-hour 6-phase simulation loop
  C. Mechanistic causal model encoded inside each tumor cell

Outputs (PNG @ 600 dpi and PDF) are written to:
  figures/figure8_simulation_architecture.{png,pdf}
  docs/manuscript/figures/figure8_simulation_architecture.{png,pdf}
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np

# ---- publication style ----
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8.5,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.labelsize": 8.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# ---- Tol-vibrant + cell-type palette (colorblind safe) ----
C_TISSUE = "#0077BB"
C_CELL = "#009988"
C_MOL = "#EE7733"
C_INPUT = "#33BBEE"
C_OUTPUT = "#CC6677"
C_GRAY = "#888888"
C_LIGHT = "#F4F5F7"
C_PANEL = "#FFFFFF"
C_INK = "#222222"
C_LINK = "#555555"

# Cell-type swatches (from default.json)
CELL_TYPE_PALETTE = [
    ("Neuron", "#4477AA"),
    ("Astrocyte", "#66CCEE"),
    ("Oligodendr.", "#228833"),
    ("Microglia", "#EE6677"),
    ("Endothelial", "#AA3377"),
    ("Pericyte", "#BBBBBB"),
    ("Tumor", "#CCBB44"),
    ("Immune", "#EE8866"),
    ("Necrotic", "#555555"),
]


def box(ax, x, y, w, h, text, *, fc=C_PANEL, ec=C_INK, lw=1.0,
        fontsize=8.5, fontweight="normal", color=C_INK, pad=0.0,
        rounded=0.04, ha="center", va="center"):
    """Rounded rectangle with centered text. pad=0 keeps interior at full size."""
    patch = FancyBboxPatch(
        (x + pad, y + pad), w - 2 * pad, h - 2 * pad,
        boxstyle=f"round,pad=0,rounding_size={rounded}",
        facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2,
    )
    ax.add_patch(patch)
    if text:
        ax.text(x + w / 2, y + h / 2, text, ha=ha, va=va,
                fontsize=fontsize, fontweight=fontweight, color=color, zorder=3)


def arrow(ax, xy_from, xy_to, *, color=C_LINK, lw=1.2, style="-",
          mut=12, alpha=1.0):
    """Single straight arrow with adjustable head."""
    p = FancyArrowPatch(
        xy_from, xy_to,
        arrowstyle="-|>", mutation_scale=mut,
        color=color, linewidth=lw, linestyle=style,
        alpha=alpha, zorder=4, shrinkA=2, shrinkB=2,
    )
    ax.add_patch(p)


def bidir_arrow(ax, xy_from, xy_to, *, color=C_LINK, lw=1.2, mut=10):
    p = FancyArrowPatch(
        xy_from, xy_to,
        arrowstyle="<|-|>", mutation_scale=mut,
        color=color, linewidth=lw, zorder=4, shrinkA=2, shrinkB=2,
    )
    ax.add_patch(p)


def label(ax, x, y, text, *, fontsize=8, color=C_INK, fontweight="normal",
          ha="center", va="center", rotation=0, alpha=1.0):
    ax.text(x, y, text, ha=ha, va=va, fontsize=fontsize, color=color,
            fontweight=fontweight, rotation=rotation, alpha=alpha, zorder=5)


def panel_setup(ax, xlim, ylim):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def panel_letter(ax, text, x=0.005, y=0.985):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=12,
            fontweight="bold", ha="left", va="top", color=C_INK)


# ===================== PANEL A: multi-scale =====================
def panel_a(ax):
    panel_setup(ax, (0, 10), (0, 5))
    panel_letter(ax, "A")
    ax.set_title("Multi-scale model state", loc="left", x=0.04, pad=2)

    # background tint
    ax.add_patch(mpatches.Rectangle((0, 0), 10, 5, facecolor=C_LIGHT,
                                    edgecolor="none", zorder=0))

    # ---- inputs (left column) ----
    ix, iw = 0.25, 1.7
    label(ax, ix + iw/2, 4.55, "INPUTS", fontsize=8, fontweight="bold",
          color=C_INPUT)
    box(ax, ix, 3.55, iw, 0.85,
        "Config (JSON)\ndomain, time,\ncell types,\nsubstrates",
        fc="white", ec=C_INPUT, lw=1.3, fontsize=7.0)
    box(ax, ix, 2.55, iw, 0.85,
        "Tumor seeds\n(x, y, ecDNA,\nphase)",
        fc="white", ec=C_INPUT, lw=1.3, fontsize=7.0)
    box(ax, ix, 1.55, iw, 0.85,
        "RNG seed\n(stochasticity)",
        fc="white", ec=C_INPUT, lw=1.3, fontsize=7.0)

    # ---- multi-scale stack (center column) ----
    sx, sw = 2.65, 4.7
    band_w = 0.55  # colored left band carrying the layer name

    def layer(y0, h, band_c, band_label, body_lines, body_fc):
        # main panel
        box(ax, sx, y0, sw, h, "", fc=body_fc, ec=band_c, lw=1.4)
        # colored band on the left
        ax.add_patch(mpatches.Rectangle(
            (sx, y0), band_w, h, facecolor=band_c, edgecolor=band_c,
            linewidth=0, zorder=3))
        label(ax, sx + band_w / 2, y0 + h / 2, band_label, fontsize=8.5,
              fontweight="bold", color="white", rotation=90)
        # body lines (left-aligned past the band, with margin)
        text_x = sx + band_w + 0.15
        n = len(body_lines)
        if n == 1:
            ys = [y0 + h / 2]
        else:
            top_y = y0 + h - 0.28
            bot_y = y0 + 0.28
            ys = np.linspace(top_y, bot_y, n)
        for line, y in zip(body_lines, ys):
            label(ax, text_x, y, line, fontsize=7.8, ha="left")

    layer(3.55, 1.05, C_TISSUE, "TISSUE",
          ["Domain 10 × 30 mm; env grids O$_2$ / glucose / VEGF / lactate",
           "Vascular network, spatial hash (20 µm buckets)"],
          "#E6F0F8")
    layer(2.30, 1.05, C_CELL, "CELL",
          ["9 cell types (tumor, vasculature, immune, normal CNS);",
           "state: type, position, cycle phase, env sample, EGFR, parent_id"],
          "#E6F5F0")
    layer(1.05, 1.05, C_MOL, "MOLECULAR",
          ["ecDNA segregation Binomial(N, 0.5) at mitosis",
           "EGFR ∝ ecDNA × HIF; modulates VEGF, migration, survival"],
          "#FBEEE2")

    # bidirectional coupling arrows between layers (right side of stack)
    cx = sx + sw + 0.10
    bidir_arrow(ax, (cx, 3.55), (cx, 3.35), color=C_LINK, lw=1.0, mut=8)
    bidir_arrow(ax, (cx, 2.30), (cx, 2.10), color=C_LINK, lw=1.0, mut=8)
    label(ax, cx + 0.35, 3.45, "couple", fontsize=6.5, color=C_LINK,
          ha="left")
    label(ax, cx + 0.35, 2.20, "couple", fontsize=6.5, color=C_LINK,
          ha="left")

    # ---- outputs (right column) ----
    ox = 8.05
    ow = 1.7
    label(ax, ox + ow/2, 4.55, "OUTPUTS", fontsize=8, fontweight="bold",
          color=C_OUTPUT)
    box(ax, ox, 3.55, ow, 0.85,
        "cells_t*.tsv\nenvironment_*\nlineage.tsv",
        fc="white", ec=C_OUTPUT, lw=1.3, fontsize=7.0)
    box(ax, ox, 2.55, ow, 0.85,
        "Vega-Lite\nspec + viewer",
        fc="white", ec=C_OUTPUT, lw=1.3, fontsize=7.0)
    box(ax, ox, 1.55, ow, 0.85,
        "report.html\nfigures/",
        fc="white", ec=C_OUTPUT, lw=1.3, fontsize=7.0)

    # input -> stack
    arrow(ax, (ix + iw + 0.05, 3.00), (sx - 0.05, 3.00), color=C_INPUT, lw=1.4)
    # stack -> output
    arrow(ax, (sx + sw + 0.55, 3.00), (ox - 0.05, 3.00), color=C_OUTPUT, lw=1.4)

    # bottom caption strip
    label(ax, 5.0, 0.45,
          "Bidirectional coupling: cells consume / secrete  ⇋  fields shape "
          "behaviors  ⇋  ecDNA gates molecular state",
          fontsize=7.5, color=C_GRAY, fontweight="normal")


# ===================== PANEL B: per-hour loop =====================
def panel_b(ax):
    panel_setup(ax, (0, 10), (0, 3.6))
    panel_letter(ax, "B")
    ax.set_title("Per-hour simulation loop (6 phases)", loc="left",
                 x=0.04, pad=2)

    ax.add_patch(mpatches.Rectangle((0, 0), 10, 3.6, facecolor=C_LIGHT,
                                    edgecolor="none", zorder=0))

    phases = [
        ("1", "Diffuse",
         "O$_2$, glucose,\nVEGF, lactate;\nThomas solver", C_TISSUE),
        ("2", "Update cells",
         "Sample env,\ndeath, divide,\nmigrate, kill", C_CELL),
        ("3", "Sprout",
         "VEGF-driven\nendothelial\nsprouting", "#AA3377"),
        ("4", "Recruit",
         "Immune cells\nfrom vasculature\n(chemokine ↑)", "#EE6677"),
        ("5", "Lyse",
         "Clear necrotic;\nupdate occupancy\nindex", C_GRAY),
        ("6", "Write",
         "Cells / env TSV,\nlineage record,\nprogress log", "#228833"),
    ]
    bw, bh = 1.46, 1.70
    gap = 0.08
    total_w = 6 * bw + 5 * gap
    x0 = (10 - total_w) / 2
    y0 = 1.05
    centers = []
    for i, (num, name, body, col) in enumerate(phases):
        x = x0 + i * (bw + gap)
        # header strip
        box(ax, x, y0 + bh - 0.42, bw, 0.42, "", fc=col, ec=col, lw=0.8)
        # number badge (circle on header left)
        ax.add_patch(mpatches.Circle((x + 0.21, y0 + bh - 0.21), 0.13,
                                     facecolor="white", edgecolor="white",
                                     linewidth=0, zorder=3))
        label(ax, x + 0.21, y0 + bh - 0.21, num, fontsize=9,
              fontweight="bold", color=col)
        # phase name centered on remaining header space
        label(ax, x + 0.21 + 0.18 + (bw - 0.21 - 0.18) / 2,
              y0 + bh - 0.21, name, fontsize=7.8,
              fontweight="bold", color="white")
        # body
        box(ax, x, y0, bw, bh - 0.42, body, fc=C_PANEL, ec=col, lw=0.8,
            fontsize=7.0)
        centers.append((x + bw / 2, y0))
        # arrow to next
        if i < 5:
            arrow(ax, (x + bw + 0.005, y0 + bh / 2),
                  (x + bw + gap - 0.005, y0 + bh / 2),
                  color=C_LINK, lw=1.4, mut=12)

    # loop-back arrow under the row
    yb = y0 - 0.40
    x_left = x0 + 0.05
    x_right = x0 + total_w - 0.05
    arrow(ax, (x_right, yb), (x_left, yb), color=C_TISSUE, lw=1.5, mut=12)
    # vertical stubs
    ax.plot([x_right, x_right], [y0, yb], color=C_TISSUE, lw=1.2)
    ax.plot([x_left, x_left], [yb, y0], color=C_TISSUE, lw=1.2)
    label(ax, 5.0, yb - 0.20, "loop next sim-hour",
          fontsize=7.5, color=C_TISSUE, fontweight="bold")


# ===================== PANEL C: causal mechanism =====================
def panel_c(ax):
    panel_setup(ax, (0, 10), (0, 5.4))
    panel_letter(ax, "C")
    ax.set_title("Mechanistic causal model (per tumor cell)", loc="left",
                 x=0.04, pad=2)

    ax.add_patch(mpatches.Rectangle((0, 0), 10, 5.4, facecolor=C_LIGHT,
                                    edgecolor="none", zorder=0))

    # central tumor-cell area (rectangle with tinted background)
    cell_x0, cell_y0 = 2.85, 0.40
    cell_w, cell_h = 4.50, 4.60
    cx = cell_x0 + cell_w / 2
    cy = cell_y0 + cell_h / 2
    box(ax, cell_x0, cell_y0, cell_w, cell_h, "",
        fc="#FDF7DA", ec="#CCBB44", lw=1.8)
    label(ax, cell_x0 + 0.15, cell_y0 + cell_h - 0.18, "Tumor cell",
          fontsize=8.5, fontweight="bold", color="#7E7320",
          ha="left", va="top")

    # ---- nodes inside the cell ----
    nodes = {
        "ecDNA":           (cx - 1.45, cy + 1.35, C_MOL),
        "EGFR":            (cx + 0.65, cy + 1.35, C_MOL),
        "HIF":             (cx + 0.65, cy + 0.20, C_GRAY),
        "VEGF\nsecretion": (cx - 1.45, cy - 0.85, "#009988"),
        "migration":       (cx + 0.65, cy - 0.85, C_CELL),
        "division\nrate":  (cx - 1.45, cy - 1.80, "#AA3377"),
        "survival":        (cx + 0.65, cy - 1.80, "#EE6677"),
    }
    node_w, node_h = 1.45, 0.55
    centers = {}
    for name, (nx, ny, col) in nodes.items():
        box(ax, nx - node_w/2, ny - node_h/2, node_w, node_h, name,
            fc=C_PANEL, ec=col, lw=1.5, fontsize=8, fontweight="bold",
            color=col)
        centers[name] = (nx, ny)

    # ---- intra-cell edges ----
    def edge(a, b, *, color=C_INK, lw=1.4, label_text=None, label_offset=(0, 0)):
        ax_, ay_ = centers[a]
        bx_, by_ = centers[b]
        # anchor at edge of node box (approximate)
        if abs(ax_ - bx_) > abs(ay_ - by_):
            sx = ax_ + (node_w/2 if bx_ > ax_ else -node_w/2)
            tx = bx_ + (-node_w/2 if bx_ > ax_ else node_w/2)
            sy, ty = ay_, by_
        else:
            sy = ay_ + (node_h/2 if by_ > ay_ else -node_h/2)
            ty = by_ + (-node_h/2 if by_ > ay_ else node_h/2)
            sx, tx = ax_, bx_
        arrow(ax, (sx, sy), (tx, ty), color=color, lw=lw, mut=10)
        if label_text:
            mx = (sx + tx) / 2 + label_offset[0]
            my = (sy + ty) / 2 + label_offset[1]
            label(ax, mx, my, label_text, fontsize=6.5, color=color,
                  fontweight="bold")

    edge("ecDNA", "EGFR", color=C_MOL, label_text="κ (dosage)", label_offset=(0, 0.16))
    edge("HIF", "EGFR", color=C_GRAY, label_text="HIF upreg", label_offset=(0.32, 0))
    edge("EGFR", "VEGF\nsecretion", color="#009988", label_text="β", label_offset=(0.15, 0.15))
    edge("EGFR", "migration", color=C_CELL, label_text="δ", label_offset=(0.13, 0))
    edge("EGFR", "division\nrate", color="#AA3377", label_text="α", label_offset=(-0.15, 0.10))
    edge("EGFR", "survival", color="#EE6677", label_text="γ", label_offset=(0.10, 0.10))

    # ---- external nodes (left side) ----
    ext_w = 2.20
    ext_x = 0.20

    # mitosis (SIV instrument) -> ecDNA
    box(ax, ext_x, cy + 1.35 - 0.40, ext_w, 0.80,
        "Mitosis event\nBinomial(N, 0.5)\nsegregation",
        fc="#FBEEE2", ec=C_MOL, lw=1.4, fontsize=7.2, fontweight="bold",
        color=C_MOL)
    arrow(ax, (ext_x + ext_w, cy + 1.35),
          (nodes["ecDNA"][0] - node_w/2, cy + 1.35),
          color=C_MOL, lw=1.4)
    label(ax, ext_x + ext_w + 0.25, cy + 1.60, "instrument Z",
          fontsize=6.8, color=C_MOL, fontweight="bold", ha="left")

    # O2 local -> HIF (hypoxia)
    box(ax, ext_x, cy + 0.20 - 0.30, ext_w, 0.60,
        "O$_2$ local\n(env grid sample)",
        fc="#E6F0F8", ec=C_TISSUE, lw=1.4, fontsize=7.5, fontweight="bold",
        color=C_TISSUE)
    arrow(ax, (ext_x + ext_w, cy + 0.20),
          (nodes["HIF"][0] - node_w/2, cy + 0.20),
          color=C_TISSUE, lw=1.4)
    label(ax, ext_x + ext_w + 0.25, cy + 0.42,
          "hypoxia indicator (O$_2$ < 18 mmHg)",
          fontsize=6.5, color=C_TISSUE, fontweight="bold", ha="left")

    # ---- external nodes (right side) ----
    rext_w = 1.85
    rext_x = 7.95

    # VEGF field -> angiogenesis
    box(ax, rext_x, cy + 1.35 - 0.40, rext_w, 0.80,
        "VEGF field (env)\ndiffusion + decay",
        fc="#E6F5F0", ec="#009988", lw=1.4, fontsize=7.2, fontweight="bold",
        color="#009988")
    arrow(ax, (nodes["VEGF\nsecretion"][0] + node_w/2, cy - 0.85),
          (rext_x, cy + 0.95), color="#009988", lw=1.2)

    box(ax, rext_x, cy + 0.20 - 0.30, rext_w, 0.60,
        "Angiogenesis\n(Phase 3 sprouting)",
        fc=C_PANEL, ec="#AA3377", lw=1.4, fontsize=7.5, fontweight="bold",
        color="#AA3377")
    arrow(ax, (rext_x + rext_w / 2, cy + 0.95),
          (rext_x + rext_w / 2, cy + 0.50), color="#AA3377", lw=1.2)

    # division -> daughters
    box(ax, rext_x, cy - 0.85 - 0.30, rext_w, 0.60,
        "Daughter cells\n(lineage.tsv)",
        fc="#FDF7DA", ec="#9C8F2E", lw=1.4, fontsize=7.5, fontweight="bold",
        color="#9C8F2E")
    arrow(ax, (nodes["division\nrate"][0] + node_w/2, cy - 1.80),
          (rext_x, cy - 0.85), color="#9C8F2E", lw=1.2)

    # survival -> apoptosis sink
    box(ax, rext_x, cy - 1.80 - 0.30, rext_w, 0.60,
        "Apoptosis /\nnecrosis sink",
        fc=C_PANEL, ec="#EE6677", lw=1.4, fontsize=7.5, fontweight="bold",
        color="#EE6677")
    arrow(ax, (nodes["survival"][0] + node_w/2, cy - 1.80),
          (rext_x, cy - 1.80), color="#EE6677", lw=1.2)


# ===================== legend =====================
def add_legend(fig):
    legend_h = 0.04  # fraction of figure
    legend_ax = fig.add_axes([0.04, 0.005, 0.92, 0.05])
    legend_ax.set_xlim(0, 10)
    legend_ax.set_ylim(0, 1)
    legend_ax.set_xticks([])
    legend_ax.set_yticks([])
    for spine in legend_ax.spines.values():
        spine.set_visible(False)

    entries = [
        ("Tissue scale", C_TISSUE, "rect"),
        ("Cell scale", C_CELL, "rect"),
        ("Molecular scale", C_MOL, "rect"),
        ("Input", C_INPUT, "rect"),
        ("Output", C_OUTPUT, "rect"),
        ("Ground-truth\ncausal edge", C_INK, "arrow"),
        ("Diffusion /\ncoupling", C_TISSUE, "dashed"),
    ]
    x = 0.2
    for name, col, kind in entries:
        if kind == "rect":
            legend_ax.add_patch(mpatches.Rectangle((x, 0.30), 0.30, 0.40,
                                                   facecolor="white",
                                                   edgecolor=col, linewidth=1.5))
        elif kind == "arrow":
            legend_ax.add_patch(FancyArrowPatch((x, 0.5), (x + 0.30, 0.5),
                                                arrowstyle="-|>",
                                                mutation_scale=10,
                                                color=col, linewidth=1.3))
        else:
            legend_ax.add_patch(FancyArrowPatch((x, 0.5), (x + 0.30, 0.5),
                                                arrowstyle="-|>",
                                                mutation_scale=10,
                                                color=col, linewidth=1.0,
                                                linestyle=":"))
        legend_ax.text(x + 0.38, 0.50, name, fontsize=7.5, va="center",
                       ha="left", color=C_INK)
        x += 1.55


# ===================== main =====================
def main():
    fig = plt.figure(figsize=(8.0, 10.5))

    gs = fig.add_gridspec(
        3, 1,
        height_ratios=[5.0, 3.6, 5.4],
        hspace=0.22,
        top=0.97, bottom=0.06, left=0.04, right=0.96,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[2, 0])

    panel_a(ax_a)
    panel_b(ax_b)
    panel_c(ax_c)
    add_legend(fig)

    for out_dir in [Path("figures"), Path("docs/manuscript/figures")]:
        out_dir.mkdir(parents=True, exist_ok=True)
        png = out_dir / "figure8_simulation_architecture.png"
        pdf = out_dir / "figure8_simulation_architecture.pdf"
        fig.savefig(png, dpi=600, bbox_inches="tight",
                    pad_inches=0.05, facecolor="white")
        fig.savefig(pdf, bbox_inches="tight", pad_inches=0.05,
                    facecolor="white")
        print(f"wrote {png}  ({png.stat().st_size/1024:.0f} KB)")
        print(f"wrote {pdf}  ({pdf.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    main()
