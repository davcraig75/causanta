"""Static export of the multi-scale simulation explainer (v21).

A 3-panel companion to ``docs/simulation_multiscale_explorer.html`` for
copy-paste into slides / the manuscript: Section (whole 6 mm tumor) → Tissue
(~400 µm proliferating-margin patch with all nine cell types) → Cell (single
tumor cell with the ecDNA → EGFR → {VEGF, migration, division, apoptosis}
causal chain instantiated). Colours, equations, and zone semantics match the
HTML explainer exactly.

Writes docs/manuscript/figures/figure_simulation_multiscale.{png,pdf}.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.family": "serif", "font.size": 10,
                     "savefig.dpi": 300, "savefig.bbox": "tight"})

ROOT = Path(__file__).resolve().parent.parent

# Cell-type palette + canonical 2 mm baseline final composition (v19 canonical).
COL = {"Neuron": "#4477AA", "Astrocyte": "#66CCEE", "Oligodendrocyte": "#228833",
       "Microglia": "#EE6677", "Endothelial": "#AA3377", "Pericyte": "#BBBBBB",
       "Tumor": "#CCBB44", "RecruitedImmune": "#EE8866", "Necrotic": "#555555"}
ORDER = ["Neuron", "Astrocyte", "Oligodendrocyte", "Microglia", "Endothelial",
         "Pericyte", "Tumor", "RecruitedImmune", "Necrotic"]
ROLE = {"Neuron": "parenchyma (post-mitotic)", "Astrocyte": "parenchyma (reactive glia)",
        "Oligodendrocyte": "parenchyma (myelinating)", "Microglia": "resident immune",
        "Endothelial": "vessel wall", "Pericyte": "vessel support (BBB)",
        "Tumor": "ecDNA-bearing tumour", "RecruitedImmune": "infiltrating T/NK cell",
        "Necrotic": "dead (hypoxic core)"}
FINAL = {"Tumor": 11259, "Endothelial": 5297, "Pericyte": 1270, "Neuron": 565,
         "Necrotic": 439, "RecruitedImmune": 393, "Astrocyte": 311,
         "Oligodendrocyte": 152, "Microglia": 89}

# Simulator structural constants (must match the explainer + the simulator).
SC = dict(X_BASE=2.89, KAPPA=1.21, KHYP=1.5, T_BASE=24.0, VEGF_BASE=600.0,
          MIGR_BASE=10.0, alpha=0.30, beta=0.10, delta=0.05, gamma=0.50,
          A_BASE=5e-5)
T_REPR = 240  # representative time for the snapshot (matches the explainer's default)


def cell_chain(ecdna: int, hypoxic: bool) -> dict:
    M = 1 if hypoxic else 0
    egfr = (SC["X_BASE"] + SC["KAPPA"] * ecdna) * (1 + SC["KHYP"] * M)
    log_ = math.log2(1 + egfr)
    return {
        "EGFR": egfr,
        "T_div": SC["T_BASE"] / (1 + SC["alpha"] * log_),
        "VEGF": SC["VEGF_BASE"] * (1 + SC["beta"] * math.sqrt(egfr)) * (0.2 + 0.8 * M),
        "v": SC["MIGR_BASE"] * (1 + SC["delta"] * egfr) * (2.0 if M else 1.0),
        "apop": SC["A_BASE"] / (1 + SC["gamma"] * log_),
        "M": M,
    }


# -------- panel A: Section (whole 6 mm tumour) -----------------
def panel_section(ax, t=T_REPR):
    ax.set_aspect("equal")
    ax.set_xlim(0, 520); ax.set_ylim(520, 0); ax.set_axis_off()
    ax.set_title("A · Section  (≈ 6 mm × 6 mm)", loc="left", fontsize=11)

    # parenchymal background dots
    rng = np.random.default_rng(7)
    for _ in range(620):
        x, y = rng.uniform(0, 520), rng.uniform(0, 520)
        tp = rng.choice(["Neuron", "Astrocyte", "Oligodendrocyte", "Microglia"], p=[0.45, 0.25, 0.20, 0.10])
        ax.scatter(x, y, s=2, c=COL[tp], alpha=0.7, linewidths=0)

    cx, cy, Rmax = 260, 260, 180
    R = Rmax * (t / 300) ** 0.8
    rNec = max(0, R - 58)
    # radial O₂ field: red core → blue rim
    field = mpatches.Circle((cx, cy), R + 30, transform=ax.transData)
    # approximate the radial gradient by 3 stacked translucent disks
    for rr, col, a in [(R + 30, "#3498db", 0.10), (R * 0.85, "#e67e22", 0.18),
                       (max(rNec * 1.05, 6), "#c0392b", 0.30)]:
        ax.add_patch(mpatches.Circle((cx, cy), rr, color=col, alpha=a, lw=0))

    # vessels
    rng2 = np.random.default_rng(7)
    for i in range(8):
        a = rng2.uniform(0, 2 * np.pi)
        x0, y0 = cx + math.cos(a) * 230, cy + math.sin(a) * 230
        pts = [(x0, y0)]
        xx, yy = x0, y0
        for _ in range(7):
            tx, ty = cx + (rng2.random() - 0.5) * 60, cy + (rng2.random() - 0.5) * 60
            xx += (tx - xx) * 0.3 + (rng2.random() - 0.5) * 26
            yy += (ty - yy) * 0.3 + (rng2.random() - 0.5) * 26
            pts.append((xx, yy))
        xs, ys = zip(*pts)
        ax.plot(xs, ys, color=COL["Endothelial"], lw=1.4, alpha=0.7)

    # tumour mass + necrotic core
    rng3 = np.random.default_rng(11)
    for _ in range(900):
        rn = math.sqrt(rng3.random()); a = rng3.uniform(0, 2 * np.pi)
        rr = rn * Rmax
        if rr > R:
            continue
        nec = rr < rNec
        ax.scatter(cx + math.cos(a) * rr, cy + math.sin(a) * rr,
                   s=4 if nec else 5, c=COL["Necrotic"] if nec else COL["Tumor"],
                   alpha=0.7 if nec else 0.9, linewidths=0)
    # immune infiltrate
    imm_frac = max(0, (t - 60) / 240)
    for _ in range(int(160 * imm_frac)):
        rn = 0.7 + rng3.random() * 0.45; a = rng3.uniform(0, 2 * np.pi)
        ax.scatter(cx + math.cos(a) * rn * Rmax, cy + math.sin(a) * rn * Rmax,
                   s=14, c=COL["RecruitedImmune"], edgecolors="#a33", linewidths=0.4)

    # "magnifier" box marking the patch the Tissue panel shows
    mw, mh = 70, 56
    bx, by = cx + 88, cy + 18
    ax.add_patch(mpatches.Rectangle((bx - mw / 2, by - mh / 2), mw, mh,
                                    fill=False, edgecolor="#111", lw=1.5, linestyle="--"))
    ax.text(bx, by - mh / 2 - 6, "neighbourhood →", ha="center", fontsize=8.5,
            color="#111", fontweight="bold")
    ax.text(8, 510, f"t = {t} h · margin zone shown at right",
            fontsize=8.5, color="#444", style="italic")


# -------- panel B: Tissue (proliferating margin, ~400 µm patch) ---
def panel_tissue(ax):
    ax.set_aspect("equal"); ax.set_xlim(0, 520); ax.set_ylim(420, 0); ax.set_axis_off()
    ax.set_title("B · Tissue  (≈ 400 µm patch · proliferating margin)", loc="left", fontsize=11)

    # green oxygenated tint
    ax.add_patch(mpatches.Rectangle((0, 0), 520, 420, color="#2ecc71", alpha=0.12, lw=0))

    # a vessel through the panel
    xs = np.linspace(20, 500, 40)
    ys = 260 + 40 * np.sin((xs - 20) / 95)
    ax.plot(xs, ys, color=COL["Endothelial"], lw=8, alpha=0.5,
            solid_capstyle="round", zorder=2)
    # local oxygenated halo around the vessel
    for rr, a in [(180, 0.10), (130, 0.10)]:
        ax.add_patch(mpatches.Circle((250, 280), rr, color="#3498db", alpha=a, lw=0))

    rng = np.random.default_rng(((480 + 18) * 131 + 6) & 0xFFFFFFFF)
    # population (matches the explainer's margin-zone counts)
    counts = [("Tumor", 34, 9), ("Endothelial", 9, 6), ("Pericyte", 7, 6),
              ("RecruitedImmune", 9, 13), ("Microglia", 3, 12),
              ("Astrocyte", 4, 9), ("Necrotic", 2, 9)]
    tumor_xy, immune_xy = [], []
    for tp, n, rad in counts:
        for _ in range(n):
            if tp in ("Endothelial", "Pericyte"):
                u = rng.random()
                xx, yy = 20 + u * 480, 260 + 40 * math.sin((xx0 := 20 + u * 480 - 20) / 95) + (rng.random() - 0.5) * 22
            else:
                xx = 30 + rng.random() * 460
                yy = 30 + rng.random() * 360
            ec = max(2, int(28 + (rng.random() - 0.5) * 22)) if tp == "Tumor" else 0
            (ax.add_patch(mpatches.Circle((xx, yy), rad,
                                          facecolor=COL[tp], edgecolor="#222",
                                          linewidth=0.6, alpha=0.7 if tp == "Necrotic" else 0.92, zorder=3)))
            if tp == "Tumor":
                tumor_xy.append((xx, yy, ec))
                for j in range(min(6, ec // 6)):
                    aa = j / max(1, min(6, ec // 6)) * 2 * np.pi
                    ax.scatter(xx + math.cos(aa) * 3.3, yy + math.sin(aa) * 3.3,
                               s=3, c="#7a5c00", zorder=4, linewidths=0)
            elif tp == "RecruitedImmune":
                immune_xy.append((xx, yy))

    # immune→tumour synapses (dashed kill lines)
    for ix, iy in immune_xy:
        best, bd = None, 140
        for tx, ty, ec in tumor_xy:
            d = math.hypot(tx - ix, ty - iy)
            if d < bd:
                bd = d; best = (tx, ty)
        if best:
            ax.plot([ix, best[0]], [iy, best[1]], color="#c0392b", lw=1.4,
                    linestyle=(0, (3, 2)), alpha=0.85, zorder=5)
            ax.add_patch(mpatches.Circle(best, 11, facecolor="none", edgecolor="#c0392b",
                                          lw=2.0, alpha=0.85, zorder=6))

    # legend swatches for the 9 types (tight 2-row grid, panel bottom)
    # y is inverted (set_ylim 420→0): large y = bottom; place legend at y ≈ 395–415.
    ax.text(8, 400, "9 cell types →", fontsize=8.5, color="#444",
            fontweight="bold", va="center")
    for i, tp in enumerate(ORDER):
        col_idx = i % 5; row_idx = i // 5
        cx = 90 + col_idx * 88; cy = 397 + row_idx * 13
        ax.add_patch(mpatches.Circle((cx, cy), 3.5, color=COL[tp], lw=0))
        ax.text(cx + 6, cy + 0.5, tp, fontsize=7.5, va="center", color="#222")


# -------- panel C: Cell (single tumour cell, ecDNA→EGFR→effects) ---
def panel_cell(ax, ecdna=24, hypoxic=True):
    ax.set_aspect("equal"); ax.set_xlim(0, 520); ax.set_ylim(420, 0); ax.set_axis_off()
    ax.set_title("C · Cell  (single tumour cell)", loc="left", fontsize=11)
    out = cell_chain(ecdna, hypoxic)

    # cell membrane + nucleus
    ax.add_patch(mpatches.Circle((130, 210), 96, facecolor="#fdf6dd",
                                  edgecolor=COL["Tumor"], lw=3))
    ax.add_patch(mpatches.Circle((130, 210), 52, facecolor="#efe2a6",
                                  edgecolor="#9a8420", lw=1.8))
    ax.text(130, 110, f"tumour cell · ecDNA = {ecdna}",
            ha="center", fontsize=10, color="#7a5c00", fontweight="bold")
    # ecDNA circles in the nucleus
    rng = np.random.default_rng(ecdna * 97 + (1 if hypoxic else 0))
    for _ in range(ecdna):
        a = rng.uniform(0, 2 * np.pi); r = rng.uniform(0, 44)
        ax.add_patch(mpatches.Circle((130 + math.cos(a) * r, 210 + math.sin(a) * r),
                                      2.4, facecolor="none", edgecolor="#7a5c00", lw=1.2))
    # EGFR line
    ax.text(130, 322,
            f"EGFR = (2.89 + 1.21·{ecdna}){' × (1+1.5)' if hypoxic else ''} = {out['EGFR']:.0f}",
            ha="center", fontsize=9.5, color="#222")
    ax.text(130, 338,
            "hypoxic: HIF-2α boosts EGFR 2.5× (confounder κ_hyp)" if hypoxic else "normoxic",
            ha="center", fontsize=8.5,
            color="#c0392b" if hypoxic else "#2980b9", style="italic")

    # four effect gauges on the right
    gauges = [
        ("α  division",  f"T = {out['T_div']:.1f} h",
            1 - out["T_div"] / SC["T_BASE"],                  "#6a1b9a"),
        ("β  VEGF",      f"{out['VEGF']:.0f} amol/h",
            min(1, out["VEGF"] / 3000),                       "#2e7d32"),
        ("δ  migration", f"{out['v']:.1f} µm/h",
            min(1, out["v"] / 40),                            "#00838f"),
        ("γ  survival",  f"apop ↓ ×{1 / (1 + SC['gamma'] * math.log2(1 + out['EGFR'])):.2f}",
            1 - 1 / (1 + SC["gamma"] * math.log2(1 + out["EGFR"])), "#c0392b"),
    ]
    for i, (label, val, frac, col) in enumerate(gauges):
        y = 70 + i * 80; x = 280
        ax.text(x, y - 6, label, fontsize=10.5, fontweight="bold", color="#222")
        ax.add_patch(mpatches.FancyBboxPatch((x, y), 210, 16,
                                              boxstyle="round,pad=0,rounding_size=8",
                                              facecolor="#eee", edgecolor="none"))
        ax.add_patch(mpatches.FancyBboxPatch((x, y), 210 * max(0.02, min(1, frac)), 16,
                                              boxstyle="round,pad=0,rounding_size=8",
                                              facecolor=col, edgecolor="none"))
        ax.text(x + 205, y + 12, val, ha="right", fontsize=9, color="#fff",
                fontweight="bold")
        ax.annotate("", xy=(x - 6, y + 8), xytext=(226, 210),
                    arrowprops=dict(arrowstyle="->", color="#bbb", lw=0.9))
    # subtitle
    ax.text(8, 400,
            "ecDNA → EGFR (gene dosage + HIF-2α) → {α division, β VEGF, δ migration, γ apoptosis}",
            fontsize=8.5, style="italic", color="#444")


# -------- compose ---------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6),
                          gridspec_kw={"width_ratios": [1.05, 1.0, 1.05]})
panel_section(axes[0], t=T_REPR)
panel_tissue(axes[1])
panel_cell(axes[2], ecdna=24, hypoxic=True)
fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.04, wspace=0.04)

outdir = ROOT / "docs" / "manuscript" / "figures"
for ext in ("png", "pdf"):
    fig.savefig(outdir / f"figure_simulation_multiscale.{ext}")
print(f"saved figure_simulation_multiscale.{{png,pdf}} (t={T_REPR}h, margin zone, ecDNA=24)")
