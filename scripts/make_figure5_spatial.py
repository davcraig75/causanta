"""Generate figure5_spatial_heterogeneity from the two canonical runs.

Panel A: spatial map of the planted-gradient run, tumor cells colored by the
         zone they fall in (core/margin/infiltrating by distance from the
         tumor-seed centroid).
Panel B: per-zone 2SLS estimate of the migration effect delta for both runs:
         the constant-delta control (flat ~0.05) and the planted-gradient run
         (recovers 0.03 / 0.07 / 0.05), with the planted values marked. This
         shows the framework recovers genuine spatial heterogeneity when it is
         present and returns a flat map when it is not (negative control).

Writes docs/manuscript/figures/figure5_spatial_heterogeneity.{png,pdf}.
"""

from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings

warnings.simplefilter("ignore")

plt.rcParams.update({
    "font.family": "serif", "font.size": 10, "axes.labelsize": 11,
    "axes.titlesize": 11, "legend.fontsize": 8, "figure.dpi": 150,
    "savefig.dpi": 300, "savefig.bbox": "tight",
})

ROOT = Path(__file__).resolve().parent.parent
TUMOR = 6
MIGR_BASE = 10.0
HYP_BOOST = 2.0
RADII = (780.0, 970.0)
PLANTED = {"core": 0.03, "margin": 0.07, "infiltrating": 0.05}
ZONES = ["core", "margin", "infiltrating"]


def _final_tumor(run_dir: Path) -> pd.DataFrame:
    files = sorted(glob.glob(f"{run_dir}/data/cells_t*.tsv"))
    df = pd.read_csv(files[-1], sep="\t", low_memory=False)
    df["cell_type"] = pd.to_numeric(df["cell_type"], errors="coerce")
    t = df[df["cell_type"] == TUMOR].copy()
    for c in ("x", "y", "ecDNA_count", "egfr_expression", "migration_rate", "is_hypoxic"):
        t[c] = pd.to_numeric(t[c], errors="coerce")
    return t.dropna(subset=["x", "y", "ecDNA_count", "egfr_expression", "migration_rate"])


def _centroid(run_dir: Path) -> tuple[float, float]:
    p = json.load(open(f"{run_dir}/data/params.json"))
    s = p.get("tumor_seeds", [])
    return float(np.mean([c["x"] for c in s])), float(np.mean([c["y"] for c in s]))


def _zone(dist: np.ndarray) -> np.ndarray:
    z = np.full(len(dist), "infiltrating", dtype=object)
    z[dist < RADII[1]] = "margin"
    z[dist < RADII[0]] = "core"
    return z


def _delta_2sls(df: pd.DataFrame) -> tuple[float, float, int]:
    Z = df["ecDNA_count"].to_numpy(float); EG = df["egfr_expression"].to_numpy(float)
    R = df["migration_rate"].to_numpy(float); M = df["is_hypoxic"].fillna(0).to_numpy(float)
    if len(Z) < 20:
        return float("nan"), float("nan"), len(Z)
    mult = np.where(M > 0.5, HYP_BOOST, 1.0)
    lhs = R / (MIGR_BASE * mult) - 1.0
    A = np.column_stack([np.ones(len(Z)), Z]); b, *_ = np.linalg.lstsq(A, EG, rcond=None)
    A2 = np.column_stack([np.ones(len(Z)), A @ b]); d, *_ = np.linalg.lstsq(A2, lhs, rcond=None)
    resid = lhs - A2 @ d
    se = float(np.sqrt((resid**2).sum() / max(len(Z) - 2, 1) * np.linalg.inv(A2.T @ A2)[1, 1]))
    return float(d[1]), se, len(Z)


def per_zone(run_dir: Path) -> dict:
    t = _final_tumor(run_dir)
    cx, cy = _centroid(run_dir)
    t["dist"] = np.sqrt((t["x"] - cx) ** 2 + (t["y"] - cy) ** 2)
    t["zone"] = _zone(t["dist"].to_numpy())
    out = {}
    for z in ZONES:
        d, se, n = _delta_2sls(t[t["zone"] == z])
        out[z] = (d, se, n)
    return out, t, (cx, cy)


def main() -> int:
    spatial_dir = ROOT / "output" / "spatial_hetero_2mm"
    const_dir = ROOT / "output" / "large_baseline"

    sp, sp_t, (cx, cy) = per_zone(spatial_dir)
    co, _, _ = per_zone(const_dir)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Panel A: spatial zone map of the planted-gradient run.
    ax = axes[0]
    colors = {"core": "#4477AA", "margin": "#EE6677", "infiltrating": "#228833"}
    for z in ZONES:
        sub = sp_t[sp_t["zone"] == z]
        ax.scatter(sub["x"], sub["y"], s=8, alpha=0.5, color=colors[z],
                   label=f"{z} (δ={PLANTED[z]})")
    # zone circles
    th = np.linspace(0, 2 * np.pi, 200)
    for r in RADII:
        ax.plot(cx + r * np.cos(th), cy + r * np.sin(th), "k--", lw=0.6, alpha=0.5)
    ax.set_aspect("equal")
    ax.set_xlabel("X (μm)"); ax.set_ylabel("Y (μm)")
    ax.set_title("A. Planted spatial zones (2 mm run)")
    ax.legend(loc="upper right", fontsize=7, framealpha=0.9)

    # Panel B: per-zone recovered delta, both runs.
    ax = axes[1]
    x = np.arange(len(ZONES)); w = 0.38
    const_vals = [co[z][0] for z in ZONES]
    const_err = [1.96 * co[z][1] for z in ZONES]
    sp_vals = [sp[z][0] for z in ZONES]
    sp_err = [1.96 * sp[z][1] for z in ZONES]

    ax.bar(x - w/2, const_vals, w, yerr=const_err, capsize=4, color="#BBBBBB",
           label="constant-δ run (control)")
    ax.bar(x + w/2, sp_vals, w, yerr=sp_err, capsize=4, color="#CC6677",
           label="planted-gradient run (2SLS)")
    # planted-truth markers
    ax.scatter(x + w/2, [PLANTED[z] for z in ZONES], marker="_", s=420,
               color="black", zorder=5, label="planted ground truth")
    ax.axhline(0.05, color="#888888", ls=":", lw=0.8)

    ax.set_xticks(x); ax.set_xticklabels([z.capitalize() for z in ZONES])
    ax.set_ylabel("Migration effect δ (per-zone 2SLS)")
    ax.set_title("B. Per-zone 2SLS recovers the planted gradient")
    ax.set_ylim(0, 0.09)
    ax.legend(loc="upper left", fontsize=7.5)

    plt.tight_layout()
    out_png = ROOT / "docs" / "manuscript" / "figures" / "figure5_spatial_heterogeneity.png"
    out_pdf = out_png.with_suffix(".pdf")
    fig.savefig(out_png); fig.savefig(out_pdf)
    plt.close(fig)

    print("Per-zone δ (constant-δ control):", {z: round(co[z][0], 4) for z in ZONES})
    print("Per-zone δ (planted-gradient run):", {z: round(sp[z][0], 4) for z in ZONES})
    print("Planted ground truth:", PLANTED)
    print(f"Saved {out_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
