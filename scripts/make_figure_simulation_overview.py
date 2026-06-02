"""New simulation-overview figure: a data-driven snapshot of what the CAUSANTA
simulation actually produces (v19, full tissue), as opposed to the schematic
architecture diagram (figure8).

Four panels, all computed from the canonical 2 mm baseline run
(output/large_baseline):
  A. Spatial tissue at the final timestep — every cell coloured by type,
     showing the populated parenchyma + vascular network + tumor mass +
     infiltrating immune + necrotic core.
  B. Oxygen microenvironment field with the hypoxic boundary and tumor outline.
  C. Population dynamics over time (tumor / necrotic / immune).
  D. The somatic instrument: ecDNA copy number -> EGFR expression in tumor
     cells, coloured by hypoxia (the two-cohort first-stage signature).

Writes docs/manuscript/figures/figure_simulation_overview.{png,pdf}.
"""
from __future__ import annotations

import glob
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy import stats
import warnings

warnings.simplefilter("ignore")
plt.rcParams.update({"font.family": "serif", "font.size": 10, "axes.titlesize": 12,
                     "savefig.dpi": 220, "savefig.bbox": "tight"})

ROOT = Path(__file__).resolve().parent.parent
RUN = ROOT / "output" / "large_baseline"
OUT = ROOT / "docs" / "manuscript" / "figures" / "figure_simulation_overview.png"

# config colours per cell type
COLORS = {
    "Neuron": "#4477AA", "Astrocyte": "#66CCEE", "Oligodendrocyte": "#228833",
    "Microglia": "#EE6677", "Endothelial": "#AA3377", "Pericyte": "#BBBBBB",
    "Tumor": "#CCBB44", "RecruitedImmune": "#EE8866", "Necrotic": "#555555",
}
# draw order: background tissue first, tumor/immune/necrotic on top
ORDER = ["Pericyte", "Endothelial", "Oligodendrocyte", "Astrocyte", "Neuron",
         "Microglia", "Necrotic", "Tumor", "RecruitedImmune"]
SIZE = {"Tumor": 5, "RecruitedImmune": 22, "Necrotic": 8, "Endothelial": 3,
        "Pericyte": 3, "Neuron": 9, "Astrocyte": 7, "Oligodendrocyte": 6, "Microglia": 16}


def main() -> int:
    cf = sorted(glob.glob(f"{RUN}/data/cells_t*.tsv"))[-1]
    ef = sorted(glob.glob(f"{RUN}/data/environment_t*.tsv"))[-1]
    df = pd.read_csv(cf, sep="\t", low_memory=False)
    env = pd.read_csv(ef, sep="\t", low_memory=False)
    summ = pd.read_csv(f"{RUN}/data/summary.log", sep="\t")
    t_final = int(cf.split("_t")[-1][:6])

    fig, axes = plt.subplots(2, 2, figsize=(13, 11.5))
    axA, axB, axC, axD = axes.ravel()

    # ---- Panel A: spatial tissue ----
    for name in ORDER:
        sub = df[df["cell_type_name"] == name]
        if len(sub):
            axA.scatter(sub["x"], sub["y"], s=SIZE[name], c=COLORS[name],
                        edgecolors="none", alpha=0.85, label=f"{name} ({len(sub)})")
    axA.set_xlim(0, 2000); axA.set_ylim(0, 2000); axA.set_aspect("equal")
    axA.set_xlabel("x (µm)"); axA.set_ylabel("y (µm)")
    axA.set_title(f"A. Populated tissue at t = {t_final} h (2 mm baseline)")
    axA.legend(loc="upper left", fontsize=6.6, markerscale=1.6, ncol=2,
               framealpha=0.9, bbox_to_anchor=(1.0, 1.0))

    # ---- Panel B: O2 field + hypoxic boundary + tumor outline ----
    nx = env["x_grid"].max() + 1; ny = env["y_grid"].max() + 1
    O2 = env.pivot(index="y_grid", columns="x_grid", values="O2").to_numpy()
    im = axB.imshow(O2, origin="lower", extent=[0, 2000, 0, 2000], cmap="viridis",
                    vmin=0, vmax=20, aspect="equal")
    # hypoxic threshold contour (18 mmHg)
    xs = np.linspace(0, 2000, O2.shape[1]); ys = np.linspace(0, 2000, O2.shape[0])
    axB.contour(xs, ys, O2, levels=[18.0], colors="white", linewidths=1.0, linestyles="--")
    tu = df[df["cell_type_name"] == "Tumor"]
    axB.scatter(tu["x"], tu["y"], s=1.5, c="#CCBB44", alpha=0.25, edgecolors="none")
    nec = df[df["cell_type_name"] == "Necrotic"]
    axB.scatter(nec["x"], nec["y"], s=3, c="#222222", alpha=0.6, edgecolors="none")
    cb = fig.colorbar(im, ax=axB, fraction=0.046, pad=0.04); cb.set_label("O₂ (mmHg)")
    axB.set_xlabel("x (µm)"); axB.set_ylabel("y (µm)")
    axB.set_title("B. Oxygen field (dashed = 18 mmHg hypoxia threshold)")
    axB.legend(handles=[Line2D([], [], ls="--", color="white", label="hypoxic boundary"),
                        Line2D([], [], marker="o", ls="", color="#CCBB44", label="tumor"),
                        Line2D([], [], marker="o", ls="", color="#222222", label="necrotic")],
               loc="upper right", fontsize=7, framealpha=0.85)

    # ---- Panel C: population dynamics ----
    axC.plot(summ["step"], summ["tumor_cells"], color=COLORS["Tumor"], lw=2, label="tumor")
    axC.plot(summ["step"], summ["necrotic_cells"], color=COLORS["Necrotic"], lw=2, label="necrotic")
    axC.plot(summ["step"], summ["immune_cells"], color=COLORS["RecruitedImmune"], lw=2, label="recruited immune")
    axC.set_yscale("symlog", linthresh=1)
    axC.set_xlabel("time (h)"); axC.set_ylabel("cell count (symlog)")
    axC.set_title("C. Population dynamics")
    axC.legend(loc="upper left", fontsize=9); axC.grid(alpha=0.2)
    axC.annotate(f"tumor → {int(summ['tumor_cells'].iloc[-1]):,}\n"
                 f"necrotic → {int(summ['necrotic_cells'].iloc[-1]):,}\n"
                 f"immune → {int(summ['immune_cells'].iloc[-1]):,}",
                 xy=(0.97, 0.05), xycoords="axes fraction", ha="right", va="bottom",
                 fontsize=8.5, bbox=dict(boxstyle="round", fc="wheat", alpha=0.6))

    # ---- Panel D: the somatic instrument (ecDNA -> EGFR) ----
    t = df[df["cell_type_name"] == "Tumor"].copy()
    Z = pd.to_numeric(t["ecDNA_count"], errors="coerce").to_numpy()
    EG = pd.to_numeric(t["egfr_expression"], errors="coerce").to_numpy()
    M = pd.to_numeric(t["is_hypoxic"], errors="coerce").fillna(0).to_numpy().astype(bool)
    ok = np.isfinite(Z) & np.isfinite(EG); Z, EG, M = Z[ok], EG[ok], M[ok]
    axD.scatter(Z[~M], EG[~M], s=8, alpha=0.5, color="#4477AA", label=f"normoxic (n={ (~M).sum() })")
    axD.scatter(Z[M], EG[M], s=8, alpha=0.25, color="#CC6677", label=f"hypoxic (n={M.sum():,})")
    xl = np.linspace(Z.min(), Z.max(), 50)
    if (~M).sum() > 2:
        s, i, *_ = stats.linregress(Z[~M], EG[~M]); axD.plot(xl, i + s*xl, color="#0077BB", lw=2, label=f"normoxic: {i:.1f}+{s:.2f}·Z")
    if M.sum() > 2:
        s, i, *_ = stats.linregress(Z[M], EG[M]); axD.plot(xl, i + s*xl, color="#AA3377", lw=2, ls="--", label=f"hypoxic: {i:.1f}+{s:.2f}·Z")
    axD.set_xlabel("ecDNA copy number Z (instrument)"); axD.set_ylabel("EGFR expression (exposure)")
    axD.set_title("D. Somatic instrument: ecDNA → EGFR (first stage)")
    axD.legend(loc="upper left", fontsize=7.5)

    fig.suptitle("CAUSANTA simulation overview — canonical 2 mm baseline run (v19, full tissue)",
                 fontsize=13, y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(OUT); fig.savefig(OUT.with_suffix(".pdf"))
    plt.close(fig)
    print(f"saved {OUT}")
    print(f"  composition: {df['cell_type_name'].value_counts().to_dict()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
