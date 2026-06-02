"""Structural-parameter recovery & identifiability figure (v20).

Visualises which of the four ecDNA→phenotype effects are recovered from the
canonical multi-seed runs, and which are not identified — the result of the
v20 α estimator + γ reframing.

Left panel  : recovered/configured ratio (truth = 1.0) with 95% CI for the
              point-identified coefficients κ, β, δ, α. γ is shown as a hatched
              "not identified" row.
Right panel : one-line rationale per parameter (identifiability class).

Reads output/structural_recovery.json (κ, β, δ, α) and writes
docs/manuscript/figures/figure_recovery_identifiability.{png,pdf}.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

plt.rcParams.update({"font.family": "serif", "font.size": 10,
                     "savefig.dpi": 300, "savefig.bbox": "tight"})

ROOT = Path(__file__).resolve().parent.parent
SR = json.loads((ROOT / "output" / "structural_recovery.json").read_text())
CELLS = SR["cells"]

TRUTH = {"kappa": 1.21, "beta": 0.10, "delta": 0.05, "alpha": 0.30, "gamma": 0.50}


def mean_sd(vals):
    a = np.array([v for v in vals if v is not None and np.isfinite(v)])
    return (float(a.mean()), float(a.std(ddof=1)) if a.size > 1 else 0.0) if a.size else (np.nan, np.nan)


# Aggregate point-identified estimates across the (scale, scenario) cells.
kappa = mean_sd([c["aggregate"]["structural_kappa_hat"]["mean"] for c in CELLS.values()])
beta = mean_sd([c["aggregate"]["beta_recovery"]["iv_coef"]["mean"]
                for c in CELLS.values() if "beta_recovery" in c["aggregate"]])
delta = mean_sd([c["aggregate"]["delta_recovery"]["iv_coef"]["mean"]
                 for c in CELLS.values() if "delta_recovery" in c["aggregate"]])
# α: only κ_hyp = 0 (removed) cells
alpha_cells = [c["aggregate"]["alpha_recovery"] for c in CELLS.values()
               if c["aggregate"].get("alpha_recovery", {}).get("identified")]
alpha = mean_sd([a["alpha_hat"]["mean"] for a in alpha_cells])

EST = {"kappa": kappa, "beta": beta, "delta": delta, "alpha": alpha}

rows = [
    ("κ", "kappa", "first stage", "ecDNA→EGFR slope, recorded per cell", "yes"),
    ("β", "beta", "2SLS (all)", "VEGF: continuous per-cell readout", "yes"),
    ("δ", "delta", "2SLS (all)", "migration: continuous per-cell readout", "yes"),
    ("α", "alpha", "κ_hyp=0", "division rate: inter-division intervals, T_base known", "yes"),
    ("γ", "gamma", "—", "apoptosis 5e-5/hr → no deaths → no survival selection", "no"),
]

C_YES, C_NO = "#2e7d32", "#c0392b"
fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.4), gridspec_kw={"width_ratios": [1.0, 1.05]})
fig.subplots_adjust(left=0.13, right=0.985, top=0.86, bottom=0.13, wspace=0.28)

# ---- Left: recovered / truth ratio ----
ys = np.arange(len(rows))[::-1]
axL.axvline(1.0, color="#444", lw=1.2, ls="--", zorder=1)
for y, (sym, key, _, _, ident) in zip(ys, rows):
    if ident == "yes":
        m, sd = EST[key]
        t = TRUTH[key]
        ratio, rsd = m / t, sd / t
        axL.errorbar(ratio, y, xerr=1.96 * rsd, fmt="o", color=C_YES, ms=8,
                     capsize=4, lw=2, zorder=3)
        axL.text(ratio, y + 0.18, f"{m:.3g} / {t:.3g}", ha="center", va="bottom",
                 fontsize=8.5, color="#1b5e20")
    else:
        axL.barh(y, 0.62, left=0.69, height=0.5, color="none",
                 edgecolor=C_NO, hatch="////", zorder=2)
        axL.text(1.0, y, "not identified", ha="center", va="center",
                 fontsize=9, color=C_NO, fontweight="bold")
axL.set_yticks(ys)
axL.set_yticklabels([f"{sym}  ({key})" for sym, key, *_ in rows], fontsize=11)
axL.set_xlim(0.55, 1.45)
axL.set_xlabel("recovered / configured  (truth = 1.0)")
axL.set_title("Structural-parameter recovery (5 seeds/cell)", fontsize=10.5, loc="left")
axL.spines[["top", "right"]].set_visible(False)

# ---- Right: identifiability rationale table ----
axR.axis("off")
axR.set_title("Identifiability class", fontsize=11, loc="left")
axR.text(0.0, 1.0, "param", fontweight="bold", fontsize=9, transform=axR.transAxes)
axR.text(0.16, 1.0, "method", fontweight="bold", fontsize=9, transform=axR.transAxes)
axR.text(0.42, 1.0, "why", fontweight="bold", fontsize=9, transform=axR.transAxes)
for i, (sym, key, method, why, ident) in enumerate(rows):
    yy = 0.88 - i * 0.17
    col = C_YES if ident == "yes" else C_NO
    axR.text(0.0, yy, sym, fontsize=13, color=col, fontweight="bold", transform=axR.transAxes)
    axR.text(0.16, yy, method, fontsize=9, transform=axR.transAxes)
    axR.text(0.42, yy, why, fontsize=8.3, transform=axR.transAxes)
axR.text(0.0, -0.04,
         "β, δ are deterministic per-cell phenotypes → read off the snapshot. "
         "α is a dynamic rate, identified only where EGFR is exactly\nreconstructable "
         "from ecDNA (κ_hyp=0). γ leaves no footprint because apoptosis is negligible "
         "— a design parameter, not a fit target.",
         fontsize=7.8, style="italic", color="#555", transform=axR.transAxes)

legend = [Patch(facecolor=C_YES, label="point-identified"),
          Patch(facecolor="none", edgecolor=C_NO, hatch="////", label="not identified")]
axL.legend(handles=legend, loc="upper left", fontsize=8, frameon=False)
for ext in ("png", "pdf"):
    out = ROOT / "docs" / "manuscript" / "figures" / f"figure_recovery_identifiability.{ext}"
    fig.savefig(out)
print("saved figure_recovery_identifiability.{png,pdf}")
print(f"  kappa={kappa[0]:.3f}  beta={beta[0]:.4f}  delta={delta[0]:.4f}  "
      f"alpha={alpha[0]:.3f} (removed)  gamma=not identified")
