"""Regenerate the 3-panel first-stage figure (Causal v15 'Figure 4') from the
seed-43 6 mm baseline run, matching the manuscript caption exactly:

  (A) ecDNA segregation: histogram of daughter fractions (of the replicated
      pool), mean ~0.501, var ~0.0071, t-test vs 0.5.
  (B) First-stage regression coloured by hypoxia: normoxic cohort
      (EGFR = 2.76 + 1.22·Z) and hypoxic cohort (EGFR = 7.25 + 3.04·Z),
      with the combined naive OLS line.
  (C) Instrument variability: distribution of ecDNA copy number per tumor
      cell (median 23, mean 24.3, IQR, 5-95%, range 0-63).

Writes docs/manuscript/figures/figure1_segregation_firststage.png (the file
embedded as image5 in the manuscript docx).
"""
from __future__ import annotations

import glob
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import warnings

warnings.simplefilter("ignore")
plt.rcParams.update({"font.family": "serif", "font.size": 9.5, "axes.titlesize": 11,
                     "savefig.dpi": 300, "savefig.bbox": "tight"})

ROOT = Path(__file__).resolve().parent.parent
RUN = ROOT / "output" / "multiseed_6mm" / "baseline_seed_43"
OUT = ROOT / "docs" / "manuscript" / "figures" / "figure1_segregation_firststage.png"
TUMOR = 6


def main() -> int:
    lin = pd.read_csv(RUN / "data" / "lineage.tsv", sep="\t")
    repl = lin["parent_ecDNA_after"] + lin["daughter_ecDNA"]
    frac = (lin["daughter_ecDNA"] / repl)[repl > 0].to_numpy()
    t, p = stats.ttest_1samp(frac, 0.5)

    df = pd.read_csv(sorted(glob.glob(f"{RUN}/data/cells_t*.tsv"))[-1], sep="\t", low_memory=False)
    df["cell_type"] = pd.to_numeric(df["cell_type"], errors="coerce")
    tu = df[df["cell_type"] == TUMOR]
    Z = pd.to_numeric(tu["ecDNA_count"], errors="coerce").to_numpy()
    EG = pd.to_numeric(tu["egfr_expression"], errors="coerce").to_numpy()
    M = pd.to_numeric(tu["is_hypoxic"], errors="coerce").fillna(0).to_numpy().astype(bool)
    ok = np.isfinite(Z) & np.isfinite(EG)
    Z, EG, M = Z[ok], EG[ok], M[ok]

    sn, in_, rn, _, _ = stats.linregress(Z[~M], EG[~M])
    sh, ih, rh, _, _ = stats.linregress(Z[M], EG[M])
    sc, ic, rc, _, _ = stats.linregress(Z, EG)
    # first-stage F (combined)
    n = len(Z)
    F = (rc**2 / 1) / ((1 - rc**2) / (n - 2))

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 5.0))  # aspect ~2.6 to match docx extent

    # Panel A
    ax = axes[0]
    ax.hist(frac, bins=30, density=True, color="#9ecae1", edgecolor="#3182bd", alpha=0.85)
    ax.axvline(0.5, color="red", ls="--", lw=2, label="Expected mean = 0.5")
    ax.axvline(frac.mean(), color="#08519c", lw=2, label=f"Observed mean = {frac.mean():.3f}")
    ax.text(0.03, 0.97, f"n = {len(frac):,}\nvar = {frac.var():.4f}\nt = {t:.2f}, p = {p:.2f}",
            transform=ax.transAxes, va="top", fontsize=8.5,
            bbox=dict(boxstyle="round", fc="wheat", alpha=0.6))
    ax.set_xlabel("Daughter ecDNA fraction"); ax.set_ylabel("Density")
    ax.set_title("A. ecDNA segregation distribution"); ax.legend(loc="upper right", fontsize=7.5); ax.set_xlim(0, 1)

    # Panel B
    ax = axes[1]
    xl = np.linspace(Z.min(), Z.max(), 100)
    ax.scatter(Z[~M], EG[~M], s=9, alpha=0.5, color="#4477AA", label=f"normoxic (n = {(~M).sum()})")
    ax.scatter(Z[M], EG[M], s=9, alpha=0.25, color="#CC6677", label=f"hypoxic (n = {M.sum():,})")
    ax.plot(xl, in_ + sn * xl, color="#0077BB", lw=2.2, label=f"normoxic: EGFR = {in_:.2f} + {sn:.2f}·Z")
    ax.plot(xl, ih + sh * xl, color="#AA3377", lw=2.2, ls="--", label=f"hypoxic: EGFR = {ih:.2f} + {sh:.2f}·Z")
    ax.plot(xl, ic + sc * xl, color="black", lw=1.4, ls=":", label=f"combined OLS (β = {sc:.2f})")
    ax.text(0.03, 0.97, f"κ̂(normoxic) = {sn:.2f}\nF = {F:,.0f}", transform=ax.transAxes, va="top",
            fontsize=8.5, bbox=dict(boxstyle="round", fc="wheat", alpha=0.6))
    ax.set_xlabel("ecDNA copy number (Z)"); ax.set_ylabel("EGFR expression")
    ax.set_title("B. First-stage regression"); ax.legend(loc="lower right", fontsize=6.8)

    # Panel C
    ax = axes[2]
    ax.hist(Z, bins=range(0, int(Z.max()) + 2), color="#a1d99b", edgecolor="#31a354", alpha=0.85)
    q1, med, q3 = np.percentile(Z, [25, 50, 75])
    p5, p95 = np.percentile(Z, [5, 95])
    ax.axvline(med, color="black", lw=1.8, label=f"median = {med:.0f}")
    ax.axvspan(q1, q3, color="grey", alpha=0.18, label=f"IQR {q1:.0f}–{q3:.0f}")
    ax.text(0.97, 0.97, f"mean = {Z.mean():.1f}\n5–95%: {p5:.0f}–{p95:.0f}\nrange {Z.min():.0f}–{Z.max():.0f}",
            transform=ax.transAxes, va="top", ha="right", fontsize=8.5,
            bbox=dict(boxstyle="round", fc="wheat", alpha=0.6))
    ax.set_xlabel("ecDNA copy number per tumor cell"); ax.set_ylabel("Cell count")
    ax.set_title("C. Instrument variability"); ax.legend(loc="center right", fontsize=7.5)

    plt.tight_layout()
    fig.savefig(OUT); fig.savefig(OUT.with_suffix(".pdf"))
    plt.close(fig)
    print(f"saved {OUT}")
    print(f"A: n_div={len(frac):,} mean={frac.mean():.4f} var={frac.var():.4f} t={t:.2f} p={p:.2f}")
    print(f"B: tumor n={n:,} normoxic={(~M).sum()} hypoxic={M.sum():,}; "
          f"normoxic {in_:.2f}+{sn:.2f}Z, hypoxic {ih:.2f}+{sh:.2f}Z, combined β={sc:.2f}, F={F:,.0f}")
    print(f"C: median={med:.0f} mean={Z.mean():.1f} IQR {q1:.0f}-{q3:.0f} range {Z.min():.0f}-{Z.max():.0f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
