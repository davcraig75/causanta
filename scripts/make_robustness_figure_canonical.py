"""Regenerate the docx 'Figure 7' robustness panel (OLS-vs-IV migration +
causal-discovery performance) with CANONICAL numbers, replacing the stale
embedded image8 in Causal.v042826.v14.docx.

Panel A: linear OLS vs 2SLS slope of migration_rate on EGFR for the three
         6 mm scenarios (baseline/reduced/removed), with OLS-bias annotations.
Panel B: PC causal-discovery precision / recall / F1 per scenario, from
         output/multiseed_full_results.json (6 mm cells).
"""
from __future__ import annotations

import glob
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings

warnings.simplefilter("ignore")
plt.rcParams.update({"font.family": "serif", "font.size": 10, "savefig.dpi": 300,
                     "savefig.bbox": "tight"})

ROOT = Path(__file__).resolve().parent.parent
TUMOR = 6
OUT = ROOT / "docs" / "manuscript" / "figures" / "robustness_figure.png"


def mig_ols_iv(run_dir: str):
    f = sorted(glob.glob(f"{run_dir}/data/cells_t*.tsv"))
    df = pd.read_csv(f[-1], sep="\t", low_memory=False)
    df["cell_type"] = pd.to_numeric(df["cell_type"], errors="coerce")
    t = df[df["cell_type"] == TUMOR]
    Z = pd.to_numeric(t["ecDNA_count"], errors="coerce").values
    EG = pd.to_numeric(t["egfr_expression"], errors="coerce").values
    R = pd.to_numeric(t["migration_rate"], errors="coerce").values
    m = np.isfinite(Z) & np.isfinite(EG) & np.isfinite(R)
    Z, EG, R = Z[m], EG[m], R[m]
    A = np.column_stack([np.ones(len(EG)), EG]); ols = np.linalg.lstsq(A, R, rcond=None)[0][1]
    Az = np.column_stack([np.ones(len(Z)), Z]); b = np.linalg.lstsq(Az, EG, rcond=None)[0]
    A2 = np.column_stack([np.ones(len(Z)), Az @ b]); iv = np.linalg.lstsq(A2, R, rcond=None)[0][1]
    return ols, iv


def main() -> int:
    scen = [("xlarge_baseline", "6mm_baseline", "Baseline\n(κ_hyp=1.5)"),
            ("xlarge_reduced", "6mm_reduced", "Reduced\n(κ_hyp=0.5)"),
            ("xlarge_removed", "6mm_removed", "Removed\n(κ_hyp=0.0)")]
    ms = json.load(open(ROOT / "output" / "multiseed_full_results.json"))

    labels, ols_v, iv_v, prec, rec, f1 = [], [], [], [], [], []
    for run, key, lab in scen:
        o, i = mig_ols_iv(str(ROOT / "output" / run))
        labels.append(lab); ols_v.append(o); iv_v.append(i)
        a = ms[key]["aggregate"]
        prec.append(a["precision_mean"]); rec.append(a["recall_mean"]); f1.append(a["f1_score_mean"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    x = np.arange(len(labels)); w = 0.36

    ax = axes[0]
    ax.bar(x - w/2, ols_v, w, color="#CC3311", label="OLS", alpha=0.85)
    ax.bar(x + w/2, iv_v, w, color="#33BBEE", label="IV (2SLS)", alpha=0.9)
    for i in range(len(labels)):
        bias = (ols_v[i] / iv_v[i] - 1) * 100
        ax.annotate(f"{bias:+.1f}% bias", (x[i], max(ols_v[i], iv_v[i]) + 0.03),
                    ha="center", fontsize=9,
                    color="#117733" if abs(bias) < 1 else "#CC3311")
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("Linear slope of migration on EGFR")
    ax.set_title("A. OLS vs IV for migration rate (6 mm)")
    ax.set_ylim(0, max(ols_v + iv_v) * 1.25); ax.legend(loc="upper right", fontsize=8)

    ax = axes[1]
    w3 = 0.26
    ax.bar(x - w3, prec, w3, label="Precision", color="#009988", alpha=0.85)
    ax.bar(x, rec, w3, label="Recall", color="#4455AA", alpha=0.85)
    ax.bar(x + w3, f1, w3, label="F1 Score", color="#EE7733", alpha=0.85)
    ax.axhline(0.7, color="grey", ls="--", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylabel("Score"); ax.set_ylim(0, 1.05)
    ax.set_title("B. PC causal-discovery performance (n=5 seeds)")
    ax.legend(loc="lower right", fontsize=8)

    plt.tight_layout()
    fig.savefig(OUT); fig.savefig(OUT.with_suffix(".pdf"))
    plt.close(fig)
    print("migration OLS/IV:", [(round(o,3), round(i,3)) for o, i in zip(ols_v, iv_v)])
    print("discovery F1:", [round(v, 3) for v in f1])
    print("saved", OUT)

    # Swap into docx image8
    docx = ROOT / "Causal.v042826.v14.docx"
    tmp = docx.with_suffix(".docx.tmp")
    with zipfile.ZipFile(docx, "r") as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for it in zin.infolist():
            data = zin.read(it.filename)
            if it.filename == "word/media/image8.png":
                data = OUT.read_bytes()
                print(f"  swapped word/media/image8.png <- robustness_figure.png ({len(data)} bytes)")
            zout.writestr(it, data)
    tmp.replace(docx)
    print("docx image8 refreshed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
