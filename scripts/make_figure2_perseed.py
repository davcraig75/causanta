"""Regenerate figure2_ols_vs_iv.png to MATCH the manuscript Figure 2 caption:
a per-seed comparison of OLS (red) and IV (blue) estimates of the EGFR -> VEGF
coefficient at the two BASELINE cells (2 mm and 6 mm, HIF = 1.5). Each of the
n = 5 seeds is a marker; a dashed line marks the per-estimator mean and a
shaded band spans mean +/- 1 SD. The OLS band is wider than the IV band,
illustrating IV's variance reduction under confounding.

(The previous figure2 was a baseline/reduced/removed bar chart, which is the
content of Figure 6, not the per-seed baseline comparison Figure 2 describes.)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.simplefilter("ignore")
plt.rcParams.update({"font.family": "serif", "font.size": 10, "axes.labelsize": 11,
                     "axes.titlesize": 11, "savefig.dpi": 300, "savefig.bbox": "tight"})

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "manuscript" / "figures" / "figure2_ols_vs_iv.png"

OLS_C, IV_C = "#CC3311", "#0077BB"


def main() -> int:
    d = json.load(open(ROOT / "output" / "multiseed_full_results.json"))
    cells = [("2mm_baseline", "2 mm baseline (HIF = 1.5)"),
             ("6mm_baseline", "6 mm baseline (HIF = 1.5)")]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=False)
    rng = np.random.default_rng(0)

    for ax, (key, title) in zip(axes, cells):
        ps = d[key]["per_seed"]; a = d[key]["aggregate"]
        ols = np.array([p["ols_beta"] for p in ps])
        iv = np.array([p["iv_beta"] for p in ps])
        rows = [("OLS", ols, OLS_C, a["ols_beta_mean"], a["ols_beta_sd"], 1.0),
                ("IV (2SLS)", iv, IV_C, a["iv_beta_mean"], a["iv_beta_sd"], 0.0)]
        for label, vals, color, mean, sd, y in rows:
            # shaded band = mean +/- 1 SD (vertical), dashed line at mean
            ax.axvspan(mean - sd, mean + sd, ymin=0.08 + (0.0 if y == 0 else 0.5),
                       ymax=0.42 + (0.0 if y == 0 else 0.5), color=color, alpha=0.15)
            ax.vlines(mean, y - 0.32, y + 0.32, color=color, ls="--", lw=1.6)
            jit = (rng.random(len(vals)) - 0.5) * 0.28
            ax.scatter(vals, np.full(len(vals), y) + jit, color=color, s=55,
                       edgecolor="white", linewidth=0.6, zorder=5,
                       label=f"{label}: {mean:.3f} ± {sd:.3f}")
        bias = a["ols_bias_pct_mean"]
        ax.set_yticks([0, 1]); ax.set_yticklabels(["IV (2SLS)", "OLS"])
        ax.set_ylim(-0.6, 1.6)
        ax.set_xlabel("Linear slope of VEGF on EGFR (β, a.u.)")
        ax.set_title(f"{title}\nOLS overshoots IV by {bias:+.1f}%")
        ax.legend(loc="lower right", fontsize=7.5, framealpha=0.9)
        ax.grid(axis="x", alpha=0.25)

    fig.suptitle("Per-seed OLS vs IV estimates of the EGFR → VEGF coefficient (n = 5 seeds per cell)",
                 fontsize=11, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT); fig.savefig(OUT.with_suffix(".pdf"))
    plt.close(fig)
    print("saved", OUT)
    for key, _ in cells:
        a = d[key]["aggregate"]
        print(f"  {key}: OLS {a['ols_beta_mean']:.3f}±{a['ols_beta_sd']:.3f}  "
              f"IV {a['iv_beta_mean']:.3f}±{a['iv_beta_sd']:.3f}  bias {a['ols_bias_pct_mean']:+.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
