"""Regenerate Figure 6 (comprehensive multi-seed robustness) from
output/multiseed_full_results.json.

(A) OLS-vs-IV bias % as a function of HIF confounder strength (kappa_hyp =
    1.5 baseline, 0.5 reduced, 0.0 removed), separately for 2 mm and 6 mm,
    with seed-to-seed SD error bars. Bias collapses toward zero as the
    confounding edge is removed.
(B) OLS beta (red) and IV beta (blue) per design cell with seed-SD error
    bars; the OLS-IV gap narrows as HIF decreases.

This generator replaces the one lost in the repository history purge; the
figure is embedded in docs/manuscript as Figure 6.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import warnings

warnings.simplefilter("ignore")
plt.rcParams.update({"font.family": "serif", "font.size": 10, "axes.titlesize": 11,
                     "savefig.dpi": 300, "savefig.bbox": "tight"})

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "manuscript" / "figures" / "comprehensive_robustness_figure.png"

HIF = {"baseline": 1.5, "reduced": 0.5, "removed": 0.0}
SCEN = ["baseline", "reduced", "removed"]
OLS_C, IV_C = "#CC3311", "#0077BB"


def main() -> int:
    d = json.load(open(ROOT / "output" / "multiseed_full_results.json"))
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.4))

    # Panel A: bias % vs HIF strength, one line per scale
    for scale, marker, ls in [("2mm", "o", "-"), ("6mm", "s", "--")]:
        xs, ys, es = [], [], []
        for sc in SCEN:
            cell = d[f"{scale}_{sc}"]
            ps = cell["per_seed"]
            bias = np.array([p["ols_bias_pct"] for p in ps])
            xs.append(HIF[sc]); ys.append(bias.mean()); es.append(bias.std(ddof=1))
        axA.errorbar(xs, ys, yerr=es, marker=marker, ls=ls, capsize=4, lw=1.8,
                     label=f"{scale}", color="#333333" if scale == "2mm" else "#888888",
                     markerfacecolor=OLS_C, markersize=7)
    axA.axhline(0, color="grey", lw=0.8, ls=":")
    axA.set_xlabel("HIF confounder strength κ$_{hyp}$"); axA.set_ylabel("OLS−IV bias (%)")
    axA.set_xticks([0.0, 0.5, 1.5]); axA.set_xticklabels(["0.0\n(removed)", "0.5\n(reduced)", "1.5\n(baseline)"])
    axA.set_title("A. OLS bias collapses when confounding removed"); axA.legend(title="scale", fontsize=8)
    axA.grid(alpha=0.2)

    # Panel B: OLS vs IV beta per cell
    labels = [f"{scl}\n{sc}" for scl in ("2mm", "6mm") for sc in SCEN]
    x = np.arange(len(labels)); w = 0.38
    ols_means, ols_sds, iv_means, iv_sds = [], [], [], []
    for scl in ("2mm", "6mm"):
        for sc in SCEN:
            a = d[f"{scl}_{sc}"]["aggregate"]
            ols_means.append(a["ols_beta_mean"]); ols_sds.append(a["ols_beta_sd"])
            iv_means.append(a["iv_beta_mean"]); iv_sds.append(a["iv_beta_sd"])
    axB.bar(x - w/2, ols_means, w, yerr=ols_sds, capsize=3, color=OLS_C, label="OLS β", alpha=0.85)
    axB.bar(x + w/2, iv_means, w, yerr=iv_sds, capsize=3, color=IV_C, label="IV β (2SLS)", alpha=0.85)
    axB.set_xticks(x); axB.set_xticklabels(labels, fontsize=8)
    axB.set_ylabel("Linear VEGF-on-EGFR coefficient (β)")
    axB.set_title("B. OLS−IV gap narrows as HIF decreases"); axB.legend(fontsize=8)
    axB.grid(axis="y", alpha=0.2)

    fig.suptitle("Multi-seed robustness across scales and confounder strengths (n = 5 seeds/cell, 30 simulations)",
                 fontsize=11, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT); fig.savefig(OUT.with_suffix(".pdf"))
    plt.close(fig)
    print(f"saved {OUT}")
    for scl in ("2mm", "6mm"):
        for sc in SCEN:
            a = d[f"{scl}_{sc}"]["aggregate"]; ps = d[f"{scl}_{sc}"]["per_seed"]
            bias = np.array([p["ols_bias_pct"] for p in ps])
            print(f"  {scl} {sc}: OLS {a['ols_beta_mean']:.3f}±{a['ols_beta_sd']:.3f} "
                  f"IV {a['iv_beta_mean']:.3f}±{a['iv_beta_sd']:.3f} bias {bias.mean():+.1f}%±{bias.std(ddof=1):.1f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
