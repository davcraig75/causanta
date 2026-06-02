"""Regenerate Figure 7 (time-resolved discovery + OLS bias + hypoxic fraction)
from output/multiseed_timeseries_results.json.

(A) PC discovery F1 vs simulation time. (B) OLS-vs-IV bias % vs time.
(C) Hypoxic fraction (fraction of tumor cells below the 18 mmHg HIF threshold)
vs time. Solid lines + circles = 2 mm runs; dashed + squares = 6 mm. Colors:
red = baseline (HIF 1.5), orange = reduced (0.5), green = removed (0.0). Error
bars = seed-to-seed SD (n = 5 seeds per cell).

Replaces the generator lost in the repository history purge; embedded in
docs/manuscript as Figure 7. Run multiseed_timeseries_analysis.py first
(requires runs executed with --keep-all-snapshots).
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
OUT = ROOT / "docs" / "manuscript" / "figures" / "timeseries_figure.png"

COLOR = {"baseline": "#CC3311", "reduced": "#EE7733", "removed": "#009988"}
STYLE = {"2mm": dict(ls="-", marker="o"), "6mm": dict(ls="--", marker="s")}


def series(cell: dict, metric: str):
    ts = sorted(cell["by_t"].keys(), key=int)
    x = [int(t) for t in ts]
    y = [cell["by_t"][t]["agg"][metric]["mean"] for t in ts]
    e = [cell["by_t"][t]["agg"][metric]["sd"] for t in ts]
    return np.array(x), np.array(y), np.array(e)


def main() -> int:
    d = json.load(open(ROOT / "output" / "multiseed_timeseries_results.json"))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    panels = [("f1", "PC discovery F1", "A. Discovery performance vs time"),
              ("ols_bias_pct", "OLS−IV bias (%)", "B. OLS bias vs time"),
              ("hypoxic_fraction", "Hypoxic fraction of tumor", "C. Hypoxic fraction vs time")]

    for ax, (metric, ylab, title) in zip(axes, panels):
        for scale in ("2mm", "6mm"):
            for sc in ("baseline", "reduced", "removed"):
                key = f"{scale}_{sc}"
                if key not in d or not d[key]["by_t"]:
                    continue
                x, y, e = series(d[key], metric)
                ax.errorbar(x, y, yerr=e, color=COLOR[sc], capsize=2, lw=1.6,
                            markersize=5, **STYLE[scale], alpha=0.9)
        if metric == "ols_bias_pct":
            ax.axhline(0, color="grey", lw=0.8, ls=":")
        ax.set_xlabel("Simulation time (hr)"); ax.set_ylabel(ylab); ax.set_title(title)
        ax.grid(alpha=0.2)

    # shared legend
    handles = [plt.Line2D([], [], color=COLOR[sc], lw=2, label=sc) for sc in ("baseline", "reduced", "removed")]
    handles += [plt.Line2D([], [], color="#555555", **STYLE[s], label=s) for s in ("2mm", "6mm")]
    axes[0].legend(handles=handles, fontsize=7.5, ncol=2, loc="lower right")

    fig.suptitle("Time-resolved IV bias, PC discovery, and hypoxic fraction (n = 5 seeds per cell)",
                 fontsize=11, y=1.02)
    plt.tight_layout()
    fig.savefig(OUT); fig.savefig(OUT.with_suffix(".pdf"))
    plt.close(fig)
    print(f"saved {OUT}")
    for scale in ("2mm", "6mm"):
        for sc in ("baseline", "reduced", "removed"):
            key = f"{scale}_{sc}"
            if key in d and d[key]["by_t"]:
                x, y, _ = series(d[key], "f1")
                xb, yb, _ = series(d[key], "ols_bias_pct")
                print(f"  {key}: F1 {y[0]:.2f}->{y[-1]:.2f}, bias {yb[0]:+.1f}%->{yb[-1]:+.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
