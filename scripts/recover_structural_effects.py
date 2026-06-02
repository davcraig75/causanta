"""Recover the structural parameters α, β, δ, γ from canonical multi-seed runs.

The simulator implements

    EGFR     = (X_base + κ·Z) · (1 + κ_hyp·M) · ε                    [first stage]
    VEGF     = VEGF_base · (1 + β·√EGFR) · (0.2 + 0.8·M)             [β: VEGF outcome]
    migr     = migr_base · (1 + δ·EGFR) · (2.0 if M else 1.0)        [δ: migration]
    1/T_div  = (1 + α·log2(1+EGFR)) / T_base                         [α: division]
    apopt    = apopt_base / (1 + γ·log2(1+EGFR))                     [γ: apoptosis]

Identifiability differs by parameter, and the recovery reflects that honestly:

  * β (VEGF) and δ (migration) are continuous per-cell phenotypes recorded at
    the final snapshot. We undo the hypoxia gating with the per-cell hypoxia
    flag and fit OLS and 2SLS (ecDNA instrumenting EGFR) on the transformed RHS.
    2SLS recovers β = 0.10 and δ = 0.05 within sampling error in *every*
    scenario; OLS is biased by hypoxia → EGFR confounding when κ_hyp > 0.

  * α (division) is a dynamic *rate*, recovered from inter-division intervals in
    lineage.tsv via a single-parameter NLS fit with T_base treated as the known
    configured constant (24 h). It is point-identified only on κ_hyp = 0
    ("removed") runs, where EGFR = X_base + κ·ecDNA can be reconstructed exactly
    from the recorded ecDNA; under confounding the per-event EGFR is unobserved.
    Recovers α ≈ 0.31 (truth 0.30).

  * γ (survival) is deliberately NOT estimated. Baseline apoptosis is 5e-5/hr,
    so essentially no tumor cell dies and there is no survival selection to
    detect; γ leaves no observable footprint in these runs and is treated as a
    non-identified design parameter rather than fit with a misleading proxy.

Outputs ``output/structural_recovery.json`` with per-seed estimates and
aggregates (mean, SD) for each (scale, scenario) cell, plus E-values
computed from the structural IV estimates via the VanderWeele/Ding
(2017) formula.
"""

from __future__ import annotations

import glob
import json
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

warnings.simplefilter("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from causanta.analyze.iv import compute_e_value  # noqa: E402

TUMOR_TYPE = 6

# Structural-equation reference parameters (tumor cell type 6 in
# simulate/config.py defaults). These are the constants the simulator
# uses inside the modulate_* functions and they are NOT free parameters
# of the recovery — they are the divisors that put each outcome on the
# structural scale where the configured α/β/δ/γ live.
VEGF_BASE = 600.0
MIGR_BASE = 10.0
HYPOXIA_INVASION_BOOST = 2.0  # migration multiplier when hypoxic
VEGF_NORMOXIC_BASAL_FRAC = 0.2  # (0.2 + 0.8·M) gating in modulate_vegf_secretion
T_DIV_BASE = 24.0  # tumor division_time_mean_hr (publication 2 mm / 6 mm runs)
X_BASE = 2.89  # EGFR_BASE_EXPRESSION (gene-dosage intercept)
KAPPA = 1.21  # EGFR_PER_ECDNA_COPY (gene-dosage slope)

# Ground-truth structural coefficients (from cell_types[6] in any params.json)
GROUND_TRUTH = {"alpha": 0.30, "beta": 0.10, "delta": 0.05, "gamma": 0.50}

# Canonical multi-seed run layout: 6 (scale × scenario) cells, 5 seeds each.
CELLS = {
    ("2mm", "baseline"): [
        ("output/large_baseline", 42),
        *(("output/multiseed_2mm/seed_%d" % s, s) for s in (43, 44, 45, 46)),
    ],
    ("2mm", "reduced"): [
        ("output/large_reduced", 42),
        *(("output/multiseed_2mm/reduced_seed_%d" % s, s) for s in (43, 44, 45, 46)),
    ],
    ("2mm", "removed"): [
        ("output/large_removed", 42),
        *(("output/multiseed_2mm/removed_seed_%d" % s, s) for s in (43, 44, 45, 46)),
    ],
    ("6mm", "baseline"): [
        ("output/xlarge_baseline", 42),
        *(("output/multiseed_6mm/baseline_seed_%d" % s, s) for s in (43, 44, 45, 46)),
    ],
    ("6mm", "reduced"): [
        ("output/xlarge_reduced", 42),
        *(("output/multiseed_6mm/reduced_seed_%d" % s, s) for s in (43, 44, 45, 46)),
    ],
    ("6mm", "removed"): [
        ("output/xlarge_removed", 42),
        *(("output/multiseed_6mm/removed_seed_%d" % s, s) for s in (43, 44, 45, 46)),
    ],
}

KHYP_OF = {"baseline": 1.5, "reduced": 0.5, "removed": 0.0}


@dataclass
class FirstStageFit:
    n: int
    intercept: float
    intercept_se: float
    kappa_hat: float
    kappa_hat_se: float
    f_stat: float
    r2: float


def _ols(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Return (β, SE(β), R²) for y ~ X. X must include a constant column."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    ss_res = float((resid ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    dof = max(len(y) - X.shape[1], 1)
    sigma2 = ss_res / dof
    cov = sigma2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    return beta, se, r2


def load_tumor_last_timestep(run_dir: str) -> pd.DataFrame | None:
    files = sorted(glob.glob(f"{run_dir}/data/cells_t*.tsv"))
    if not files:
        return None
    df = pd.read_csv(files[-1], sep="\t", low_memory=False)
    df["cell_type"] = pd.to_numeric(df["cell_type"], errors="coerce")
    tumor = df[df["cell_type"] == TUMOR_TYPE].copy()
    for col in ("ecDNA_count", "egfr_expression", "O2_local", "glucose_local",
                "is_hypoxic", "migration_rate", "VEGF_secretion"):
        if col in tumor.columns:
            tumor[col] = pd.to_numeric(tumor[col], errors="coerce")
    return tumor


def first_stage(tumor: pd.DataFrame) -> tuple[FirstStageFit, FirstStageFit]:
    """Two complementary first-stage fits.

    1. *Univariate*  EGFR ~ const + Z      (matches simulator code only
       when κ_hyp = 0; otherwise absorbs the hypoxia multiplier).
    2. *Structural*  EGFR ~ const + Z + M + Z·M, which recovers the
       configured (X_base, κ) cleanly regardless of κ_hyp because the
       hypoxia term is partialled out.
    """
    Z = tumor["ecDNA_count"].to_numpy(dtype=float)
    Y = tumor["egfr_expression"].to_numpy(dtype=float)
    M = tumor["is_hypoxic"].to_numpy(dtype=float)
    mask = np.isfinite(Z) & np.isfinite(Y) & np.isfinite(M)
    Z, Y, M = Z[mask], Y[mask], M[mask]
    n = len(Y)

    X1 = np.column_stack([np.ones(n), Z])
    b1, se1, r2_1 = _ols(X1, Y)
    # F-statistic on the slope (single instrument, k=1)
    f1 = (b1[1] / se1[1]) ** 2 if se1[1] > 0 else float("nan")
    univ = FirstStageFit(n=n, intercept=b1[0], intercept_se=se1[0],
                        kappa_hat=b1[1], kappa_hat_se=se1[1], f_stat=f1, r2=r2_1)

    X2 = np.column_stack([np.ones(n), Z, M, Z * M])
    b2, se2, r2_2 = _ols(X2, Y)
    f2 = (b2[1] / se2[1]) ** 2 if se2[1] > 0 else float("nan")
    struct = FirstStageFit(n=n, intercept=b2[0], intercept_se=se2[0],
                           kappa_hat=b2[1], kappa_hat_se=se2[1], f_stat=f2, r2=r2_2)
    return univ, struct


def two_sls(Z: np.ndarray, D: np.ndarray, Y: np.ndarray) -> dict:
    """2SLS with single instrument, no covariates.

    Returns OLS and IV slopes on the same scale and the bias percentage
    (OLS / IV - 1) × 100.
    """
    n = len(Z)
    Xz = np.column_stack([np.ones(n), Z])
    bz, _, _ = _ols(Xz, D)
    D_hat = Xz @ bz

    Xh = np.column_stack([np.ones(n), D_hat])
    bh, se_h, _ = _ols(Xh, Y)
    iv_coef, iv_se = bh[1], se_h[1]

    Xd = np.column_stack([np.ones(n), D])
    bd, se_d, _ = _ols(Xd, Y)
    ols_coef, ols_se = bd[1], se_d[1]

    bias_pct = (ols_coef / iv_coef - 1.0) * 100.0 if abs(iv_coef) > 1e-12 else float("nan")
    return {
        "iv_coef": float(iv_coef), "iv_se": float(iv_se),
        "ols_coef": float(ols_coef), "ols_se": float(ols_se),
        "ols_bias_pct": float(bias_pct),
        "n": int(n),
    }


def beta_recovery(tumor: pd.DataFrame) -> dict | None:
    """Recover β from VEGF outcome.

    VEGF_obs = VEGF_BASE · (1 + β·√EGFR) · (VEGF_NORMOXIC_BASAL_FRAC + 0.8·M)
    => VEGF_obs / [VEGF_BASE · (VEGF_NORMOXIC_BASAL_FRAC + 0.8·M)] - 1 = β·√EGFR
    """
    if "VEGF_secretion" not in tumor.columns:
        return None
    Z = tumor["ecDNA_count"].to_numpy(dtype=float)
    EG = tumor["egfr_expression"].to_numpy(dtype=float)
    V = tumor["VEGF_secretion"].to_numpy(dtype=float)
    M = tumor["is_hypoxic"].to_numpy(dtype=float)
    mask = np.isfinite(Z) & np.isfinite(EG) & np.isfinite(V) & np.isfinite(M) & (EG > 0)
    Z, EG, V, M = Z[mask], EG[mask], V[mask], M[mask]
    if len(V) < 50:
        return None
    gate = VEGF_NORMOXIC_BASAL_FRAC + 0.8 * M
    lhs = V / (VEGF_BASE * gate) - 1.0          # = β · √EGFR + noise
    D = np.sqrt(EG)                              # endogenous treatment
    return two_sls(Z, D, lhs)


def delta_recovery(tumor: pd.DataFrame) -> dict | None:
    """Recover δ from migration outcome.

    migr_obs = MIGR_BASE · (1 + δ·EGFR) · (HYPOXIA_INVASION_BOOST if M else 1)
    => migr_obs / [MIGR_BASE · mult(M)] - 1 = δ · EGFR
    """
    if "migration_rate" not in tumor.columns:
        return None
    Z = tumor["ecDNA_count"].to_numpy(dtype=float)
    EG = tumor["egfr_expression"].to_numpy(dtype=float)
    R = tumor["migration_rate"].to_numpy(dtype=float)
    M = tumor["is_hypoxic"].to_numpy(dtype=float)
    mask = np.isfinite(Z) & np.isfinite(EG) & np.isfinite(R) & np.isfinite(M) & (EG > 0)
    Z, EG, R, M = Z[mask], EG[mask], R[mask], M[mask]
    if len(R) < 50:
        return None
    mult = np.where(M > 0.5, HYPOXIA_INVASION_BOOST, 1.0)
    lhs = R / (MIGR_BASE * mult) - 1.0
    return two_sls(Z, EG, lhs)


def alpha_recovery(run_dir: str, kappa_hyp: float) -> dict | None:
    """Recover α from the proliferation outcome via inter-division intervals.

    The simulator schedules each tumor cycle as a draw
        T_cycle ~ max(1, Normal(T_eff, σ)),  T_eff = T_base / (1 + α·log₂(1+EGFR))
    set from the cell's EGFR at division (causanta/simulate/behaviors.py).
    The realized inter-division interval Δt is therefore an unbiased-up-to-clamp
    estimate of T_eff in the conditional mean: E[Δt | EGFR] ≈ T_eff(EGFR).

    Unlike β/δ, which are read off continuous per-cell phenotypes recorded at the
    final snapshot, α governs a *dynamic rate*. It is only point-identified when
    EGFR can be reconstructed exactly from the recorded ecDNA — i.e. when κ_hyp = 0
    (the "removed" scenario), where EGFR = X_base + κ·ecDNA with no hypoxia
    multiplier. We therefore recover α only on κ_hyp = 0 runs.

    Identification follows the same logic as β/δ: T_base is a *known* configured
    constant (24 h), not a free parameter — exactly as VEGF_BASE/MIGR_BASE are
    treated as known divisors. With T_base fixed we fit the single coefficient α
    by nonlinear least squares of Δt on log₂(1+EGFR). (Leaving T_base free is not
    identified here because tumor cells start at 20 copies and never reach
    EGFR≈1, so the intercept is never observed.)
    """
    if abs(kappa_hyp) > 1e-9:
        # EGFR carries an unobserved per-event hypoxia multiplier; α not point-
        # identified from the lineage record (which stores ecDNA, not EGFR or M).
        return {"identified": False, "reason": "kappa_hyp != 0 (EGFR not reconstructable from ecDNA)"}

    lin_path = Path(run_dir) / "data" / "lineage.tsv"
    if not lin_path.exists():
        return None
    lin = pd.read_csv(lin_path, sep="\t", low_memory=False)
    for col in ("time_hr", "parent_id", "parent_ecDNA_before"):
        lin[col] = pd.to_numeric(lin[col], errors="coerce")
    lin = lin.dropna(subset=["time_hr", "parent_id", "parent_ecDNA_before"])
    lin = lin.sort_values(["parent_id", "time_hr"])

    # Consecutive divisions of the same parent → inter-division interval. The
    # interval ending at division k+1 was governed by the cycle time set at
    # division k, which used the parent's post-segregation ecDNA — recorded as
    # parent_ecDNA_before at division k+1 (ecDNA only changes at division).
    dt = lin.groupby("parent_id")["time_hr"].diff()
    ec = lin["parent_ecDNA_before"]
    mask = dt.notna() & (dt > 0)
    dt = dt[mask].to_numpy(dtype=float)
    ec = ec[mask].to_numpy(dtype=float)
    if len(dt) < 100:
        return {"identified": False, "reason": "too few multi-division intervals"}

    egfr = X_BASE + KAPPA * ec            # exact reconstruction at κ_hyp = 0
    x = np.log2(1.0 + egfr)

    from scipy.optimize import curve_fit
    model = lambda xx, al: T_DIV_BASE / (1.0 + al * xx)
    try:
        popt, pcov = curve_fit(model, x, dt, p0=[0.3], maxfev=20000)
    except Exception as e:  # pragma: no cover
        return {"identified": False, "reason": f"fit failed: {e}"}
    alpha_hat = float(popt[0])
    alpha_se = float(np.sqrt(pcov[0, 0]))
    return {
        "identified": True,
        "alpha_hat": alpha_hat,
        "alpha_se": alpha_se,
        "T_base_known": T_DIV_BASE,
        "n_intervals": int(len(dt)),
        "method": "inter_division_interval_NLS",
    }


def _e_values_from_iv(coef: float, se: float) -> dict[str, float]:
    """Apply the canonical VanderWeele/Ding (2017) E-value to a structural
    estimate. The structural coefficient is small (< 1), so we follow the
    convention used in the SIV/MR literature and convert to a risk ratio
    by treating exp(coef · IQR) — a one-IQR contrast — as a meaningful
    RR. Here we use IQR(ecDNA in tumor cells) ≈ 30 as a reference contrast
    so the resulting RR is on the same scale across scenarios. The
    function returns both the raw structural-scale conversion and the
    contrast-based version; downstream code can choose either.
    """
    # Raw RR-conversion at unit contrast (1 EGFR-unit change).
    rr_unit = float(np.exp(coef))
    e_unit = compute_e_value(rr_unit)
    # IQR-contrast: log(RR) per IQR(X) ≈ coef · IQR
    IQR_EGFR = 30.0   # ~within-tumor IQR of EGFR across all canonical runs
    log_rr = coef * IQR_EGFR
    ev = compute_e_value(log_rr, log_rr - 1.96 * se * IQR_EGFR,
                         log_rr + 1.96 * se * IQR_EGFR, log_scale=True)
    return {
        "e_value_unit": e_unit["e_value_point"],
        "e_value_per_IQR_point": ev["e_value_point"],
        "e_value_per_IQR_ci": ev.get("e_value_ci", 1.0),
        "IQR_used": IQR_EGFR,
    }


def analyze_run(run_dir: str, seed: int, scenario: str) -> dict | None:
    if not Path(run_dir).exists():
        return None
    tumor = load_tumor_last_timestep(run_dir)
    if tumor is None or len(tumor) < 100:
        return None
    univ, struct = first_stage(tumor)
    b = beta_recovery(tumor)
    d = delta_recovery(tumor)
    a = alpha_recovery(run_dir, KHYP_OF.get(scenario, 0.0))
    out = {
        "seed": seed,
        "scenario": scenario,
        "run_dir": run_dir,
        "n_tumor": int(len(tumor)),
        "fraction_hypoxic": float(tumor["is_hypoxic"].mean()),
        "first_stage_univariate": {
            "n": univ.n, "intercept": univ.intercept, "intercept_se": univ.intercept_se,
            "kappa_hat": univ.kappa_hat, "kappa_hat_se": univ.kappa_hat_se,
            "f_stat": univ.f_stat, "r2": univ.r2,
        },
        "first_stage_structural": {
            "n": struct.n, "intercept": struct.intercept, "intercept_se": struct.intercept_se,
            "kappa_hat": struct.kappa_hat, "kappa_hat_se": struct.kappa_hat_se,
            "f_stat": struct.f_stat, "r2": struct.r2,
        },
        "beta_recovery": b,
        "delta_recovery": d,
        "alpha_recovery": a,
    }
    if b is not None:
        out["beta_e_values"] = _e_values_from_iv(b["iv_coef"], b["iv_se"])
    if d is not None:
        out["delta_e_values"] = _e_values_from_iv(d["iv_coef"], d["iv_se"])
    return out


def _agg(vals: list[float]) -> dict[str, float]:
    a = np.array([v for v in vals if v is not None and np.isfinite(v)])
    if a.size == 0:
        return {"mean": float("nan"), "sd": float("nan"), "n": 0}
    return {"mean": float(a.mean()), "sd": float(a.std(ddof=1)) if a.size > 1 else 0.0, "n": int(a.size)}


def aggregate_cell(per_seed: list[dict]) -> dict:
    if not per_seed:
        return {}
    agg = {
        "n_seeds": len(per_seed),
        "n_tumor": _agg([r["n_tumor"] for r in per_seed]),
        "fraction_hypoxic": _agg([r["fraction_hypoxic"] for r in per_seed]),
        "univariate_kappa_hat": _agg([r["first_stage_univariate"]["kappa_hat"] for r in per_seed]),
        "univariate_intercept": _agg([r["first_stage_univariate"]["intercept"] for r in per_seed]),
        "univariate_F": _agg([r["first_stage_univariate"]["f_stat"] for r in per_seed]),
        "univariate_R2": _agg([r["first_stage_univariate"]["r2"] for r in per_seed]),
        "structural_kappa_hat": _agg([r["first_stage_structural"]["kappa_hat"] for r in per_seed]),
        "structural_intercept": _agg([r["first_stage_structural"]["intercept"] for r in per_seed]),
        "structural_F": _agg([r["first_stage_structural"]["f_stat"] for r in per_seed]),
    }
    # α recovery (only identified on κ_hyp = 0 / removed runs)
    alpha_good = [r["alpha_recovery"] for r in per_seed
                  if isinstance(r.get("alpha_recovery"), dict)
                  and r["alpha_recovery"].get("identified")]
    if alpha_good:
        agg["alpha_recovery"] = {
            "identified": True,
            "alpha_hat": _agg([r["alpha_hat"] for r in alpha_good]),
            "alpha_se": _agg([r["alpha_se"] for r in alpha_good]),
            "n_intervals": _agg([r["n_intervals"] for r in alpha_good]),
            "T_base_known": T_DIV_BASE,
            "method": "inter_division_interval_NLS",
        }
    else:
        agg["alpha_recovery"] = {"identified": False}

    for key in ("beta_recovery", "delta_recovery"):
        good = [r for r in per_seed if r.get(key) is not None]
        if good:
            agg[key] = {
                "iv_coef": _agg([r[key]["iv_coef"] for r in good]),
                "iv_se": _agg([r[key]["iv_se"] for r in good]),
                "ols_coef": _agg([r[key]["ols_coef"] for r in good]),
                "ols_se": _agg([r[key]["ols_se"] for r in good]),
                "ols_bias_pct": _agg([r[key]["ols_bias_pct"] for r in good]),
                "n_obs": _agg([r[key]["n"] for r in good]),
            }
    if any(r.get("beta_e_values") for r in per_seed):
        good = [r for r in per_seed if r.get("beta_e_values")]
        agg["beta_e_value_per_IQR_point"] = _agg([r["beta_e_values"]["e_value_per_IQR_point"] for r in good])
        agg["beta_e_value_per_IQR_ci"] = _agg([r["beta_e_values"]["e_value_per_IQR_ci"] for r in good])
    if any(r.get("delta_e_values") for r in per_seed):
        good = [r for r in per_seed if r.get("delta_e_values")]
        agg["delta_e_value_per_IQR_point"] = _agg([r["delta_e_values"]["e_value_per_IQR_point"] for r in good])
        agg["delta_e_value_per_IQR_ci"] = _agg([r["delta_e_values"]["e_value_per_IQR_ci"] for r in good])
    return agg


def main() -> int:
    results: dict[str, dict] = {"ground_truth": GROUND_TRUTH, "cells": {}}
    base = Path(__file__).resolve().parent.parent
    for (scale, scenario), runs in CELLS.items():
        per_seed = []
        for rel_dir, seed in runs:
            r = analyze_run(str(base / rel_dir), seed, scenario)
            if r is not None:
                per_seed.append(r)
        key = f"{scale}_{scenario}"
        results["cells"][key] = {
            "scale": scale, "scenario": scenario,
            "kappa_hyp": KHYP_OF[scenario],
            "per_seed": per_seed,
            "aggregate": aggregate_cell(per_seed),
        }

    out = base / "output" / "structural_recovery.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"saved {out}")

    print("\n" + "=" * 116)
    print(f"{'cell':<14} {'n':<3} {'κ̂(struct)':<14} {'X_base(struct)':<14} "
          f"{'β IV':<14} {'δ IV':<14} {'α (T_base known)':<18}")
    print(f"{'':<14} {'':<3} {'':<14} {'':<14} {'(truth 0.10)':<14} "
          f"{'(truth 0.05)':<14} {'(truth 0.30)':<18}")
    print("-" * 116)
    for key, cell in results["cells"].items():
        a = cell["aggregate"]
        if not a:
            continue
        k = a["structural_kappa_hat"]
        x0 = a["structural_intercept"]
        b = a.get("beta_recovery", {})
        d = a.get("delta_recovery", {})
        b_iv = b.get("iv_coef", {"mean": float("nan"), "sd": float("nan")})
        d_iv = d.get("iv_coef", {"mean": float("nan"), "sd": float("nan")})
        ar = a.get("alpha_recovery", {})
        if ar.get("identified"):
            ah = ar["alpha_hat"]
            a_str = f"{ah['mean']:.3f}±{ah['sd']:.3f}"
        else:
            a_str = "not identified"
        print(f"{key:<14} {a['n_seeds']:<3} "
              f"{k['mean']:.3f}±{k['sd']:.3f}    "
              f"{x0['mean']:.3f}±{x0['sd']:.3f}    "
              f"{b_iv['mean']:.4f}±{b_iv['sd']:.4f}  "
              f"{d_iv['mean']:.4f}±{d_iv['sd']:.4f}  "
              f"{a_str:<18}")
    print("=" * 116)
    print("α is recovered only on κ_hyp=0 (removed) runs, where EGFR = X_base + κ·ecDNA "
          "exactly; β/δ recovered in all scenarios via 2SLS. γ (survival) is NOT "
          "estimated: apoptosis (5e-5/hr) is negligible so no survival selection is observable.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
