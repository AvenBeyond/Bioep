"""Round-1 clinical association analysis utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from scipy.stats import chi2_contingency, kruskal


def _normalize_sample_id(sample_id: str) -> str:
    sid = str(sample_id).strip().upper()
    return sid[:16] if len(sid) >= 16 else sid


def _safe_numeric(x: pd.Series) -> pd.Series:
    return pd.to_numeric(x, errors="coerce")


def _find_first_col(cols: list[str], keywords: list[str]) -> str | None:
    low = [c.lower() for c in cols]
    for kw in keywords:
        for i, c in enumerate(low):
            if kw in c:
                return cols[i]
    return None


def evaluate_clinical_association(
    labels: pd.Series,
    clinical_df: pd.DataFrame,
    survival_df: pd.DataFrame,
    endpoint: str,
) -> dict[str, Any]:
    _ = clinical_df
    s = survival_df.copy()
    s["normalized_sample_id"] = s["sample"].astype(str).map(_normalize_sample_id)
    lbl = labels.copy().rename("cluster").to_frame()
    lbl["normalized_sample_id"] = lbl.index

    merged = lbl.merge(s, on="normalized_sample_id", how="inner")
    time_col = f"{endpoint}.time"
    event_col = endpoint
    if time_col not in merged.columns or event_col not in merged.columns:
        return {
            "km_ready": False,
            "logrank_p": np.nan,
            "cox_hr_summary": "missing_endpoint_columns",
            "notes": f"endpoint {endpoint} unavailable",
        }

    merged[time_col] = _safe_numeric(merged[time_col])
    merged[event_col] = _safe_numeric(merged[event_col])
    merged = merged.dropna(subset=[time_col, event_col, "cluster"])
    merged[event_col] = merged[event_col].astype(int)
    merged["cluster"] = merged["cluster"].astype(int)

    if merged.empty or merged["cluster"].nunique() < 2:
        return {
            "km_ready": False,
            "logrank_p": np.nan,
            "cox_hr_summary": "insufficient_groups",
            "notes": "not enough clusters with survival data",
        }

    logrank = multivariate_logrank_test(
        event_durations=merged[time_col],
        groups=merged["cluster"],
        event_observed=merged[event_col],
    )

    cox_summary = "cox_not_run"
    try:
        cox_df = merged[[time_col, event_col, "cluster"]].copy()
        cox_df["cluster"] = cox_df["cluster"].astype("category")
        cph = CoxPHFitter()
        cph.fit(cox_df, duration_col=time_col, event_col=event_col, formula="cluster")
        hr_values = np.exp(cph.params_).values
        cox_summary = f"hr_range=[{hr_values.min():.3f},{hr_values.max():.3f}]"
    except (ValueError, np.linalg.LinAlgError, TypeError):
        cox_summary = "cox_failed_singular_or_low_events"

    return {
        "km_ready": True,
        "logrank_p": float(logrank.p_value),
        "cox_hr_summary": cox_summary,
        "notes": "round1_discovery_only",
    }


def plot_km_curve(
    labels: pd.Series,
    survival_df: pd.DataFrame,
    endpoint: str,
    out_path: Path,
) -> None:
    s = survival_df.copy()
    s["normalized_sample_id"] = s["sample"].astype(str).map(_normalize_sample_id)
    lbl = labels.rename("cluster").to_frame()
    lbl["normalized_sample_id"] = lbl.index
    merged = lbl.merge(s, on="normalized_sample_id", how="inner")

    time_col = f"{endpoint}.time"
    event_col = endpoint
    merged[time_col] = _safe_numeric(merged[time_col])
    merged[event_col] = _safe_numeric(merged[event_col])
    merged = merged.dropna(subset=[time_col, event_col, "cluster"])

    if merged.empty:
        return

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(7, 5))
    for cid in sorted(merged["cluster"].astype(int).unique()):
        part = merged[merged["cluster"].astype(int) == cid]
        kmf.fit(part[time_col], event_observed=part[event_col], label=f"cluster {cid}")
        kmf.plot(ci_show=False)
    plt.title(f"KM Curve ({endpoint})")
    plt.xlabel("Time")
    plt.ylabel("Survival probability")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=220)
    plt.close()


def evaluate_clinical_variables(
    labels: pd.Series,
    clinical_df: pd.DataFrame,
) -> list[dict[str, Any]]:
    c = clinical_df.copy()
    sample_col = "sampleID" if "sampleID" in c.columns else c.columns[0]
    c["normalized_sample_id"] = c[sample_col].astype(str).map(_normalize_sample_id)

    merged = labels.rename("cluster").to_frame()
    merged["normalized_sample_id"] = merged.index
    merged = merged.merge(c, on="normalized_sample_id", how="inner")

    results: list[dict[str, Any]] = []
    cols = list(merged.columns)
    age_col = _find_first_col(cols, ["age_at_initial", "age"]) 
    sex_col = _find_first_col(cols, ["gender", "sex"])
    stage_col = _find_first_col(cols, ["pathologic_stage", "tumor_stage", "stage"]) 
    lauren_col = _find_first_col(cols, ["lauren"])

    # Age: Kruskal-Wallis
    if age_col is not None:
        age = _safe_numeric(merged[age_col])
        tmp = merged.assign(_age=age).dropna(subset=["_age", "cluster"])
        groups = [g["_age"].values for _, g in tmp.groupby("cluster") if len(g) > 0]
        if len(groups) >= 2:
            stat, p = kruskal(*groups)
            results.append({"clinical_variable": "age", "association_p": float(p), "notes": f"kruskal_stat={stat:.3f}"})

    for var_name, col in [("sex", sex_col), ("stage", stage_col), ("lauren", lauren_col)]:
        if col is None:
            continue
        tab = pd.crosstab(merged["cluster"], merged[col].astype(str))
        if tab.shape[0] >= 2 and tab.shape[1] >= 2:
            chi2, p, _, _ = chi2_contingency(tab)
            results.append({"clinical_variable": var_name, "association_p": float(p), "notes": f"chi2={chi2:.3f}"})

    return results


def plot_clinical_association_heatmap(df: pd.DataFrame, out_path: Path) -> None:
    if df.empty:
        return
    pivot = df.pivot_table(index="method", columns="clinical_variable", values="association_p", aggfunc="min")
    if pivot.empty:
        return
    sns.set_theme(style="white")
    plt.figure(figsize=(8, 4))
    sns.heatmap(-np.log10(pivot.clip(lower=1e-10)), annot=True, cmap="YlOrRd", fmt=".2f")
    plt.title("Clinical Association Heatmap (-log10 p)")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=220)
    plt.close()
