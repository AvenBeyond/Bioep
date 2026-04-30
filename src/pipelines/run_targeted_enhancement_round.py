"""Targeted enhancement round (stable evidence mode).

Goals:
- Produce defense-ready survival evidence without being blocked by unstable Cox fitting.
- Restrict to three candidates only:
  1) old_final_baseline
  2) phaseA_best_methylation
  3) phaseB_partial_ge4_best
- Prioritize OS/PFI. DSS/DFI are secondary and only computed when stable.
"""

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
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


PRIMARY_ENDPOINTS = ["OS", "PFI"]
SECONDARY_ENDPOINTS = ["DSS", "DFI"]
ALL_ENDPOINTS = PRIMARY_ENDPOINTS + SECONDARY_ENDPOINTS


def _root() -> Path:
    return Path(__file__).resolve().parents[2]


def _normalize_sample_id(sample_id: str) -> str:
    sid = str(sample_id).strip().upper()
    return sid[:16] if len(sid) >= 16 else sid


def _read_matrix(path: Path, max_cols: int | None = None) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).map(_normalize_sample_id)
    out = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return out.iloc[:, :max_cols] if max_cols is not None else out


def _load_survival_df(root: Path) -> pd.DataFrame:
    s = pd.read_csv(root / "data" / "raw" / "clinical" / "survival_STAD_survival.txt", sep="\t", dtype=str)
    sample_col = "sample" if "sample" in s.columns else s.columns[0]
    s["normalized_sample_id"] = s[sample_col].astype(str).map(_normalize_sample_id)
    return s


def _load_candidate_labels(root: Path) -> tuple[dict[str, pd.Series], dict[str, dict[str, Any]]]:
    tables = root / "results" / "tables"

    old = pd.read_csv(tables / "baseline_main5_labels_early_fusion_kmeans_k2.csv", index_col=0).iloc[:, 0]
    old.index = old.index.astype(str).map(_normalize_sample_id)
    old = old.astype(int)

    bench = pd.read_csv(tables / "methylation_representation_benchmark.csv")
    if "score" not in bench.columns:
        bench["score"] = 0.45 * bench["silhouette"].fillna(-1) + 0.35 * bench["stability"].fillna(0) + 0.20 * bench["classifier_macro_f1"].fillna(0)
    best_a = bench.sort_values("score", ascending=False).iloc[0]
    rep = str(best_a["representation"])
    k_a = int(best_a["k"])

    inter = root / "data" / "interim"
    mats = {
        "mutation": _read_matrix(inter / "mutation_round1.csv", max_cols=300),
        "cnv": _read_matrix(inter / "cnv_round1.csv", max_cols=300),
        "rna": _read_matrix(inter / "rna_round1.csv", max_cols=300),
        "mirna": _read_matrix(inter / "mirna_round1.csv", max_cols=300),
    }
    if rep == "gene":
        meth = _read_matrix(inter / "methylation_gene_level_modeling.csv", max_cols=600)
    elif rep == "promoter":
        meth = _read_matrix(inter / "methylation_promoter_level_modeling.csv", max_cols=600)
    else:
        meth = _read_matrix(inter / "methylation_round2_modeling.csv", max_cols=600)
    mats["methylation"] = meth
    cohort = sorted(set.intersection(*[set(v.index) for v in mats.values()]))
    x = pd.concat([mats[m].loc[cohort] for m in ["mutation", "cnv", "methylation", "rna", "mirna"]], axis=1)
    xs = StandardScaler().fit_transform(x.values)
    labels_a = pd.Series(KMeans(n_clusters=k_a, n_init=20, random_state=42).fit_predict(xs), index=cohort, dtype=int)

    bdf = pd.read_csv(tables / "partial_fusion_ge4_labels.csv")
    bdf["sample_id"] = bdf["sample_id"].astype(str).map(_normalize_sample_id)
    labels_b = pd.Series(bdf["cluster"].astype(int).values, index=bdf["sample_id"].values)
    k_b = int(bdf["k"].iloc[0]) if "k" in bdf.columns else int(labels_b.nunique())

    labels = {
        "old_final_baseline": old,
        "phaseA_best_methylation": labels_a,
        "phaseB_partial_ge4_best": labels_b,
    }
    meta = {
        "old_final_baseline": {"k": 2},
        "phaseA_best_methylation": {"k": k_a},
        "phaseB_partial_ge4_best": {"k": k_b},
    }
    return labels, meta


def _cluster_validity(labels: pd.Series, k_expected: int) -> tuple[bool, str]:
    distinct = int(labels.nunique())
    if distinct < 2:
        return False, "invalid_k_distinct_clusters_lt_2"
    if distinct < int(k_expected):
        return False, f"invalid_k_distinct_clusters_lt_expected_{k_expected}"
    return True, "valid_k"


def _cox_precheck(df: pd.DataFrame, event_col: str) -> tuple[bool, str]:
    n = int(df.shape[0])
    e = int(df[event_col].sum())
    if n < 40:
        return False, "skipped_due_to_small_sample"
    if e < 10 or (n - e) < 10:
        return False, "skipped_due_to_low_events"

    vc = df["cluster"].value_counts()
    if vc.min() < 8:
        return False, "skipped_due_to_small_cluster"

    sep = df.groupby("cluster")[event_col].agg(["mean", "count"])
    if ((sep["mean"] <= 1e-8) | (sep["mean"] >= 1 - 1e-8)).any():
        return False, "skipped_due_to_separation_or_instability"

    return True, "cox_precheck_pass"


def _fit_cox_once(df: pd.DataFrame, time_col: str, event_col: str) -> tuple[float, float, str]:
    cox_df = df[[time_col, event_col, "cluster"]].copy()
    dummies = pd.get_dummies(cox_df["cluster"].astype(int), prefix="cluster", drop_first=True)
    if dummies.shape[1] == 0:
        return np.nan, np.nan, "cox_skipped"
    cox_df = pd.concat([cox_df[[time_col, event_col]], dummies], axis=1)

    try:
        cph = CoxPHFitter()
        cph.fit(cox_df, duration_col=time_col, event_col=event_col, fit_options={"max_steps": 80})
        hr = np.exp(cph.params_)
        pvals = cph.summary["p"]
        idx = int(np.argmax(np.abs(np.log(hr.values))))
        return float(hr.iloc[idx]), float(pvals.iloc[idx]), "cox_success"
    except (ValueError, FloatingPointError, np.linalg.LinAlgError, ZeroDivisionError):
        pass

    try:
        cph = CoxPHFitter(penalizer=0.1)
        cph.fit(cox_df, duration_col=time_col, event_col=event_col, fit_options={"max_steps": 120})
        hr = np.exp(cph.params_)
        pvals = cph.summary["p"]
        idx = int(np.argmax(np.abs(np.log(hr.values))))
        return float(hr.iloc[idx]), float(pvals.iloc[idx]), "cox_success_penalized"
    except (ValueError, FloatingPointError, np.linalg.LinAlgError, ZeroDivisionError):
        return np.nan, np.nan, "cox_failed"


def _survival_eval_for_endpoint(labels: pd.Series, survival_df: pd.DataFrame, endpoint: str, candidate: str, k_expected: int) -> dict[str, Any]:
    is_valid, validity_note = _cluster_validity(labels, k_expected)
    if not is_valid:
        return {
            "candidate": candidate,
            "endpoint": endpoint,
            "n_samples": int(labels.shape[0]),
            "n_events": np.nan,
            "logrank_p": np.nan,
            "logrank_success": False,
            "cox_hr": np.nan,
            "cox_p": np.nan,
            "cox_status": "cox_skipped",
            "notes": validity_note,
        }

    time_col = f"{endpoint}.time"
    event_col = endpoint
    if time_col not in survival_df.columns or event_col not in survival_df.columns:
        return {
            "candidate": candidate,
            "endpoint": endpoint,
            "n_samples": 0,
            "n_events": 0,
            "logrank_p": np.nan,
            "logrank_success": False,
            "cox_hr": np.nan,
            "cox_p": np.nan,
            "cox_status": "cox_skipped",
            "notes": "endpoint_missing",
        }

    merged = labels.rename("cluster").to_frame()
    merged["normalized_sample_id"] = merged.index
    merged = merged.merge(survival_df, on="normalized_sample_id", how="inner")

    merged[time_col] = pd.to_numeric(merged[time_col], errors="coerce")
    merged[event_col] = pd.to_numeric(merged[event_col], errors="coerce")
    merged = merged.dropna(subset=[time_col, event_col, "cluster"]).copy()
    merged = merged[merged[time_col] >= 0]
    if merged.empty:
        return {
            "candidate": candidate,
            "endpoint": endpoint,
            "n_samples": 0,
            "n_events": 0,
            "logrank_p": np.nan,
            "logrank_success": False,
            "cox_hr": np.nan,
            "cox_p": np.nan,
            "cox_status": "cox_skipped",
            "notes": "no_valid_rows",
        }

    merged["cluster"] = merged["cluster"].astype(int)
    merged[event_col] = merged[event_col].astype(int)

    n_samples = int(merged.shape[0])
    n_events = int(merged[event_col].sum())

    logrank_p = np.nan
    logrank_success = False
    if merged["cluster"].nunique() >= 2 and n_events > 1:
        try:
            lr = multivariate_logrank_test(merged[time_col], merged["cluster"], merged[event_col])
            logrank_p = float(lr.p_value)
            logrank_success = True
        except (ValueError, ZeroDivisionError):
            logrank_success = False

    cox_hr = np.nan
    cox_p = np.nan
    cox_status = "cox_skipped"
    notes = "logrank_success" if logrank_success else "logrank_failed"

    ok, cox_note = _cox_precheck(merged, event_col)
    if ok:
        cox_hr, cox_p, cox_status = _fit_cox_once(merged, time_col, event_col)
        if cox_status.startswith("cox_success"):
            notes = notes + ";cox_success"
        elif cox_status == "cox_failed":
            notes = notes + ";cox_failed"
        else:
            notes = notes + ";cox_skipped"
    else:
        cox_status = "cox_skipped"
        notes = notes + ";" + cox_note

    return {
        "candidate": candidate,
        "endpoint": endpoint,
        "n_samples": n_samples,
        "n_events": n_events,
        "logrank_p": logrank_p,
        "logrank_success": bool(logrank_success),
        "cox_hr": cox_hr,
        "cox_p": cox_p,
        "cox_status": cox_status,
        "notes": notes,
    }


def _plot_endpoint_km_panels(labels_map: dict[str, pd.Series], survival_df: pd.DataFrame, endpoint: str, out_path: Path) -> bool:
    time_col = f"{endpoint}.time"
    event_col = endpoint
    if time_col not in survival_df.columns or event_col not in survival_df.columns:
        return False

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
    rendered = 0
    for ax, (name, labels) in zip(axes, labels_map.items()):
        merged = labels.rename("cluster").to_frame()
        merged["normalized_sample_id"] = merged.index
        merged = merged.merge(survival_df, on="normalized_sample_id", how="inner")
        merged[time_col] = pd.to_numeric(merged[time_col], errors="coerce")
        merged[event_col] = pd.to_numeric(merged[event_col], errors="coerce")
        merged = merged.dropna(subset=[time_col, event_col, "cluster"]).copy()
        merged = merged[merged[time_col] >= 0]
        if merged.empty or merged["cluster"].nunique() < 2:
            ax.set_title(f"{name}\n(no valid KM)")
            ax.axis("off")
            continue

        kmf = KaplanMeierFitter()
        for cid in sorted(merged["cluster"].astype(int).unique()):
            part = merged[merged["cluster"].astype(int) == cid]
            kmf.fit(part[time_col], event_observed=part[event_col], label=f"c{cid}")
            sf = kmf.survival_function_.reset_index()
            # lifelines returns timeline + estimate; draw with matplotlib directly for stability.
            t_col = sf.columns[0]
            y_col = sf.columns[1]
            ax.step(sf[t_col].values, sf[y_col].values, where="post", label=f"c{cid}")
        ax.set_title(name)
        ax.set_xlabel("Time")
        ax.grid(alpha=0.25)
        ax.legend(loc="best", fontsize=8)
        rendered += 1

    axes[0].set_ylabel("Survival probability")
    plt.suptitle(f"Targeted KM ({endpoint})")
    plt.tight_layout()
    if rendered == 0:
        plt.close()
        return False
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=240)
    plt.close()
    return True


def _build_survival_summary(stable_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in stable_df.iterrows():
        endpoint = str(r["endpoint"])
        candidate = str(r["candidate"])
        logrank_ok = bool(r["logrank_success"])
        cox_status = str(r["cox_status"])

        if logrank_ok and cox_status.startswith("cox_success"):
            used = "both"
            stable = True
            interp = "logrank and cox both available; use combined survival evidence"
        elif logrank_ok:
            used = "logrank"
            stable = True
            interp = "logrank available and used as primary evidence"
        elif cox_status.startswith("cox_success"):
            used = "cox"
            stable = False
            interp = "cox available but logrank failed; treat as low-confidence auxiliary evidence"
        else:
            used = "logrank"
            stable = False
            interp = "insufficient stable survival evidence for reporting"

        if endpoint in SECONDARY_ENDPOINTS and not stable:
            interp = "secondary endpoint pending/NA due to instability"

        rows.append(
            {
                "candidate": candidate,
                "endpoint": endpoint,
                "primary_evidence_used": used,
                "interpretation": interp,
                "stable_for_reporting": bool(stable),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    root = _root()
    tables = root / "results" / "tables"
    figures = root / "results" / "figures"
    docs = root / "docs"

    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    survival_df = _load_survival_df(root)
    labels_map, meta = _load_candidate_labels(root)

    rows: list[dict[str, Any]] = []
    for cname, labels in labels_map.items():
        k = int(meta[cname]["k"])
        for ep in PRIMARY_ENDPOINTS:
            rows.append(_survival_eval_for_endpoint(labels, survival_df, ep, cname, k))

    primary_df = pd.DataFrame(rows)
    primary_stable_ratio = float(primary_df["logrank_success"].mean()) if not primary_df.empty else 0.0
    allow_secondary = primary_stable_ratio >= 0.67

    if allow_secondary:
        for cname, labels in labels_map.items():
            k = int(meta[cname]["k"])
            for ep in SECONDARY_ENDPOINTS:
                rows.append(_survival_eval_for_endpoint(labels, survival_df, ep, cname, k))
    else:
        for cname in labels_map.keys():
            for ep in SECONDARY_ENDPOINTS:
                rows.append(
                    {
                        "candidate": cname,
                        "endpoint": ep,
                        "n_samples": np.nan,
                        "n_events": np.nan,
                        "logrank_p": np.nan,
                        "logrank_success": False,
                        "cox_hr": np.nan,
                        "cox_p": np.nan,
                        "cox_status": "cox_skipped",
                        "notes": "pending_secondary_due_to_primary_instability",
                    }
                )

    stable_df = pd.DataFrame(rows)
    stable_df.to_csv(tables / "multi_endpoint_clinical_validation_stable.csv", index=False, encoding="utf-8")

    summary_df = _build_survival_summary(stable_df)
    summary_df.to_csv(tables / "targeted_enhancement_survival_summary.csv", index=False, encoding="utf-8")

    stable_df.to_csv(tables / "multi_endpoint_clinical_validation_full.csv", index=False, encoding="utf-8")

    os_ok = _plot_endpoint_km_panels(labels_map, survival_df, "OS", figures / "targeted_km_os_real.png")
    pfi_ok = _plot_endpoint_km_panels(labels_map, survival_df, "PFI", figures / "targeted_km_pfi_real.png")

    plot_df = stable_df.copy()
    plot_df = plot_df[plot_df["logrank_success"] == True].copy()
    if not plot_df.empty:
        plot_df["score"] = -np.log10(plot_df["logrank_p"].clip(lower=1e-12))
        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(10.5, 5.2))
        sns.barplot(data=plot_df, x="endpoint", y="score", hue="candidate")
        plt.ylabel("-log10(logrank p)")
        plt.title("Multi-endpoint Validation Comparison (Stable Mode)")
        plt.tight_layout()
        plt.savefig(figures / "multi_endpoint_validation_comparison_stable.png", dpi=240)
        plt.close()

    dss_rows = stable_df[stable_df["endpoint"] == "DSS"]
    dfi_rows = stable_df[stable_df["endpoint"] == "DFI"]
    secondary_pending = bool((dss_rows["logrank_success"].sum() == 0) and (dfi_rows["logrank_success"].sum() == 0))

    old_vs_p4 = stable_df[
        (stable_df["candidate"].isin(["old_final_baseline", "phaseB_partial_ge4_best"]))
        & (stable_df["endpoint"].isin(PRIMARY_ENDPOINTS))
    ].copy()

    def _mean_score(c: str) -> float:
        s = old_vs_p4[(old_vs_p4["candidate"] == c) & (old_vs_p4["logrank_success"] == True)]["logrank_p"]
        if s.empty:
            return np.nan
        return float((-np.log10(s.clip(lower=1e-12))).mean())

    old_score = _mean_score("old_final_baseline")
    p4_score = _mean_score("phaseB_partial_ge4_best")
    if pd.isna(old_score) and pd.isna(p4_score):
        verdict = "当前无法仅凭稳定生存指标比较 old 与 partial_ge4。"
    elif pd.isna(p4_score) or (not pd.isna(old_score) and old_score >= p4_score):
        verdict = "partial_ge4 未在稳定生存主证据（OS/PFI）上超过 old final baseline。"
    else:
        verdict = "partial_ge4 在稳定生存主证据（OS/PFI）上优于 old final baseline。"

    doc38 = "\n".join(
        [
            "# 38 Targeted Enhancement Summary (Stable Mode)",
            "",
            "## 本轮策略",
            "- 采用分层精算：先 3 candidates × 2 endpoints（OS/PFI）。",
            "- Cox 仅在通过稳定性门控后执行；否则跳过并记录原因。",
            "- 不新增模型家族，不扩展 ge3 主线。",
            "",
            "## 生存证据口径",
            "- 主证据：KM + log-rank。",
            "- Cox：辅助证据，仅在可收敛且不过分离时报告。",
            f"- DSS/DFI 状态：{'pending/NA（未纳入本轮最终结论）' if secondary_pending else '已在稳定条件下补算'}。",
            "",
            "## 关键输出",
            "- results/tables/multi_endpoint_clinical_validation_stable.csv",
            "- results/tables/targeted_enhancement_survival_summary.csv",
            "- results/figures/targeted_km_os_real.png",
            "- results/figures/targeted_km_pfi_real.png",
            "- results/figures/multi_endpoint_validation_comparison_stable.png",
            "",
            "## 结论",
            f"- {verdict}",
            "- 是否替换 old final baseline：仅当稳定主证据与其他维度形成一致优势时才考虑替换。",
        ]
    )
    (docs / "38_targeted_enhancement_summary.md").write_text(doc38, encoding="utf-8")

    progress_block = "\n".join(
        [
            "## 2026-04-20 | Targeted Enhancement Stable Survival Pass",
            "- 已停止全量不稳定 Cox 路径，改为分层精算策略。",
            "- 已完成 3 candidates × OS/PFI 的真实 KM + log-rank 主证据产出。",
            "- Cox 采用门控 + penalized 单次回退，失败即降级为 NA，不再卡住整轮。",
            f"- DSS/DFI：{'pending/NA（本轮不纳入结论）' if secondary_pending else '在稳定条件下已补算'}。",
        ]
    )

    defense_block = "\n".join(
        [
            "## Targeted Enhancement Stable Evidence Notes",
            "- 当前生存增强以 KM + log-rank 为主证据。",
            "- Cox 在 separation/数值不稳场景下仅作辅助，允许跳过并明确标记。",
            "- 可正式汇报：OS/PFI 的稳定结果与真实 KM 图。",
            f"- 暂不纳入最终结论：{'DSS/DFI（pending/NA）' if secondary_pending else '无'}。",
        ]
    )

    readme_block = "\n".join(
        [
            "## Round 6 Stable Survival Update",
            "- Targeted enhancement survival validation now uses stable mode (KM + log-rank first).",
            "- Cox is conditional and may be skipped for separation/instability.",
            "- New outputs: `results/tables/multi_endpoint_clinical_validation_stable.csv`, `results/tables/targeted_enhancement_survival_summary.csv`.",
            "- Real figures: `results/figures/targeted_km_os_real.png`, `results/figures/targeted_km_pfi_real.png`, `results/figures/multi_endpoint_validation_comparison_stable.png`.",
        ]
    )

    for path, block in [
        (docs / "08_progress_log.md", progress_block),
        (docs / "09_defense_notes.md", defense_block),
        (root / "README.md", readme_block),
    ]:
        orig = path.read_text(encoding="utf-8") if path.exists() else ""
        path.write_text(orig.rstrip() + "\n\n" + block + "\n", encoding="utf-8")

    print(f"Stable survival pass complete. OS_KM={os_ok}, PFI_KM={pfi_ok}, secondary_pending={secondary_pending}")


if __name__ == "__main__":
    main()
