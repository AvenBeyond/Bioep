"""Final system entrypoint for STAD multi-omics subtype discovery and assignment.

Modes:
- --summary: print final system summary
- --report: generate dashboard table/figures and report-ready checks
- --predict-subtype: prototype subtype assignment for new samples
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Support both `python -m src.pipelines.run_final_system` and
# `python src/pipelines/run_final_system.py` invocation styles.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.evaluation.clinical_association import plot_km_curve


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _normalize_sample_id(sample_id: str) -> str:
    sid = str(sample_id).strip().upper()
    return sid[:16] if len(sid) >= 16 else sid


def _read_interim(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).map(_normalize_sample_id)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def _load_final_tables(root: Path) -> dict[str, pd.DataFrame]:
    tables = root / "results" / "tables"
    return {
        "final_cmp": pd.read_csv(tables / "final_model_comparison_summary.csv"),
        "key": pd.read_csv(tables / "final_key_numbers_for_ppt.csv"),
        "baseline": pd.read_csv(tables / "baseline_main5_metrics.csv"),
        "proposed": pd.read_csv(tables / "proposed_main5_metrics.csv"),
        "improve": pd.read_csv(tables / "improvement_model_comparison.csv"),
        "clf": pd.read_csv(tables / "subtype_classifier_cv_results.csv"),
        "ablation_rank": pd.read_csv(tables / "ablation_rank_summary.csv"),
    }


def _best_row(df: pd.DataFrame, metric: str = "silhouette") -> pd.Series:
    tmp = df.copy()
    tmp[metric] = pd.to_numeric(tmp[metric], errors="coerce")
    return tmp.sort_values(metric, ascending=False).iloc[0]


def _get_final_main_result(final_cmp: pd.DataFrame) -> pd.Series:
    s = final_cmp[final_cmp["selected_as_final"] == True]
    if s.empty:
        return _best_row(final_cmp, metric="silhouette")
    return s.iloc[0]


def _build_master_dashboard(root: Path, tables_in: dict[str, pd.DataFrame]) -> pd.DataFrame:
    final_cmp = tables_in["final_cmp"]
    proposed = tables_in["proposed"]
    improve = tables_in["improve"]
    clf = tables_in["clf"]

    main_row = _get_final_main_result(final_cmp)
    proposed_best = _best_row(proposed, metric="silhouette")
    improve_best = improve[improve["selected_as_best_improvement"] == True]
    if improve_best.empty:
        improve_best = pd.DataFrame([_best_row(improve, metric="silhouette")])
    improve_best = improve_best.iloc[0]

    rows: list[dict[str, Any]] = []

    def add(category: str, method: str, modality_set: str, k: Any, metric_name: str, metric_value: Any, selected: bool, interpretation: str, notes: str) -> None:
        rows.append(
            {
                "category": category,
                "method": method,
                "modality_set": modality_set,
                "k": k,
                "metric_name": metric_name,
                "metric_value": metric_value,
                "selected_as_main_result": selected,
                "interpretation": interpretation,
                "notes": notes,
            }
        )

    # Clustering quality
    add("clustering quality", str(main_row["method"]), str(main_row["modality_set"]), int(main_row["k"]), "silhouette", float(main_row["silhouette"]), True, "higher is better; final main result", "final_main")
    add("clustering quality", "proposed_weighted_fusion", "main5", int(proposed_best["k"]), "silhouette", float(proposed_best["silhouette"]), False, "comparison to final main", "proposed_baseline_compare")
    add("clustering quality", str(improve_best["method"]), str(improve_best["cohort"]) if "cohort" in improve_best else "main5", int(improve_best["k"]), "silhouette", float(improve_best["silhouette"]), False, "best improvement branch candidate", "improvement_compare")

    # Robustness
    add("robustness", str(main_row["method"]), str(main_row["modality_set"]), int(main_row["k"]), "stability", float(main_row["stability"]), True, "higher is better", "from final comparison")
    add("robustness", str(improve_best["method"]), str(improve_best.get("cohort", "main5")), int(improve_best["k"]), "consensus_stability", float(improve_best.get("consensus_stability", np.nan)), False, "higher is better", "from improvement comparison")

    # Clinical relevance
    add("clinical relevance", str(main_row["method"]), str(main_row["modality_set"]), int(main_row["k"]), "logrank_os_p", float(main_row["clinical_separation_os"]), True, "lower is better; significant when <0.05", "from final comparison")
    add("clinical relevance", str(main_row["method"]), str(main_row["modality_set"]), int(main_row["k"]), "logrank_pfi_p", float(main_row["clinical_separation_pfi"]), True, "lower is better; significant when <0.05", "from final comparison")
    add("clinical relevance", str(improve_best["method"]), str(improve_best.get("cohort", "main5")), int(improve_best["k"]), "logrank_os_p", float(improve_best.get("logrank_os_p", np.nan)), False, "lower is better", "from improvement comparison")

    # Classification performance
    clf_valid = clf[pd.to_numeric(clf["macro_f1"], errors="coerce").notna()].copy()
    clf_valid["macro_f1"] = pd.to_numeric(clf_valid["macro_f1"], errors="coerce")
    clf_valid["accuracy"] = pd.to_numeric(clf_valid["accuracy"], errors="coerce")
    clf_best = clf_valid.sort_values("macro_f1", ascending=False).iloc[0]

    add("classification performance", str(clf_best["model"]), str(clf_best["feature_set"]), "NA", "macro_f1", float(clf_best["macro_f1"]), True, "higher is better", "internal validation only")
    add("classification performance", str(clf_best["model"]), str(clf_best["feature_set"]), "NA", "accuracy", float(clf_best["accuracy"]), True, "higher is better", "internal validation only")

    out = pd.DataFrame(rows)
    out_path = root / "results" / "tables" / "final_dashboard_master.csv"
    out.to_csv(out_path, index=False, encoding="utf-8")
    return out


def _make_dashboard_overview(root: Path, dashboard: pd.DataFrame) -> None:
    fig_path = root / "results" / "figures" / "final_dashboard_overview.png"
    category_metric = {
        "clustering quality": "silhouette",
        "robustness": "stability",
        "clinical relevance": "logrank_os_p",
        "classification performance": "macro_f1",
    }

    vals = []
    for cat, m in category_metric.items():
        x = dashboard[(dashboard["category"] == cat) & (dashboard["metric_name"] == m) & (dashboard["selected_as_main_result"] == True)]
        if x.empty:
            vals.append(np.nan)
        else:
            vals.append(float(x.iloc[0]["metric_value"]))

    # Normalize to [0,1] and invert p-value axis.
    arr = np.array(vals, dtype=float)
    arr_plot = arr.copy()
    if np.isfinite(arr_plot[2]):
        arr_plot[2] = -np.log10(max(arr_plot[2], 1e-12))
    for i in range(arr_plot.size):
        if not np.isfinite(arr_plot[i]):
            arr_plot[i] = 0.0
    mn, mx = arr_plot.min(), arr_plot.max()
    if mx - mn > 1e-12:
        arr_norm = (arr_plot - mn) / (mx - mn)
    else:
        arr_norm = np.zeros_like(arr_plot)

    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(9, 5))
    cats = list(category_metric.keys())
    plt.bar(cats, arr_norm, color=["#3a86ff", "#2a9d8f", "#f4a261", "#e76f51"])
    plt.ylim(0, 1)
    plt.ylabel("normalized score")
    plt.title("Final Dashboard Overview: Four Evaluation Dimensions")
    plt.xticks(rotation=12, ha="right")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=240)
    plt.close()


def _make_method_comparison_compact(root: Path, tables_in: dict[str, pd.DataFrame]) -> None:
    fig_path = root / "results" / "figures" / "final_method_comparison_compact.png"
    final_cmp = tables_in["final_cmp"]
    improve = tables_in["improve"]
    clf = tables_in["clf"]

    main_row = _get_final_main_result(final_cmp)
    proposed = final_cmp[final_cmp["method"] == "proposed_weighted_fusion"]
    proposed_row = proposed.iloc[0] if not proposed.empty else None
    best_imp = improve[improve["selected_as_best_improvement"] == True]
    best_imp = best_imp.iloc[0] if not best_imp.empty else _best_row(improve)

    clf_valid = clf[pd.to_numeric(clf["macro_f1"], errors="coerce").notna()].copy()
    clf_valid["macro_f1"] = pd.to_numeric(clf_valid["macro_f1"], errors="coerce")
    clf_best = clf_valid.sort_values("macro_f1", ascending=False).iloc[0]

    rows = [
        {
            "candidate": "final baseline",
            "silhouette": float(main_row["silhouette"]),
            "stability": float(main_row["stability"]),
            "clinical": -np.log10(max(float(main_row["clinical_separation_os"]), 1e-12)),
            "classification": float(clf_best["macro_f1"]),
        },
        {
            "candidate": "proposed weighted fusion",
            "silhouette": float(proposed_row["silhouette"]) if proposed_row is not None else np.nan,
            "stability": float(proposed_row["stability"]) if proposed_row is not None and pd.notna(proposed_row["stability"]) else 0.0,
            "clinical": -np.log10(max(float(proposed_row["clinical_separation_os"]), 1e-12)) if proposed_row is not None else np.nan,
            "classification": np.nan,
        },
        {
            "candidate": "improvement best",
            "silhouette": float(best_imp["silhouette"]),
            "stability": float(best_imp.get("consensus_stability", np.nan)),
            "clinical": -np.log10(max(float(best_imp.get("logrank_os_p", np.nan)), 1e-12)) if pd.notna(best_imp.get("logrank_os_p", np.nan)) else np.nan,
            "classification": np.nan,
        },
    ]
    cmp = pd.DataFrame(rows).set_index("candidate")
    cmp = cmp.fillna(0.0)

    # Min-max by column for compact visual comparison.
    norm = cmp.copy()
    for c in norm.columns:
        v = norm[c].values.astype(float)
        lo, hi = float(v.min()), float(v.max())
        norm[c] = (v - lo) / (hi - lo) if hi - lo > 1e-12 else 0.0

    plt.figure(figsize=(10, 5))
    norm.plot(kind="bar", ax=plt.gca(), width=0.75)
    plt.title("Compact Method Comparison (normalized across key dimensions)")
    plt.ylabel("normalized score")
    plt.xticks(rotation=12, ha="right")
    plt.legend(loc="upper right", ncol=2)
    plt.tight_layout()
    plt.savefig(fig_path, dpi=240)
    plt.close()


def _make_pipeline_flowchart(root: Path) -> None:
    fig_path = root / "results" / "figures" / "final_pipeline_flowchart.png"

    plt.figure(figsize=(14, 4))
    ax = plt.gca()
    ax.axis("off")

    boxes = [
        "raw data",
        "preprocess",
        "subtype discovery",
        "model comparison",
        "final subtype labels",
        "subtype classifier",
        "new sample assignment",
    ]
    x_positions = np.linspace(0.05, 0.95, len(boxes))
    y = 0.5

    for i, (x, txt) in enumerate(zip(x_positions, boxes)):
        ax.text(
            x,
            y,
            txt,
            ha="center",
            va="center",
            fontsize=11,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "#edf2f4", "edgecolor": "#2b2d42"},
            transform=ax.transAxes,
        )
        if i < len(boxes) - 1:
            ax.annotate(
                "",
                xy=(x_positions[i + 1] - 0.05, y),
                xytext=(x + 0.05, y),
                arrowprops={"arrowstyle": "->", "lw": 1.8, "color": "#2b2d42"},
                xycoords=ax.transAxes,
                textcoords=ax.transAxes,
            )

    plt.title("Final System Pipeline: Discovery + Assignment", fontsize=13)
    plt.tight_layout()
    plt.savefig(fig_path, dpi=240)
    plt.close()


def _ensure_final_km_curve(root: Path) -> Path | None:
    labels_path = root / "results" / "tables" / "baseline_main5_labels_early_fusion_kmeans_k2.csv"
    survival_path = root / "data" / "raw" / "clinical" / "survival_STAD_survival.txt"
    out = root / "results" / "figures" / "final_km_curve_os.png"

    if not labels_path.exists() or not survival_path.exists():
        return None

    labels = pd.read_csv(labels_path, index_col=0).iloc[:, 0]
    labels.index = labels.index.astype(str).map(_normalize_sample_id)
    labels.name = "cluster"
    survival_df = pd.read_csv(survival_path, sep="\t", dtype=str)

    plot_km_curve(labels=labels, survival_df=survival_df, endpoint="OS", out_path=out)
    return out


def _paste_image(ax, path: Path, title: str) -> None:
    ax.axis("off")
    if path.exists():
        img = plt.imread(path)
        ax.imshow(img)
        ax.set_title(title, fontsize=10)
    else:
        ax.text(0.5, 0.5, f"missing\n{path.name}", ha="center", va="center", fontsize=10)
        ax.set_title(title, fontsize=10)


def _make_final_storyboard(root: Path) -> None:
    fig_path = root / "results" / "figures" / "final_result_storyboard.png"
    figs = root / "results" / "figures"
    km_path = _ensure_final_km_curve(root)

    fig = plt.figure(figsize=(14, 8), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, wspace=0.05, hspace=0.20)

    ax1 = fig.add_subplot(gs[0, 0])
    _paste_image(ax1, figs / "main_vs_with_rppa_intersection.png", "Data Scale")

    ax2 = fig.add_subplot(gs[0, 1])
    _paste_image(ax2, figs / "baseline_main5_metric_comparison.png", "Main Model Context")

    ax3 = fig.add_subplot(gs[0, 2])
    _paste_image(ax3, km_path if km_path is not None else figs / "km_early_fusion_kmeans_main4_bestk_os.png", "KM Curve")

    ax4 = fig.add_subplot(gs[1, 0])
    _paste_image(ax4, figs / "ablation_comparison.png", "Ablation")

    ax5 = fig.add_subplot(gs[1, 1])
    _paste_image(ax5, figs / "subtype_classifier_performance.png", "Classifier")

    ax6 = fig.add_subplot(gs[1, 2])
    _paste_image(ax6, figs / "final_method_comparison_compact.png", "System Summary")

    fig.suptitle("Final Result Storyboard: Discovery and Assignment System", fontsize=14)
    fig.savefig(fig_path, dpi=240)
    plt.close(fig)


def generate_report_assets(root: Path) -> dict[str, Any]:
    tables = _load_final_tables(root)
    dashboard = _build_master_dashboard(root, tables)
    _make_dashboard_overview(root, dashboard)
    _make_method_comparison_compact(root, tables)
    _make_pipeline_flowchart(root)
    _make_final_storyboard(root)

    return {
        "dashboard_rows": int(dashboard.shape[0]),
        "final_main_method": str(_get_final_main_result(tables["final_cmp"])["method"]),
    }


def _load_main5_concat(root: Path) -> pd.DataFrame:
    interim = root / "data" / "interim"
    mats = {
        "mutation": _read_interim(interim / "mutation_round1.csv"),
        "cnv": _read_interim(interim / "cnv_round1.csv"),
        "methylation": _read_interim(interim / "methylation_round2_modeling.csv"),
        "rna": _read_interim(interim / "rna_round1.csv"),
        "mirna": _read_interim(interim / "mirna_round1.csv"),
    }
    cohort = sorted(set.intersection(*[set(m.index) for m in mats.values()]))

    # Prefix columns to avoid collisions.
    aligned = []
    for mod, df in mats.items():
        x = df.loc[cohort].copy()
        x.columns = [f"{mod}::{c}" for c in x.columns]
        aligned.append(x)
    return pd.concat(aligned, axis=1)


def predict_subtype_prototype(root: Path, new_input_csv: Path, out_csv: Path | None) -> pd.DataFrame:
    x_train = _load_main5_concat(root)
    labels_path = root / "results" / "tables" / "baseline_main5_labels_early_fusion_kmeans_k2.csv"
    y = pd.read_csv(labels_path, index_col=0).iloc[:, 0]
    y.index = y.index.astype(str).map(_normalize_sample_id)

    # Align train/labels index.
    common = sorted(set(x_train.index) & set(y.index))
    x_train = x_train.loc[common]
    y = y.loc[common]

    x_new = pd.read_csv(new_input_csv, index_col=0)
    x_new.index = x_new.index.astype(str).map(_normalize_sample_id)

    common_cols = [c for c in x_train.columns if c in x_new.columns]
    if len(common_cols) < 50:
        raise ValueError("predict-subtype prototype requires at least 50 overlapping prefixed features.")

    x_tr = x_train[common_cols].values
    x_ne = x_new[common_cols].fillna(0.0).values

    centroids = {}
    for cls in sorted(y.unique()):
        centroids[int(cls)] = x_tr[y.values == cls].mean(axis=0)

    pred = []
    for i in range(x_ne.shape[0]):
        d = {cls: float(np.linalg.norm(x_ne[i] - c)) for cls, c in centroids.items()}
        best = sorted(d.items(), key=lambda kv: kv[1])[0][0]
        pred.append(int(best))

    out = pd.DataFrame(
        {
            "sample_id": x_new.index,
            "predicted_subtype": pred,
            "model_used": "prototype_nearest_centroid_from_final_discovery_labels",
            "disclaimer": "research_prototype_only_not_for_clinical_diagnosis",
        }
    )

    if out_csv is None:
        out_csv = root / "results" / "tables" / "prototype_predicted_subtypes.csv"
    out.to_csv(out_csv, index=False, encoding="utf-8")
    return out


def print_summary(root: Path) -> None:
    tables = _load_final_tables(root)
    final_cmp = tables["final_cmp"]
    key = tables["key"]
    clf = tables["clf"]

    main_row = _get_final_main_result(final_cmp)
    clf_valid = clf[pd.to_numeric(clf["macro_f1"], errors="coerce").notna()].copy()
    clf_valid["macro_f1"] = pd.to_numeric(clf_valid["macro_f1"], errors="coerce")
    clf_best = clf_valid.sort_values("macro_f1", ascending=False).iloc[0]

    key_map = {str(r["key"]): r["value"] for _, r in key.iterrows()}

    print("=== Final System Summary ===")
    print("System name: 基于多组学数据的胃腺癌分型与分型预测系统")
    print(f"Final main model: {main_row['method']} (modality_set={main_row['modality_set']}, k={int(main_row['k'])})")
    print(f"Main metrics: silhouette={float(main_row['silhouette']):.4f}, stability={float(main_row['stability']):.4f}")
    print(f"Clinical separation: OS p={float(main_row['clinical_separation_os']):.6g}, PFI p={float(main_row['clinical_separation_pfi']):.6g}")
    print(f"Classifier best (internal): {clf_best['model']} macro_f1={float(clf_best['macro_f1']):.4f}, accuracy={float(clf_best['accuracy']):.4f}")
    print(f"Main cohort size: {key_map.get('main5_intersection_samples', 'NA')}")
    print("Known limitations: no external validation; warning-prone spectral connectivity in some branches; research prototype assignment.")


def run_report(root: Path) -> None:
    info = generate_report_assets(root)

    required_tables = [
        root / "results" / "tables" / "final_model_comparison_summary.csv",
        root / "results" / "tables" / "final_key_numbers_for_ppt.csv",
        root / "results" / "tables" / "final_dashboard_master.csv",
    ]
    required_figs = [
        root / "results" / "figures" / "final_dashboard_overview.png",
        root / "results" / "figures" / "final_method_comparison_compact.png",
        root / "results" / "figures" / "final_pipeline_flowchart.png",
        root / "results" / "figures" / "final_result_storyboard.png",
    ]
    required_labels = [
        root / "results" / "tables" / "baseline_main5_labels_early_fusion_kmeans_k2.csv",
    ]

    def _check(paths: list[Path]) -> tuple[int, list[str]]:
        missing = [str(p.relative_to(root)) for p in paths if not p.exists()]
        return len(missing), missing

    mt, missing_tables = _check(required_tables)
    mf, missing_figs = _check(required_figs)
    ml, missing_labels = _check(required_labels)

    print("=== Report-Ready Asset Check ===")
    print(f"dashboard_rows={info['dashboard_rows']}")
    print(f"final_main_method={info['final_main_method']}")
    print(f"missing_tables={mt}")
    print(f"missing_figures={mf}")
    print(f"missing_label_files={ml}")
    if missing_tables:
        print("tables_missing_list=")
        for x in missing_tables:
            print("-", x)
    if missing_figs:
        print("figures_missing_list=")
        for x in missing_figs:
            print("-", x)
    if missing_labels:
        print("labels_missing_list=")
        for x in missing_labels:
            print("-", x)


def _build_innovation_report_tables(root: Path, tables_in: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, pd.DataFrame]:
    final_cmp = tables_in["final_cmp"].copy()
    proposed = tables_in["proposed"].copy()
    improve = tables_in["improve"].copy()
    clf = tables_in["clf"].copy()

    final_main = _get_final_main_result(final_cmp)

    proposed["silhouette"] = pd.to_numeric(proposed["silhouette"], errors="coerce")
    proposed_best = proposed.sort_values("silhouette", ascending=False).iloc[0]

    improve_best = improve[improve["selected_as_best_improvement"] == True]
    if improve_best.empty:
        improve_best = pd.DataFrame([_best_row(improve, metric="silhouette")])
    improve_best = improve_best.iloc[0]

    clf_valid = clf[pd.to_numeric(clf["macro_f1"], errors="coerce").notna()].copy()
    clf_valid["macro_f1"] = pd.to_numeric(clf_valid["macro_f1"], errors="coerce")
    clf_valid["accuracy"] = pd.to_numeric(clf_valid["accuracy"], errors="coerce")
    clf_best = clf_valid.sort_values("macro_f1", ascending=False).iloc[0]

    rows = [
        {
            "model_track": "final_main_baseline",
            "method": str(final_main["method"]),
            "modality_set": str(final_main["modality_set"]),
            "k": int(final_main["k"]),
            "silhouette": float(final_main["silhouette"]),
            "stability": float(final_main["stability"]),
            "os_p": float(final_main["clinical_separation_os"]),
            "pfi_p": float(final_main["clinical_separation_pfi"]),
            "cluster_balance": float(final_main["cluster_balance"]),
            "classifier_macro_f1": float(clf_best["macro_f1"]),
            "selected_in_innovation_report": True,
            "interpretation": "reference main result for comparison",
        },
        {
            "model_track": "innovation_proposed_weighted_fusion",
            "method": "proposed_weighted_fusion",
            "modality_set": str(proposed_best.get("modality_set", "main5")),
            "k": int(proposed_best["k"]),
            "silhouette": float(proposed_best["silhouette"]),
            "stability": np.nan,
            "os_p": float(final_cmp.loc[final_cmp["method"] == "proposed_weighted_fusion", "clinical_separation_os"].iloc[0]) if (final_cmp["method"] == "proposed_weighted_fusion").any() else np.nan,
            "pfi_p": float(final_cmp.loc[final_cmp["method"] == "proposed_weighted_fusion", "clinical_separation_pfi"].iloc[0]) if (final_cmp["method"] == "proposed_weighted_fusion").any() else np.nan,
            "cluster_balance": float(proposed_best.get("cluster_balance_score", np.nan)),
            "classifier_macro_f1": np.nan,
            "selected_in_innovation_report": True,
            "interpretation": "innovation mechanism validated, but not replacing final main",
        },
        {
            "model_track": "improvement_best",
            "method": str(improve_best["method"]),
            "modality_set": str(improve_best.get("cohort", "main5")),
            "k": int(improve_best["k"]),
            "silhouette": float(improve_best["silhouette"]),
            "stability": float(improve_best.get("consensus_stability", np.nan)),
            "os_p": float(improve_best.get("logrank_os_p", np.nan)),
            "pfi_p": float(improve_best.get("logrank_pfi_p", np.nan)),
            "cluster_balance": float(improve_best.get("cluster_balance_score", np.nan)),
            "classifier_macro_f1": np.nan,
            "selected_in_innovation_report": True,
            "interpretation": "best balanced candidate in improvement round",
        },
    ]

    out_df = pd.DataFrame(rows)
    out_df.to_csv(root / "results" / "tables" / "final_innovation_model_report.csv", index=False, encoding="utf-8")

    takeaways = pd.DataFrame(
        [
            {"item": "innovation_best_candidate", "value": str(improve_best["method"])},
            {"item": "proposed_weighted_fusion_role", "value": "innovation evidence track; not final replacement"},
            {"item": "does_innovation_replace_main", "value": "no"},
            {"item": "final_statement", "value": "system keeps baseline main result while preserving innovation/improvement evidence"},
        ]
    )
    takeaways.to_csv(root / "results" / "tables" / "final_innovation_key_takeaways.csv", index=False, encoding="utf-8")
    return out_df, takeaways


def _make_innovation_comparison_figure(root: Path, report_df: pd.DataFrame) -> None:
    fig_path = root / "results" / "figures" / "final_innovation_comparison.png"
    plot_df = report_df[["model_track", "silhouette", "cluster_balance"]].copy()
    for c in ["silhouette", "cluster_balance"]:
        v = pd.to_numeric(plot_df[c], errors="coerce").fillna(0.0)
        lo, hi = float(v.min()), float(v.max())
        plot_df[c] = (v - lo) / (hi - lo) if hi - lo > 1e-12 else 0.0

    plt.figure(figsize=(9, 4.5))
    x = np.arange(len(plot_df))
    w = 0.35
    plt.bar(x - w / 2, plot_df["silhouette"], width=w, label="silhouette(norm)", color="#3a86ff")
    plt.bar(x + w / 2, plot_df["cluster_balance"], width=w, label="balance(norm)", color="#2a9d8f")
    plt.xticks(x, plot_df["model_track"], rotation=10, ha="right")
    plt.ylim(0, 1)
    plt.title("Innovation-Focused Comparison (normalized)")
    plt.ylabel("normalized score")
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_path, dpi=240)
    plt.close()


def _write_innovation_report_doc(root: Path, report_df: pd.DataFrame, _takeaways_df: pd.DataFrame) -> None:
    out = root / "docs" / "33_innovation_improvement_final_report.md"
    main_row = report_df.loc[report_df["model_track"] == "final_main_baseline"].iloc[0]
    imp_row = report_df.loc[report_df["model_track"] == "improvement_best"].iloc[0]
    final_statement = str(
        _takeaways_df.loc[_takeaways_df["item"] == "final_statement", "value"].iloc[0]
    ) if (_takeaways_df["item"] == "final_statement").any() else "innovation evidence retained as supporting track"
    text = f"""# 33 创新改进方法版结果与最终报告

## 报告目的
本报告给出“创新/改进方法视角”的结果整合版本，用于与主线最终结果并行展示。

## 核心结论
- 主线最终结果仍为：{main_row['method']} (modality_set={main_row['modality_set']}, k={int(main_row['k'])})。
- 创新改进分支最佳候选：{imp_row['method']} (cohort={imp_row['modality_set']}, k={int(imp_row['k'])})。
- 结论口径：{final_statement}。

## 结果文件
- results/tables/final_innovation_model_report.csv
- results/tables/final_innovation_key_takeaways.csv
- results/figures/final_innovation_comparison.png

## 使用建议（答辩）
1. 先报告主线最终系统结论（保证稳健与可复现）。
2. 再报告创新/改进分支的增益与边界（强调真实对照与不强改结论）。
3. 明确 assignment 当前为研究原型而非临床部署。
"""
    out.write_text(text, encoding="utf-8")


def print_innovation_summary(root: Path) -> None:
    tables = _load_final_tables(root)
    rep, tk = _build_innovation_report_tables(root, tables)
    _make_innovation_comparison_figure(root, rep)
    _write_innovation_report_doc(root, rep, tk)

    print("=== Innovation-Focused Summary ===")
    print("Generated: results/tables/final_innovation_model_report.csv")
    print("Generated: results/tables/final_innovation_key_takeaways.csv")
    print("Generated: results/figures/final_innovation_comparison.png")
    print("Generated: docs/33_innovation_improvement_final_report.md")


def run_report_innovation(root: Path) -> None:
    print_innovation_summary(root)
    required = [
        root / "results" / "tables" / "final_innovation_model_report.csv",
        root / "results" / "tables" / "final_innovation_key_takeaways.csv",
        root / "results" / "figures" / "final_innovation_comparison.png",
        root / "docs" / "33_innovation_improvement_final_report.md",
    ]
    missing = [str(p.relative_to(root)) for p in required if not p.exists()]
    print("=== Innovation Report Asset Check ===")
    print(f"missing_count={len(missing)}")
    if missing:
        for x in missing:
            print("-", x)


def main() -> None:
    parser = argparse.ArgumentParser(description="Final STAD system entrypoint")
    parser.add_argument("--summary", action="store_true", help="print final system summary")
    parser.add_argument("--report", action="store_true", help="generate dashboard and run report-ready checks")
    parser.add_argument("--summary-innovation", action="store_true", help="print and generate innovation-focused summary artifacts")
    parser.add_argument("--report-innovation", action="store_true", help="generate innovation-focused report files and checks")
    parser.add_argument("--predict-subtype", action="store_true", help="prototype subtype assignment")
    parser.add_argument("--input", type=str, default=None, help="input csv for --predict-subtype")
    parser.add_argument("--output", type=str, default=None, help="output csv for --predict-subtype")
    args = parser.parse_args()

    root = _project_root()

    if args.summary:
        print_summary(root)
        return

    if args.report:
        run_report(root)
        return

    if args.summary_innovation:
        print_innovation_summary(root)
        return

    if args.report_innovation:
        run_report_innovation(root)
        return

    if args.predict_subtype:
        if args.input is None:
            raise ValueError("--predict-subtype requires --input <csv>")
        out_df = predict_subtype_prototype(
            root=root,
            new_input_csv=Path(args.input),
            out_csv=Path(args.output) if args.output else None,
        )
        print("=== Predict Subtype Prototype Output ===")
        print(out_df.head().to_string(index=False))
        return

    parser.print_help()


if __name__ == "__main__":
    main()
