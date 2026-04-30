"""Phase 11-15 runner.

This script executes:
1) modeling readiness audit
2) similarity building
3) baseline round1 runs
4) stability and clinical association round1
5) proposed-method readiness preview
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

from src.clustering.baseline_early_fusion import evaluate_early_fusion_round1
from src.clustering.baseline_snf import evaluate_equal_weight_fusion_round1
from src.clustering.proposed_weighted_fusion import build_weight_input_preview, dry_run_weight_estimation
from src.evaluation.clinical_association import (
    evaluate_clinical_association,
    evaluate_clinical_variables,
    plot_clinical_association_heatmap,
    plot_km_curve,
)
from src.evaluation.cluster_metrics import estimate_stability_by_subsampling
from src.feature_engineering.build_similarity_matrices import build_similarity_round1
from src.utils.config_utils import load_project_config
from src.utils.io_utils import ensure_dir


def _read_interim(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).str.upper().str[:16]
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def build_modeling_readiness_summary(project_root: Path) -> tuple[pd.DataFrame, list[str], str]:
    paths_cfg, _ = load_project_config(project_root)
    interim_root = Path(paths_cfg["data"]["interim_root"])
    tables_dir = Path(paths_cfg["results"]["tables"])
    ensure_dir(tables_dir)

    files = {
        "mutation": interim_root / "mutation_round1.csv",
        "cnv": interim_root / "cnv_round1.csv",
        "methylation": interim_root / "methylation_round1_summary.csv",
        "rna": interim_root / "rna_round1.csv",
        "mirna": interim_root / "mirna_round1.csv",
        "rppa": interim_root / "rppa_round1.csv",
    }

    rows: list[dict[str, Any]] = []
    ready_modalities: list[str] = []
    for m, p in files.items():
        exists = p.exists()
        ready = False
        sample_count = 0
        feature_count = 0
        reason = ""
        missingness_level = "unknown"
        norm_status = "unknown"
        alignment_status = "unknown"

        if exists:
            if m == "methylation":
                # Current output is summary/preview, not full modeling matrix.
                summary = pd.read_csv(p)
                sample_count = int(summary.loc[0, "sample_count_kept"])
                feature_count = int(summary.loc[0, "top_var_preview_features"])
                missingness_level = f"{float(summary.loc[0, 'missing_rate']) * 100:.3f}%"
                ready = False
                reason = "round1 only has chunked summary + preview, no full sample-feature matrix"
                norm_status = "not_applicable"
                alignment_status = "partial_ready"
            else:
                df = _read_interim(p)
                sample_count, feature_count = df.shape
                miss = df.isna().mean().mean() * 100
                missingness_level = f"{miss:.3f}%"
                norm_status = "zscore_or_processed"
                alignment_status = "sample_id_normalized"
                ready = True
                reason = ""
                ready_modalities.append(m)

        rows.append(
            {
                "modality": m,
                "interim_file": str(p),
                "exists": bool(exists),
                "matrix_ready_for_modeling": bool(ready),
                "sample_count": sample_count,
                "feature_count": feature_count,
                "missingness_level": missingness_level,
                "normalization_status": norm_status,
                "patient_alignment_status": alignment_status,
                "can_be_used_in_round1_baseline": bool(ready and m in ["mutation", "cnv", "rna", "mirna"]),
                "reason_if_not_ready": reason,
                "notes": "phase11_audit",
            }
        )

    readiness_df = pd.DataFrame(rows)
    readiness_df.to_csv(tables_dir / "modeling_readiness_summary.csv", index=False, encoding="utf-8")

    if all(m in ready_modalities for m in ["mutation", "cnv", "rna", "mirna", "methylation"]):
        baseline_modalities = ["mutation", "cnv", "methylation", "rna", "mirna"]
        plan = "main5"
    else:
        baseline_modalities = ["mutation", "cnv", "rna", "mirna"]
        plan = "main4_without_methylation"

    return readiness_df, baseline_modalities, plan


def _plot_baseline_comparison(metrics: pd.DataFrame, out_fig: Path) -> None:
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(9, 5))
    sns.lineplot(data=metrics, x="k", y="silhouette", hue="method", marker="o")
    plt.title("Baseline Round1 Silhouette Comparison")
    plt.tight_layout()
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_fig, dpi=220)
    plt.close()


def _plot_cluster_size(cluster_sizes: pd.DataFrame, out_fig: Path) -> None:
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 5))
    sns.barplot(data=cluster_sizes, x="k", y="cluster_size", hue="method")
    plt.title("Cluster Size Comparison (Round1)")
    plt.tight_layout()
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_fig, dpi=220)
    plt.close()


def _method_recluster_fn(k: int, random_state: int = 42):
    # Round1 simplified stability: unified KMeans reclustering for speed/robustness.
    def fn(x_sub: np.ndarray) -> np.ndarray:
        return KMeans(n_clusters=k, n_init=20, random_state=random_state).fit_predict(x_sub)

    return fn


def _run_stability_round1(
    metrics_df: pd.DataFrame,
    labels_map: dict[tuple[str, int], pd.Series],
    x_scaled: np.ndarray,
    modality_set_name: str,
    out_table: Path,
) -> pd.DataFrame:
    pca_dim = min(50, x_scaled.shape[1], x_scaled.shape[0] - 1)
    x_reduced = PCA(n_components=pca_dim, random_state=42).fit_transform(x_scaled)

    rows: list[dict[str, Any]] = []
    for _, r in metrics_df.iterrows():
        method = str(r["method"])
        k = int(r["k"])
        key = (method, k)
        if key not in labels_map:
            continue
        labels = labels_map[key].values
        stats = estimate_stability_by_subsampling(
            x=x_reduced,
            full_labels=labels,
            recluster_fn=_method_recluster_fn(k),
            runs=6,
            subsample_ratio=0.8,
            random_state=42,
        )
        rows.append(
            {
                "method": method,
                "modality_set": modality_set_name,
                "k": k,
                "resampling_runs": int(stats["resampling_runs"]),
                "mean_nmi": stats["mean_nmi"],
                "mean_ari": stats["mean_ari"],
                "consensus_stability": stats["consensus_stability"],
                "notes": "round1_subsampling_stability",
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(out_table, index=False, encoding="utf-8")
    return out


def _load_clinical_tables(project_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    clinical = pd.read_csv(project_root / "data" / "raw" / "clinical" / "TCGA.STAD.sampleMap_STAD_clinicalMatrix", sep="\t", dtype=str)
    survival = pd.read_csv(project_root / "data" / "raw" / "clinical" / "survival_STAD_survival.txt", sep="\t", dtype=str)
    return clinical, survival


def _run_clinical_round1(
    metrics_df: pd.DataFrame,
    labels_map: dict[tuple[str, int], pd.Series],
    modality_set_name: str,
    figures_dir: Path,
    out_table: Path,
    clinical_df: pd.DataFrame,
    survival_df: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    # best-k per method by silhouette
    best_per_method = (
        metrics_df.sort_values(["method", "silhouette"], ascending=[True, False])
        .groupby("method", as_index=False)
        .head(1)
    )

    for _, rec in best_per_method.iterrows():
        method = str(rec["method"])
        k = int(rec["k"])
        labels = labels_map[(method, k)]

        for endpoint in ["OS", "PFI"]:
            stats = evaluate_clinical_association(labels, clinical_df, survival_df, endpoint=endpoint)
            rows.append(
                {
                    "method": method,
                    "modality_set": modality_set_name,
                    "k": k,
                    "endpoint": endpoint,
                    "logrank_p": stats["logrank_p"],
                    "cox_hr_summary": stats["cox_hr_summary"],
                    "clinical_variable": "survival_endpoint",
                    "association_p": stats["logrank_p"],
                    "notes": stats["notes"],
                }
            )
            plot_km_curve(
                labels,
                survival_df,
                endpoint=endpoint,
                out_path=figures_dir / f"km_{method}_{modality_set_name}_bestk_{endpoint.lower()}.png",
            )

        var_rows = evaluate_clinical_variables(labels, clinical_df)
        for vr in var_rows:
            rows.append(
                {
                    "method": method,
                    "modality_set": modality_set_name,
                    "k": k,
                    "endpoint": "clinical_variable",
                    "logrank_p": np.nan,
                    "cox_hr_summary": "n/a",
                    "clinical_variable": vr["clinical_variable"],
                    "association_p": vr["association_p"],
                    "notes": vr["notes"],
                }
            )

    out = pd.DataFrame(rows)
    out.to_csv(out_table, index=False, encoding="utf-8")
    plot_clinical_association_heatmap(out[out["endpoint"] == "clinical_variable"], figures_dir / "clinical_association_heatmap_round1.png")
    return out


def _run_proposed_readiness(
    project_root: Path,
    modalities: list[str],
    metrics_df: pd.DataFrame,
    readiness_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    tables_dir = project_root / "results" / "tables"

    q_scores: dict[str, float] = {}
    c_scores: dict[str, float] = {}
    s_scores: dict[str, float] = {}
    m_scores: dict[str, float] = {}

    # Build lightweight score components from available round1 evidence.
    for m in modalities:
        cov = float(readiness_df.loc[readiness_df["modality"] == m, "sample_count"].iloc[0])
        cov_ratio = cov / 364.0

        feat = float(readiness_df.loc[readiness_df["modality"] == m, "feature_count"].iloc[0])
        q = float(1.0 / np.log1p(max(feat, 1.0)))
        c = float(1.0 - min(1.0, metrics_df["logrank_p"].dropna().min() if "logrank_p" in metrics_df.columns else 0.5))
        s = float(cov_ratio)
        miss_raw = readiness_df.loc[readiness_df["modality"] == m, "missingness_level"].iloc[0]
        miss = 0.0
        if isinstance(miss_raw, str) and miss_raw.endswith("%"):
            miss = float(miss_raw[:-1]) / 100.0
        q_scores[m] = q
        c_scores[m] = c
        s_scores[m] = s
        m_scores[m] = miss

    preview_df = build_weight_input_preview(modalities, q_scores, c_scores, s_scores, m_scores)
    score_table = preview_df.set_index("modality")[["Q_m", "C_m", "S_m", "M_m"]]
    est_df = dry_run_weight_estimation(score_table, {"alpha": 1.0, "beta": 1.0, "delta": 1.0, "gamma": 1.0})
    est_df["notes"] = "round0_softmax_estimate_not_final"

    preview_df.to_csv(tables_dir / "proposed_weight_inputs_preview.csv", index=False, encoding="utf-8")
    est_df.to_csv(tables_dir / "proposed_weight_estimates_round0.csv", index=False, encoding="utf-8")
    return preview_df, est_df


def run_phase11_15(project_root: Path) -> None:
    paths_cfg, default_cfg = load_project_config(project_root)
    tables_dir = Path(paths_cfg["results"]["tables"])
    figures_dir = Path(paths_cfg["results"]["figures"])
    ensure_dir(tables_dir)
    ensure_dir(figures_dir)

    # Phase 11
    readiness_df, baseline_modalities, plan = build_modeling_readiness_summary(project_root)
    modality_set_name = "main5" if plan == "main5" else "main4"

    # similarity build on selected baseline modalities
    sim_out = build_similarity_round1(project_root, modalities=baseline_modalities)

    aligned_mats = {
        m: pd.read_csv(Path(sim_out["aligned_dir"]) / f"{m}_aligned.csv", index_col=0)
        for m in baseline_modalities
    }

    # Phase 12 baselines
    k_values = list(default_cfg.get("k_range", [2, 3, 4, 5, 6]))
    early = evaluate_early_fusion_round1(
        modality_tables=aligned_mats,
        k_values=k_values,
        algorithms=["kmeans", "spectral"],
        out_tables_dir=tables_dir,
        out_figures_dir=figures_dir,
        modality_set_name=modality_set_name,
        random_state=int(default_cfg.get("random_seed", 42)),
    )

    # Build one combined aligned matrix for SNF metric calculations
    combined_x = pd.concat([aligned_mats[m] for m in baseline_modalities], axis=1)
    snf = evaluate_equal_weight_fusion_round1(
        affinity_mats=sim_out["affinities"],
        aligned_feature_matrix=combined_x,
        k_values=k_values,
        out_tables_dir=tables_dir,
        out_figures_dir=figures_dir,
        modality_set_name=modality_set_name,
        random_state=int(default_cfg.get("random_seed", 42)),
    )

    # Merge baseline metrics/sizes
    metrics_df = pd.concat([early["metrics"], snf["metrics"]], ignore_index=True)
    sizes_df = pd.concat([early["cluster_sizes"], snf["cluster_sizes"]], ignore_index=True)
    metrics_df.to_csv(tables_dir / "baseline_round1_metrics.csv", index=False, encoding="utf-8")
    sizes_df.to_csv(tables_dir / "baseline_round1_cluster_sizes.csv", index=False, encoding="utf-8")
    snf["fused_summary"].to_csv(tables_dir / "fused_similarity_summary.csv", index=False, encoding="utf-8")

    _plot_baseline_comparison(metrics_df, figures_dir / "baseline_metric_comparison.png")
    _plot_cluster_size(sizes_df, figures_dir / "cluster_size_comparison.png")

    # Phase 13 stability
    labels_map: dict[tuple[str, int], pd.Series] = {}
    for (algo, k), s in early["labels"].items():
        labels_map[(f"early_fusion_{algo}", k)] = s
    for k, s in snf["labels"].items():
        labels_map[("equal_weight_fusion", k)] = s

    stability_df = _run_stability_round1(
        metrics_df=metrics_df,
        labels_map=labels_map,
        x_scaled=early["x_scaled"],
        modality_set_name=modality_set_name,
        out_table=tables_dir / "cluster_stability_round1.csv",
    )

    # Inject consensus to metrics table
    if not stability_df.empty:
        metrics_df = metrics_df.merge(
            stability_df[["method", "k", "consensus_stability"]],
            on=["method", "k"],
            how="left",
            suffixes=("", "_stability"),
        )
        metrics_df["consensus_stability"] = metrics_df["consensus_stability_stability"].combine_first(metrics_df["consensus_stability"])
        metrics_df = metrics_df.drop(columns=["consensus_stability_stability"])
        metrics_df.to_csv(tables_dir / "baseline_round1_metrics.csv", index=False, encoding="utf-8")

    # Clinical associations
    clinical_df, survival_df = _load_clinical_tables(project_root)
    _run_clinical_round1(
        metrics_df=metrics_df,
        labels_map=labels_map,
        modality_set_name=modality_set_name,
        figures_dir=figures_dir,
        out_table=tables_dir / "clinical_association_round1.csv",
        clinical_df=clinical_df,
        survival_df=survival_df,
    )

    # Phase 14 proposed readiness preview
    _run_proposed_readiness(project_root, baseline_modalities, metrics_df, readiness_df)


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[2]
    run_phase11_15(root)
    print("Phase 11-15 run completed.")
