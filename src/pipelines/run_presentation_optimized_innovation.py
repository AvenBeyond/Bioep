"""Build a presentation-focused optimized innovation branch.

This script does not replace the final baseline. It creates an independent
optimized shared-embedding v2 branch and compares it with the selected
traditional final baseline on report-friendly metrics.
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
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    normalized_mutual_info_score,
    silhouette_score,
)
from sklearn.preprocessing import StandardScaler


MAIN5_MODALITIES = ["mutation", "cnv", "methylation", "rna", "mirna"]


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _normalize_sample_id(sample_id: str) -> str:
    sid = str(sample_id).strip().upper()
    return sid[:16] if len(sid) >= 16 else sid


def _read_interim(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).map(_normalize_sample_id)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def _load_main5_matrices(root: Path) -> dict[str, pd.DataFrame]:
    interim = root / "data" / "interim"
    mats = {
        "mutation": _read_interim(interim / "mutation_round1.csv"),
        "cnv": _read_interim(interim / "cnv_round1.csv"),
        "methylation": _read_interim(interim / "methylation_round2_modeling.csv"),
        "rna": _read_interim(interim / "rna_round1.csv"),
        "mirna": _read_interim(interim / "mirna_round1.csv"),
    }
    cohort = sorted(set.intersection(*[set(df.index) for df in mats.values()]))
    return {m: df.loc[cohort].copy() for m, df in mats.items()}


def _load_survival(root: Path) -> pd.DataFrame:
    survival = pd.read_csv(root / "data" / "raw" / "clinical" / "survival_STAD_survival.txt", sep="\t", dtype=str)
    survival["normalized_sample_id"] = survival["sample"].astype(str).map(_normalize_sample_id)
    return survival


def _build_full_feature_matrix(mats: dict[str, pd.DataFrame]) -> pd.DataFrame:
    aligned = []
    for modality in MAIN5_MODALITIES:
        x = mats[modality].copy()
        x.columns = [f"{modality}::{c}" for c in x.columns]
        aligned.append(x)
    return pd.concat(aligned, axis=1)


def _build_optimized_shared_embedding(
    mats: dict[str, pd.DataFrame],
    per_modality_dim: int = 40,
    shared_dim: int = 50,
    random_state: int = 42,
) -> tuple[np.ndarray, dict[str, int]]:
    blocks: list[np.ndarray] = []
    dims: dict[str, int] = {}
    for modality in MAIN5_MODALITIES:
        x = StandardScaler().fit_transform(mats[modality].values)
        dim = min(per_modality_dim, x.shape[1], x.shape[0] - 1)
        dim = max(2, dim)
        dims[modality] = int(dim)
        blocks.append(PCA(n_components=dim, random_state=random_state).fit_transform(x))

    concat = np.concatenate(blocks, axis=1)
    out_dim = min(shared_dim, concat.shape[1], concat.shape[0] - 1)
    z = PCA(n_components=max(2, out_dim), random_state=random_state).fit_transform(concat)
    dims["shared_embedding"] = int(z.shape[1])
    return z, dims


def _cluster_balance(labels: np.ndarray) -> float:
    counts = np.bincount(labels.astype(int))
    counts = counts[counts > 0]
    return float(counts.min() / counts.max()) if counts.size else 0.0


def _clinical_p(labels: pd.Series, survival_df: pd.DataFrame, endpoint: str) -> float:
    time_col = f"{endpoint}.time"
    event_col = endpoint
    if time_col not in survival_df.columns or event_col not in survival_df.columns:
        return np.nan

    merged = labels.rename("cluster").to_frame()
    merged["normalized_sample_id"] = merged.index
    merged = merged.merge(survival_df, on="normalized_sample_id", how="inner")
    merged[time_col] = pd.to_numeric(merged[time_col], errors="coerce")
    merged[event_col] = pd.to_numeric(merged[event_col], errors="coerce")
    merged = merged.dropna(subset=[time_col, event_col, "cluster"]).copy()
    merged = merged[merged[time_col] >= 0]
    if merged.empty or merged["cluster"].nunique() < 2 or merged[event_col].sum() <= 1:
        return np.nan

    stat = multivariate_logrank_test(merged[time_col], merged["cluster"], merged[event_col])
    return float(stat.p_value)


def _stability_for_embedding(
    embedding: np.ndarray,
    labels: np.ndarray,
    k: int,
    runs: int = 30,
    subsample_ratio: float = 0.8,
    random_state: int = 42,
) -> tuple[float, float, float]:
    rng = np.random.default_rng(random_state)
    n = embedding.shape[0]
    sub_n = max(2, int(n * subsample_ratio))
    nmis: list[float] = []
    aris: list[float] = []
    for _ in range(runs):
        idx = np.sort(rng.choice(n, size=sub_n, replace=False))
        sub_labels = KMeans(n_clusters=k, n_init=30, random_state=random_state).fit_predict(embedding[idx])
        nmis.append(float(normalized_mutual_info_score(labels[idx], sub_labels)))
        aris.append(float(adjusted_rand_score(labels[idx], sub_labels)))
    nmi = float(np.mean(nmis))
    ari = float(np.mean(aris))
    return nmi, ari, float((nmi + ari) / 2.0)


def _clinical_score(os_p: float, pfi_p: float) -> float:
    vals = [p for p in [os_p, pfi_p] if pd.notna(p) and np.isfinite(p)]
    if not vals:
        return np.nan
    return float(np.mean([-np.log10(max(float(p), 1e-12)) for p in vals]))


def _make_metric_row(
    method: str,
    labels: pd.Series,
    metric_embedding: np.ndarray,
    full_scaled: np.ndarray,
    survival_df: pd.DataFrame,
    stability_embedding: np.ndarray,
    notes: str,
) -> dict[str, Any]:
    labels_arr = labels.values.astype(int)
    nmi, ari, cons = _stability_for_embedding(
        embedding=stability_embedding,
        labels=labels_arr,
        k=int(labels.nunique()),
    )
    os_p = _clinical_p(labels, survival_df, "OS")
    pfi_p = _clinical_p(labels, survival_df, "PFI")
    dss_p = _clinical_p(labels, survival_df, "DSS")
    dfi_p = _clinical_p(labels, survival_df, "DFI")

    return {
        "method": method,
        "k": int(labels.nunique()),
        "n_samples": int(labels.shape[0]),
        "full_feature_silhouette": float(silhouette_score(full_scaled, labels_arr)),
        "method_space_silhouette": float(silhouette_score(metric_embedding, labels_arr)),
        "calinski_harabasz_full": float(calinski_harabasz_score(full_scaled, labels_arr)),
        "davies_bouldin_full": float(davies_bouldin_score(full_scaled, labels_arr)),
        "stability_nmi": nmi,
        "stability_ari": ari,
        "consensus_stability": cons,
        "cluster_balance": _cluster_balance(labels_arr),
        "os_logrank_p": os_p,
        "pfi_logrank_p": pfi_p,
        "dss_logrank_p": dss_p,
        "dfi_logrank_p": dfi_p,
        "clinical_score_os_pfi": _clinical_score(os_p, pfi_p),
        "notes": notes,
    }


def _plot_metric_comparison(comparison: pd.DataFrame, out_path: Path) -> None:
    plot_df = comparison[
        ["method", "full_feature_silhouette", "consensus_stability", "cluster_balance", "clinical_score_os_pfi"]
    ].copy()
    plot_df = plot_df.melt(id_vars="method", var_name="metric", value_name="value")

    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10.5, 5.2))
    sns.barplot(data=plot_df, x="metric", y="value", hue="method")
    plt.title("Optimized Innovation vs Traditional Baseline")
    plt.ylabel("raw metric value")
    plt.xticks(rotation=12, ha="right")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=240)
    plt.close()


def _plot_embedding(z: np.ndarray, labels: pd.Series, out_path: Path) -> None:
    emb2 = PCA(n_components=2, random_state=42).fit_transform(z)
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(7, 6))
    plt.scatter(emb2[:, 0], emb2[:, 1], c=labels.values, cmap="tab10", s=24, alpha=0.9)
    plt.title("Optimized Shared Embedding v2")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=240)
    plt.close()


def _plot_km(labels: pd.Series, survival_df: pd.DataFrame, endpoint: str, out_path: Path) -> None:
    time_col = f"{endpoint}.time"
    event_col = endpoint
    merged = labels.rename("cluster").to_frame()
    merged["normalized_sample_id"] = merged.index
    merged = merged.merge(survival_df, on="normalized_sample_id", how="inner")
    merged[time_col] = pd.to_numeric(merged[time_col], errors="coerce")
    merged[event_col] = pd.to_numeric(merged[event_col], errors="coerce")
    merged = merged.dropna(subset=[time_col, event_col, "cluster"]).copy()
    merged = merged[merged[time_col] >= 0]
    if merged.empty:
        return

    kmf = KaplanMeierFitter()
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(7, 5))
    for cid in sorted(merged["cluster"].astype(int).unique()):
        part = merged[merged["cluster"].astype(int) == cid]
        kmf.fit(part[time_col], event_observed=part[event_col], label=f"cluster {cid}")
        kmf.plot(ci_show=False)
    plt.title(f"Optimized Shared Embedding KM ({endpoint})")
    plt.xlabel("Time")
    plt.ylabel("Survival probability")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=240)
    plt.close()


def _build_win_table(comparison: pd.DataFrame) -> pd.DataFrame:
    base = comparison[comparison["method"] == "traditional_early_fusion_kmeans"].iloc[0]
    opt = comparison[comparison["method"] == "optimized_shared_embedding_v2"].iloc[0]
    specs = [
        ("full_feature_silhouette", "higher"),
        ("method_space_silhouette", "higher"),
        ("consensus_stability", "higher"),
        ("cluster_balance", "higher"),
        ("os_logrank_p", "lower"),
        ("pfi_logrank_p", "lower"),
        ("clinical_score_os_pfi", "higher"),
    ]
    rows: list[dict[str, Any]] = []
    for metric, direction in specs:
        b = float(base[metric])
        o = float(opt[metric])
        if direction == "higher":
            winner = "optimized_shared_embedding_v2" if o > b else "traditional_early_fusion_kmeans"
        else:
            winner = "optimized_shared_embedding_v2" if o < b else "traditional_early_fusion_kmeans"
        rows.append(
            {
                "metric": metric,
                "direction": direction,
                "traditional_baseline": b,
                "optimized_innovation": o,
                "winner": winner,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    root = _project_root()
    tables = root / "results" / "tables"
    figures = root / "results" / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    mats = _load_main5_matrices(root)
    full = _build_full_feature_matrix(mats)
    full_scaled = StandardScaler().fit_transform(full.values)
    full_pca = PCA(n_components=min(50, full_scaled.shape[1], full_scaled.shape[0] - 1), random_state=42).fit_transform(full_scaled)
    survival = _load_survival(root)

    z, dims = _build_optimized_shared_embedding(mats, per_modality_dim=40, shared_dim=50, random_state=42)
    opt_labels_arr = KMeans(n_clusters=2, n_init=80, random_state=42).fit_predict(z)
    opt_labels = pd.Series(opt_labels_arr, index=full.index, name="cluster")

    base_path = tables / "baseline_main5_labels_early_fusion_kmeans_k2.csv"
    base_labels = pd.read_csv(base_path, index_col=0).iloc[:, 0]
    base_labels.index = base_labels.index.astype(str).map(_normalize_sample_id)
    base_labels = base_labels.loc[full.index].astype(int)
    base_labels.name = "cluster"

    comparison = pd.DataFrame(
        [
            _make_metric_row(
                method="traditional_early_fusion_kmeans",
                labels=base_labels,
                metric_embedding=full_scaled,
                full_scaled=full_scaled,
                survival_df=survival,
                stability_embedding=full_pca,
                notes="selected_final_traditional_baseline",
            ),
            _make_metric_row(
                method="optimized_shared_embedding_v2",
                labels=opt_labels,
                metric_embedding=z,
                full_scaled=full_scaled,
                survival_df=survival,
                stability_embedding=z,
                notes=f"per_modality_pca=40;shared_pca=50;dims={dims}",
            ),
        ]
    )

    opt_labels.to_frame().to_csv(tables / "optimized_shared_embedding_v2_labels_k2.csv", encoding="utf-8")
    comparison.to_csv(tables / "optimized_innovation_comparison.csv", index=False, encoding="utf-8")
    _build_win_table(comparison).to_csv(tables / "optimized_innovation_metric_wins.csv", index=False, encoding="utf-8")

    _plot_metric_comparison(comparison, figures / "optimized_innovation_comparison.png")
    _plot_embedding(z, opt_labels, figures / "optimized_shared_embedding_v2_pca.png")
    _plot_km(opt_labels, survival, "OS", figures / "optimized_shared_embedding_v2_km_os.png")
    _plot_km(opt_labels, survival, "PFI", figures / "optimized_shared_embedding_v2_km_pfi.png")

    print("Optimized innovation artifacts generated.")
    print(comparison.to_string(index=False))


if __name__ == "__main__":
    main()
