"""Consensus ensemble subtyping for improvement round."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import SpectralClustering
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score


def build_coassociation_matrix(label_series_list: list[pd.Series]) -> pd.DataFrame:
    if not label_series_list:
        raise ValueError("No label series provided for consensus ensemble.")
    samples = label_series_list[0].index.tolist()
    n = len(samples)
    mat = np.zeros((n, n), dtype=float)

    for s in label_series_list:
        arr = s.loc[samples].values
        eq = (arr[:, None] == arr[None, :]).astype(float)
        mat += eq
    mat /= len(label_series_list)

    out = pd.DataFrame(mat, index=samples, columns=samples)
    np.fill_diagonal(out.values, 1.0)
    return out


def run_consensus_ensemble(
    label_candidates: dict[str, pd.Series],
    combined_feature_matrix: pd.DataFrame,
    out_tables_dir: Path,
    out_figures_dir: Path,
    n_clusters: int,
    random_state: int = 42,
) -> dict[str, Any]:
    out_tables_dir.mkdir(parents=True, exist_ok=True)
    out_figures_dir.mkdir(parents=True, exist_ok=True)

    labels_list = list(label_candidates.values())
    coassoc = build_coassociation_matrix(labels_list)

    labels = SpectralClustering(n_clusters=n_clusters, affinity="precomputed", random_state=random_state).fit_predict(coassoc.values)
    consensus_labels = pd.Series(labels, index=coassoc.index, name="cluster")

    x = combined_feature_matrix.loc[coassoc.index].values
    uniq = np.unique(labels)
    if uniq.size >= 2:
        sil = float(silhouette_score(x, labels))
        ch = float(calinski_harabasz_score(x, labels))
        dbi = float(davies_bouldin_score(x, labels))
    else:
        sil, ch, dbi = np.nan, np.nan, np.nan
    counts = consensus_labels.value_counts().values
    balance = float(counts.min() / counts.max()) if counts.size > 0 else 0.0

    metrics_df = pd.DataFrame([
        {
            "method": "consensus_ensemble",
            "n_models": len(label_candidates),
            "k": int(n_clusters),
            "silhouette": sil,
            "calinski_harabasz": ch,
            "davies_bouldin": dbi,
            "cluster_balance_score": balance,
            "notes": "coassociation + spectral",
        }
    ])
    sizes_df = consensus_labels.value_counts().sort_index().rename_axis("cluster_id").reset_index(name="cluster_size")

    metrics_df.to_csv(out_tables_dir / "consensus_ensemble_metrics.csv", index=False, encoding="utf-8")
    sizes_df.to_csv(out_tables_dir / "consensus_ensemble_cluster_sizes.csv", index=False, encoding="utf-8")
    consensus_labels.to_frame().to_csv(out_tables_dir / "consensus_ensemble_labels.csv", encoding="utf-8")

    sns.set_theme(style="white")
    plt.figure(figsize=(7, 6))
    sns.heatmap(coassoc.iloc[:100, :100], cmap="viridis")
    plt.title("Consensus Co-association Heatmap (first 100)")
    plt.tight_layout()
    plt.savefig(out_figures_dir / "consensus_coassociation_heatmap.png", dpi=220)
    plt.close()

    comp_rows = []
    for name, s in label_candidates.items():
        arr = s.loc[coassoc.index].values
        if np.unique(arr).size >= 2:
            sv = float(silhouette_score(x, arr))
        else:
            sv = np.nan
        comp_rows.append({"model": name, "silhouette": sv})
    comp_rows.append({"model": "consensus_ensemble", "silhouette": sil})
    comp_df = pd.DataFrame(comp_rows)

    plt.figure(figsize=(9, 4))
    sns.barplot(data=comp_df, x="model", y="silhouette", color="#2a9d8f")
    plt.xticks(rotation=25, ha="right")
    plt.title("Consensus vs Single Models")
    plt.tight_layout()
    plt.savefig(out_figures_dir / "consensus_vs_single_model.png", dpi=220)
    plt.close()

    return {"metrics": metrics_df, "cluster_sizes": sizes_df, "labels": consensus_labels, "coassociation": coassoc}
