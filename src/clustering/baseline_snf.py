"""Baseline 2: equal-weight similarity network fusion (round1 executable)."""

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
from sklearn.decomposition import PCA
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score


def _cluster_balance_score(labels: np.ndarray) -> float:
    counts = pd.Series(labels).value_counts().values
    return float(counts.min() / counts.max()) if counts.size > 0 else 0.0


def run_equal_weight_fusion(similarity_mats: dict[str, pd.DataFrame]) -> pd.DataFrame:
    keys = list(similarity_mats.keys())
    fused = sum(similarity_mats[k] for k in keys) / len(keys)
    return fused


def evaluate_equal_weight_fusion_round1(
    affinity_mats: dict[str, pd.DataFrame],
    aligned_feature_matrix: pd.DataFrame,
    k_values: list[int],
    out_tables_dir: Path,
    out_figures_dir: Path,
    modality_set_name: str,
    random_state: int = 42,
) -> dict[str, Any]:
    out_tables_dir.mkdir(parents=True, exist_ok=True)
    out_figures_dir.mkdir(parents=True, exist_ok=True)

    fused = run_equal_weight_fusion(affinity_mats)
    fused = ((fused + fused.T) / 2.0).clip(lower=0)

    sns.set_theme(style="white")
    plt.figure(figsize=(7, 6))
    sns.heatmap(fused.iloc[:100, :100], cmap="YlGnBu")
    plt.title("Fused Network Heatmap (first 100 samples)")
    plt.xlabel("Samples")
    plt.ylabel("Samples")
    plt.tight_layout()
    plt.savefig(out_figures_dir / "fused_network_heatmap.png", dpi=220)
    plt.close()

    metric_rows: list[dict[str, Any]] = []
    cluster_size_rows: list[dict[str, Any]] = []
    labels_out: dict[int, pd.Series] = {}

    x_for_metrics = aligned_feature_matrix.values

    for k in k_values:
        model = SpectralClustering(n_clusters=k, affinity="precomputed", random_state=random_state)
        labels = model.fit_predict(fused.values)
        s = pd.Series(labels, index=fused.index, name="cluster")
        labels_out[k] = s
        s.to_frame().to_csv(
            out_tables_dir / f"baseline_round1_labels_equal_weight_fusion_{modality_set_name}_k{k}.csv",
            encoding="utf-8",
        )

        metric_rows.append(
            {
                "method": "equal_weight_fusion",
                "modality_set": modality_set_name,
                "k": k,
                "n_samples": int(len(s)),
                "silhouette": float(silhouette_score(x_for_metrics, labels)),
                "calinski_harabasz": float(calinski_harabasz_score(x_for_metrics, labels)),
                "davies_bouldin": float(davies_bouldin_score(x_for_metrics, labels)),
                "cluster_balance_score": _cluster_balance_score(labels),
                "consensus_stability": np.nan,
                "notes": "round1_real_run",
            }
        )

        for cid, csize in s.value_counts().sort_index().items():
            cluster_size_rows.append(
                {
                    "method": "equal_weight_fusion",
                    "modality_set": modality_set_name,
                    "k": k,
                    "cluster_id": int(cid),
                    "cluster_size": int(csize),
                }
            )

    best_k = sorted(metric_rows, key=lambda r: r["silhouette"], reverse=True)[0]["k"]
    best_labels = labels_out[best_k].values
    emb = PCA(n_components=2, random_state=random_state).fit_transform(x_for_metrics)
    plt.figure(figsize=(7, 6))
    plt.scatter(emb[:, 0], emb[:, 1], c=best_labels, cmap="tab10", s=20)
    plt.title(f"Embedding (round1): equal_weight_fusion, best k={best_k}")
    plt.xlabel("dim1")
    plt.ylabel("dim2")
    plt.tight_layout()
    plt.savefig(out_figures_dir / f"embedding_umap_equal_weight_fusion_{modality_set_name}_bestk.png", dpi=220)
    plt.close()

    fused_summary = pd.DataFrame(
        [
            {
                "method": "equal_weight_fusion",
                "modality_set": modality_set_name,
                "n_samples": fused.shape[0],
                "matrix_mean": float(fused.values.mean()),
                "matrix_std": float(fused.values.std()),
                "matrix_min": float(fused.values.min()),
                "matrix_max": float(fused.values.max()),
            }
        ]
    )

    return {
        "fused": fused,
        "metrics": pd.DataFrame(metric_rows),
        "cluster_sizes": pd.DataFrame(cluster_size_rows),
        "labels": labels_out,
        "fused_summary": fused_summary,
    }
