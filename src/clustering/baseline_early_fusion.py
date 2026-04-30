"""Baseline 1: executable early-fusion clustering for round1.

Input:
- aligned modality matrices (sample x feature)

Output:
- labels per method/k
- metrics and cluster-size summaries
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score
from sklearn.preprocessing import StandardScaler


def _cluster_balance_score(labels: np.ndarray) -> float:
    counts = pd.Series(labels).value_counts().values
    if counts.size == 0:
        return 0.0
    return float(counts.min() / counts.max())


def run_early_fusion_clustering(
    modality_tables: dict[str, pd.DataFrame],
    k: int = 3,
    algorithm: str = "kmeans",
    random_state: int = 42,
) -> pd.Series:
    modalities = sorted(modality_tables.keys())
    cohort = sorted(set.intersection(*[set(modality_tables[m].index) for m in modalities]))
    x = pd.concat([modality_tables[m].loc[cohort] for m in modalities], axis=1)
    x_scaled = StandardScaler().fit_transform(x.values)

    if algorithm == "kmeans":
        model = KMeans(n_clusters=k, n_init=30, random_state=random_state)
        labels = model.fit_predict(x_scaled)
    elif algorithm == "spectral":
        model = SpectralClustering(n_clusters=k, affinity="nearest_neighbors", random_state=random_state)
        labels = model.fit_predict(x_scaled)
    else:
        raise ValueError(f"Unsupported algorithm: {algorithm}")

    return pd.Series(labels, index=cohort, name="cluster")


def evaluate_early_fusion_round1(
    modality_tables: dict[str, pd.DataFrame],
    k_values: list[int],
    algorithms: list[str],
    out_tables_dir: Path,
    out_figures_dir: Path,
    modality_set_name: str,
    random_state: int = 42,
) -> dict[str, Any]:
    out_tables_dir.mkdir(parents=True, exist_ok=True)
    out_figures_dir.mkdir(parents=True, exist_ok=True)

    modalities = sorted(modality_tables.keys())
    cohort = sorted(set.intersection(*[set(modality_tables[m].index) for m in modalities]))
    x = pd.concat([modality_tables[m].loc[cohort] for m in modalities], axis=1)
    x_scaled = StandardScaler().fit_transform(x.values)

    metric_rows: list[dict[str, Any]] = []
    cluster_size_rows: list[dict[str, Any]] = []
    labels_out: dict[tuple[str, int], pd.Series] = {}

    for algorithm in algorithms:
        for k in k_values:
            labels = run_early_fusion_clustering(
                modality_tables,
                k=k,
                algorithm=algorithm,
                random_state=random_state,
            )

            arr = labels.values
            metric_rows.append(
                {
                    "method": f"early_fusion_{algorithm}",
                    "modality_set": modality_set_name,
                    "k": k,
                    "n_samples": len(labels),
                    "silhouette": float(silhouette_score(x_scaled, arr)),
                    "calinski_harabasz": float(calinski_harabasz_score(x_scaled, arr)),
                    "davies_bouldin": float(davies_bouldin_score(x_scaled, arr)),
                    "cluster_balance_score": _cluster_balance_score(arr),
                    "consensus_stability": np.nan,
                    "notes": "round1_real_run",
                }
            )

            for cid, csize in labels.value_counts().sort_index().items():
                cluster_size_rows.append(
                    {
                        "method": f"early_fusion_{algorithm}",
                        "modality_set": modality_set_name,
                        "k": k,
                        "cluster_id": int(cid),
                        "cluster_size": int(csize),
                    }
                )

            labels_out[(algorithm, k)] = labels
            labels.to_frame().to_csv(
                out_tables_dir / f"baseline_round1_labels_early_fusion_{algorithm}_{modality_set_name}_k{k}.csv",
                encoding="utf-8",
            )

        # Embedding plot for best-k by silhouette
        algo_metrics = [m for m in metric_rows if m["method"] == f"early_fusion_{algorithm}"]
        best_k = sorted(algo_metrics, key=lambda r: r["silhouette"], reverse=True)[0]["k"]
        best_labels = labels_out[(algorithm, best_k)].values

        if x_scaled.shape[0] > 300:
            emb = PCA(n_components=2, random_state=random_state).fit_transform(x_scaled)
        else:
            emb = TSNE(n_components=2, init="pca", learning_rate="auto", random_state=random_state).fit_transform(x_scaled)

        plt.figure(figsize=(7, 6))
        plt.scatter(emb[:, 0], emb[:, 1], c=best_labels, cmap="tab10", s=20)
        plt.title(f"Embedding (round1): early_fusion_{algorithm}, best k={best_k}")
        plt.xlabel("dim1")
        plt.ylabel("dim2")
        plt.tight_layout()
        plt.savefig(
            out_figures_dir / f"embedding_umap_early_fusion_{algorithm}_{modality_set_name}_bestk.png",
            dpi=220,
        )
        plt.close()

    return {
        "metrics": pd.DataFrame(metric_rows),
        "cluster_sizes": pd.DataFrame(cluster_size_rows),
        "labels": labels_out,
        "x_scaled": x_scaled,
    }
