"""Lightweight shared-embedding prototype (exploratory)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score
from sklearn.preprocessing import StandardScaler


def run_shared_embedding_prototype(
    modality_tables: dict[str, pd.DataFrame],
    k_values: list[int],
    out_tables_dir: Path,
    out_figures_dir: Path,
    random_state: int = 42,
) -> dict[str, Any]:
    out_tables_dir.mkdir(parents=True, exist_ok=True)
    out_figures_dir.mkdir(parents=True, exist_ok=True)

    modalities = sorted(modality_tables.keys())
    cohort = sorted(set.intersection(*[set(modality_tables[m].index) for m in modalities]))

    reduced_blocks = []
    for m in modalities:
        x = modality_tables[m].loc[cohort].values
        x = StandardScaler().fit_transform(x)
        dim = min(20, x.shape[1], x.shape[0] - 1)
        reduced = PCA(n_components=max(2, dim), random_state=random_state).fit_transform(x)
        reduced_blocks.append(reduced)

    concat = np.concatenate(reduced_blocks, axis=1)
    shared_dim = min(32, concat.shape[1], concat.shape[0] - 1)
    z = PCA(n_components=max(2, shared_dim), random_state=random_state).fit_transform(concat)

    rows = []
    labels_map: dict[int, pd.Series] = {}
    for k in k_values:
        labels = KMeans(n_clusters=k, n_init=30, random_state=random_state).fit_predict(z)
        s = pd.Series(labels, index=cohort, name="cluster")
        labels_map[k] = s
        s.to_frame().to_csv(out_tables_dir / f"shared_embedding_labels_k{k}.csv", encoding="utf-8")

        uniq = np.unique(labels)
        if uniq.size >= 2:
            sil = float(silhouette_score(z, labels))
            ch = float(calinski_harabasz_score(z, labels))
            dbi = float(davies_bouldin_score(z, labels))
        else:
            sil, ch, dbi = np.nan, np.nan, np.nan
        counts = s.value_counts().values
        bal = float(counts.min() / counts.max()) if counts.size > 0 else 0.0
        rows.append(
            {
                "method": "shared_embedding_prototype",
                "variant": "pca_per_modality + shared_pca",
                "k": int(k),
                "silhouette": sil,
                "calinski_harabasz": ch,
                "davies_bouldin": dbi,
                "cluster_balance_score": bal,
                "notes": "exploratory_lightweight",
            }
        )

    metrics_df = pd.DataFrame(rows)
    metrics_df.to_csv(out_tables_dir / "shared_embedding_prototype_metrics.csv", index=False, encoding="utf-8")

    best_k = int(metrics_df.sort_values("silhouette", ascending=False).iloc[0]["k"])
    emb2 = PCA(n_components=2, random_state=random_state).fit_transform(z)
    plt.figure(figsize=(7, 6))
    plt.scatter(emb2[:, 0], emb2[:, 1], c=labels_map[best_k].values, cmap="tab10", s=20)
    plt.title(f"Shared Embedding Prototype (best k={best_k})")
    plt.tight_layout()
    plt.savefig(out_figures_dir / "shared_embedding_umap.png", dpi=220)
    plt.close()

    return {"metrics": metrics_df, "labels": labels_map, "shared_embedding": z}
