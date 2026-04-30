"""Graph connectivity-aware tuning utilities (improvement round)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score
from sklearn.metrics.pairwise import euclidean_distances


def _knn_affinity(
    x: np.ndarray,
    n_neighbors: int,
    mutual_knn: bool,
    local_scaling: bool,
    threshold_q: float | None,
) -> np.ndarray:
    d = euclidean_distances(x)
    n = d.shape[0]
    idx = np.argsort(d, axis=1)[:, 1 : n_neighbors + 1]

    w = np.zeros((n, n), dtype=float)
    if local_scaling:
        # Zelnik-Manor local scaling proxy
        sigma = np.take_along_axis(d, idx[:, -1][:, None], axis=1).reshape(-1)
        sigma = np.clip(sigma, 1e-6, None)
    else:
        sigma = np.full(n, np.std(d[np.triu_indices(n, k=1)]) + 1e-6)

    for i in range(n):
        for j in idx[i]:
            if local_scaling:
                val = np.exp(-(d[i, j] ** 2) / (sigma[i] * sigma[j]))
            else:
                val = np.exp(-(d[i, j] ** 2) / (2 * sigma[i] ** 2))
            w[i, j] = val

    if mutual_knn:
        w = np.minimum(w, w.T)
    else:
        w = np.maximum(w, w.T)

    if threshold_q is not None:
        nz = w[w > 0]
        if nz.size > 0:
            th = np.quantile(nz, threshold_q)
            w[w < th] = 0.0

    np.fill_diagonal(w, 1.0)
    return w


def run_graph_connectivity_tuning(
    x_df: pd.DataFrame,
    k_cluster: int,
    out_tables_dir: Path,
    out_figures_dir: Path,
    random_state: int = 42,
) -> dict[str, Any]:
    out_tables_dir.mkdir(parents=True, exist_ok=True)
    out_figures_dir.mkdir(parents=True, exist_ok=True)

    x = x_df.values
    configs = []
    for n_neighbors in [5, 10, 15, 20]:
        for mutual in [False, True]:
            for local in [False, True]:
                for q in [None, 0.2]:
                    configs.append((n_neighbors, mutual, local, q))

    rows: list[dict[str, Any]] = []
    best_tuple = None
    best_score = -np.inf

    for n_neighbors, mutual, local, q in configs:
        a = _knn_affinity(x, n_neighbors=n_neighbors, mutual_knn=mutual, local_scaling=local, threshold_q=q)
        g = csr_matrix((a > 0).astype(int))
        n_comp, labels_comp = connected_components(g, directed=False)
        comp_sizes = pd.Series(labels_comp).value_counts().values
        isolated = int((comp_sizes == 1).sum())
        density = float((a > 0).sum() / (a.shape[0] * a.shape[1]))

        labels = SpectralClustering(n_clusters=k_cluster, affinity="precomputed", random_state=random_state).fit_predict(a)
        if np.unique(labels).size >= 2:
            sil = float(silhouette_score(x, labels))
        else:
            sil = np.nan

        # Composite objective: fewer components and higher silhouette.
        obj = (0.7 * (sil if np.isfinite(sil) else -1.0)) - 0.2 * (n_comp - 1) - 0.1 * isolated
        rows.append(
            {
                "n_neighbors": n_neighbors,
                "mutual_knn": bool(mutual),
                "local_scaling": bool(local),
                "threshold_strategy": "quantile_0.2" if q is not None else "none",
                "n_components": int(n_comp),
                "isolated_nodes": int(isolated),
                "graph_density": density,
                "silhouette": sil,
                "objective_score": obj,
                "notes": "graph_connectivity_tuning",
            }
        )
        if obj > best_score:
            best_score = obj
            best_tuple = (n_neighbors, mutual, local, q)

    res_df = pd.DataFrame(rows).sort_values("objective_score", ascending=False)
    res_df.to_csv(out_tables_dir / "graph_connectivity_tuning_results.csv", index=False, encoding="utf-8")
    res_df.head(1).to_csv(out_tables_dir / "graph_connectivity_best_settings.csv", index=False, encoding="utf-8")

    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(11, 5))
    plot_df = res_df.copy().head(12)
    plot_df["cfg"] = plot_df.apply(
        lambda r: f"kNN={int(r['n_neighbors'])}|M={int(bool(r['mutual_knn']))}|L={int(bool(r['local_scaling']))}|T={r['threshold_strategy']}",
        axis=1,
    )
    sns.barplot(data=plot_df, x="cfg", y="objective_score", color="#4361ee")
    plt.xticks(rotation=45, ha="right")
    plt.title("Graph Connectivity Tuning (Top Configurations)")
    plt.tight_layout()
    plt.savefig(out_figures_dir / "graph_connectivity_tuning.png", dpi=220)
    plt.close()

    plt.figure(figsize=(10, 4))
    sns.scatterplot(data=res_df, x="n_components", y="silhouette", hue="mutual_knn", style="local_scaling", s=90)
    plt.title("Graph Components vs Clustering Quality")
    plt.tight_layout()
    plt.savefig(out_figures_dir / "graph_component_summary_v2.png", dpi=220)
    plt.close()

    return {"results": res_df, "best": res_df.head(1)}
