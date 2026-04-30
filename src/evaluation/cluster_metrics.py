"""Clustering metrics and round-1 stability utilities."""

from __future__ import annotations

from typing import Any

import numpy as np
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    normalized_mutual_info_score,
    silhouette_score,
)


def compute_basic_cluster_metrics(x: np.ndarray, labels: np.ndarray) -> dict[str, Any]:
    return {
        "silhouette": float(silhouette_score(x, labels)),
        "ch": float(calinski_harabasz_score(x, labels)),
        "dbi": float(davies_bouldin_score(x, labels)),
    }


def cluster_balance_score(labels: np.ndarray) -> float:
    counts = np.bincount(labels)
    counts = counts[counts > 0]
    if counts.size == 0:
        return 0.0
    return float(counts.min() / counts.max())


def estimate_stability_by_subsampling(
    x: np.ndarray,
    full_labels: np.ndarray,
    recluster_fn,
    runs: int = 20,
    subsample_ratio: float = 0.8,
    random_state: int = 42,
) -> dict[str, float]:
    rng = np.random.default_rng(random_state)
    n = x.shape[0]
    sub_n = max(2, int(n * subsample_ratio))

    nmis: list[float] = []
    aris: list[float] = []
    for _ in range(runs):
        idx = np.sort(rng.choice(n, size=sub_n, replace=False))
        x_sub = x[idx, :]
        baseline_sub = full_labels[idx]
        new_sub = recluster_fn(x_sub)
        nmis.append(float(normalized_mutual_info_score(baseline_sub, new_sub)))
        aris.append(float(adjusted_rand_score(baseline_sub, new_sub)))

    return {
        "resampling_runs": float(runs),
        "mean_nmi": float(np.mean(nmis)) if nmis else np.nan,
        "mean_ari": float(np.mean(aris)) if aris else np.nan,
        "consensus_stability": float((np.mean(nmis) + np.mean(aris)) / 2.0) if nmis and aris else np.nan,
    }
