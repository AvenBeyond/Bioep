"""Proposed weighted-fusion implementation for formal main5 experiments.

Input:
- modality similarity matrices
- modality-level quality scores

Output:
- weighted fused matrix
- candidate K scores

Purpose:
- implement three innovations:
  1) quality + clinical weighted fusion
  2) weakly paired partial-fusion fallback
  3) multi-objective K selection

TODO:
- extend to more advanced partial-fusion graph optimization.
"""

from __future__ import annotations

import math
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


def compute_modality_weights(score_table: pd.DataFrame, alpha: float, beta: float, delta: float, gamma: float) -> pd.Series:
    # score_table columns: Q_m, C_m, S_m, M_m
    raw = alpha * score_table["Q_m"] + beta * score_table["C_m"] + delta * score_table["S_m"] - gamma * score_table["M_m"]
    raw = raw - raw.max()
    exps = raw.apply(math.exp)
    return exps / exps.sum()


def build_weight_input_preview(
    modalities: list[str],
    q_scores: dict[str, float],
    c_scores: dict[str, float],
    s_scores: dict[str, float],
    m_scores: dict[str, float],
) -> pd.DataFrame:
    rows = []
    for m in modalities:
        rows.append(
            {
                "modality": m,
                "Q_m": float(q_scores.get(m, np.nan)),
                "C_m": float(c_scores.get(m, np.nan)),
                "S_m": float(s_scores.get(m, np.nan)),
                "M_m": float(m_scores.get(m, np.nan)),
                "notes": "round0_input_preview",
            }
        )
    return pd.DataFrame(rows)


def partial_fusion(similarity_mats: dict[str, pd.DataFrame], weights: pd.Series) -> pd.DataFrame:
    # Placeholder: weighted sum on currently available modalities.
    fused = None
    for m, mat in similarity_mats.items():
        w = float(weights.get(m, 0.0))
        fused = mat * w if fused is None else fused + mat * w
    return fused


def score_k(silhouette: float, consensus: float, clinical_sep: float, balance: float) -> float:
    return 0.35 * silhouette + 0.25 * consensus + 0.20 * clinical_sep + 0.20 * balance


def select_best_k(metrics_by_k: dict[int, dict[str, float]]) -> tuple[int, dict[int, float]]:
    scores: dict[int, float] = {}
    for k, m in metrics_by_k.items():
        scores[k] = score_k(m["silhouette"], m["consensus"], m["clinical_sep"], m["balance"])
    best_k = max(scores, key=scores.get)
    return best_k, scores


def run_proposed_weighted_fusion(
    similarity_mats: dict[str, pd.DataFrame],
    score_table: pd.DataFrame,
    params: dict[str, Any],
) -> pd.DataFrame:
    weights = compute_modality_weights(
        score_table,
        alpha=float(params.get("alpha", 1.0)),
        beta=float(params.get("beta", 1.0)),
        delta=float(params.get("delta", 1.0)),
        gamma=float(params.get("gamma", 1.0)),
    )
    return partial_fusion(similarity_mats, weights)


def build_proposed_weight_components(
    modalities: list[str],
    q_scores: dict[str, float],
    c_scores: dict[str, float],
    s_scores: dict[str, float],
    m_scores: dict[str, float],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for m in modalities:
        rows.append(
            {
                "modality": m,
                "Q_m": float(q_scores.get(m, np.nan)),
                "C_m": float(c_scores.get(m, np.nan)),
                "S_m": float(s_scores.get(m, np.nan)),
                "M_m": float(m_scores.get(m, np.nan)),
            }
        )
    return pd.DataFrame(rows)


def estimate_final_weights(
    component_df: pd.DataFrame,
    alpha: float,
    beta: float,
    delta: float,
    gamma: float,
) -> pd.DataFrame:
    score_table = component_df.set_index("modality")[["Q_m", "C_m", "S_m", "M_m"]]
    weights = compute_modality_weights(score_table, alpha=alpha, beta=beta, delta=delta, gamma=gamma)
    out = component_df.copy()
    out["final_weight"] = out["modality"].map(weights.to_dict())
    return out.sort_values("final_weight", ascending=False).reset_index(drop=True)


def run_weighted_fusion_clustering(
    affinity_mats: dict[str, pd.DataFrame],
    weight_table: pd.DataFrame,
    combined_feature_matrix: pd.DataFrame,
    k_values: list[int],
    modality_set_name: str,
    out_tables_dir: Path,
    out_figures_dir: Path,
    label_prefix: str = "proposed_main5",
    random_state: int = 42,
) -> dict[str, Any]:
    out_tables_dir.mkdir(parents=True, exist_ok=True)
    out_figures_dir.mkdir(parents=True, exist_ok=True)

    weights = weight_table.set_index("modality")["final_weight"].to_dict()
    fused = None
    for m, mat in affinity_mats.items():
        if m not in weights:
            continue
        w = float(weights[m])
        fused = mat * w if fused is None else fused + mat * w
    if fused is None:
        raise ValueError("No overlap between affinity matrices and weight table modalities.")

    fused = ((fused + fused.T) / 2.0).clip(lower=0)

    sns.set_theme(style="white")
    plt.figure(figsize=(7, 6))
    sns.heatmap(fused.iloc[:100, :100], cmap="YlOrBr")
    plt.title("Proposed Weighted Fused Network (first 100 samples)")
    plt.tight_layout()
    plt.savefig(out_figures_dir / "proposed_fused_network_heatmap.png", dpi=220)
    plt.close()

    metric_rows: list[dict[str, Any]] = []
    cluster_size_rows: list[dict[str, Any]] = []
    labels_map: dict[int, pd.Series] = {}

    x = combined_feature_matrix.loc[fused.index].values
    for k in k_values:
        labels = SpectralClustering(n_clusters=k, affinity="precomputed", random_state=random_state).fit_predict(fused.values)
        s = pd.Series(labels, index=fused.index, name="cluster")
        labels_map[k] = s
        s.to_frame().to_csv(out_tables_dir / f"{label_prefix}_labels_k{k}.csv", encoding="utf-8")

        counts = s.value_counts().values
        balance = float(counts.min() / counts.max()) if counts.size > 0 else 0.0
        if np.unique(labels).size >= 2:
            sil = float(silhouette_score(x, labels))
            ch = float(calinski_harabasz_score(x, labels))
            dbi = float(davies_bouldin_score(x, labels))
        else:
            sil = np.nan
            ch = np.nan
            dbi = np.nan
        metric_rows.append(
            {
                "method": "proposed_weighted_fusion",
                "modality_set": modality_set_name,
                "k": int(k),
                "n_samples": int(len(s)),
                "silhouette": sil,
                "calinski_harabasz": ch,
                "davies_bouldin": dbi,
                "cluster_balance_score": balance,
                "notes": "formal_weighted_fusion_run",
            }
        )

        for cid, csize in s.value_counts().sort_index().items():
            cluster_size_rows.append(
                {
                    "method": "proposed_weighted_fusion",
                    "modality_set": modality_set_name,
                    "k": int(k),
                    "cluster_id": int(cid),
                    "cluster_size": int(csize),
                }
            )

    metrics_df = pd.DataFrame(metric_rows)
    best_k = int(metrics_df.sort_values("silhouette", ascending=False).iloc[0]["k"])
    emb = PCA(n_components=2, random_state=random_state).fit_transform(x)
    plt.figure(figsize=(7, 6))
    plt.scatter(emb[:, 0], emb[:, 1], c=labels_map[best_k].values, cmap="tab10", s=20)
    plt.title(f"Proposed Weighted Fusion Embedding (best k={best_k})")
    plt.xlabel("dim1")
    plt.ylabel("dim2")
    plt.tight_layout()
    plt.savefig(out_figures_dir / "proposed_embedding_umap_bestk.png", dpi=220)
    plt.close()

    weight_plot = weight_table[["modality", "final_weight"]].copy().sort_values("final_weight", ascending=False)
    plt.figure(figsize=(7, 4))
    sns.barplot(data=weight_plot, x="modality", y="final_weight", color="#457b9d")
    plt.title("Proposed Modality Weights")
    plt.tight_layout()
    plt.savefig(out_figures_dir / "proposed_weight_barplot.png", dpi=220)
    plt.close()

    return {
        "fused": fused,
        "metrics": metrics_df,
        "cluster_sizes": pd.DataFrame(cluster_size_rows),
        "labels": labels_map,
        "best_k": best_k,
        "fused_summary": pd.DataFrame(
            [
                {
                    "method": "proposed_weighted_fusion",
                    "modality_set": modality_set_name,
                    "n_samples": int(fused.shape[0]),
                    "matrix_mean": float(fused.values.mean()),
                    "matrix_std": float(fused.values.std()),
                    "matrix_min": float(fused.values.min()),
                    "matrix_max": float(fused.values.max()),
                }
            ]
        ),
    }


def dry_run_weight_estimation(score_table: pd.DataFrame, params: dict[str, Any]) -> pd.DataFrame:
    weights = compute_modality_weights(
        score_table,
        alpha=float(params.get("alpha", 1.0)),
        beta=float(params.get("beta", 1.0)),
        delta=float(params.get("delta", 1.0)),
        gamma=float(params.get("gamma", 1.0)),
    )
    return pd.DataFrame({"modality": weights.index, "weight": weights.values, "stage": "round0_dry_run"})
