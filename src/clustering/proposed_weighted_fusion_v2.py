"""Interaction-aware weighted fusion v2 (frontier-inspired, lightweight)."""

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


def _safe_minmax(x: pd.Series) -> pd.Series:
    lo = float(x.min())
    hi = float(x.max())
    if hi - lo < 1e-12:
        return pd.Series(np.zeros(len(x)), index=x.index)
    return (x - lo) / (hi - lo)


def _build_interaction_matrix(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    # Geometric mean fusion encourages agreement across modalities.
    a_clip = a.clip(lower=0)
    b_clip = b.clip(lower=0)
    inter = np.sqrt(a_clip.values * b_clip.values)
    out = pd.DataFrame(inter, index=a.index, columns=a.columns)
    return ((out + out.T) / 2.0).clip(lower=0)


def _interaction_coupling_score(a: pd.DataFrame, b: pd.DataFrame) -> float:
    # Correlation on upper triangle as a robust lightweight CKA-like proxy.
    iu = np.triu_indices(a.shape[0], k=1)
    va = a.values[iu]
    vb = b.values[iu]
    if np.std(va) < 1e-12 or np.std(vb) < 1e-12:
        return 0.0
    c = np.corrcoef(va, vb)[0, 1]
    if not np.isfinite(c):
        return 0.0
    return float(max(0.0, c))


def run_weighted_fusion_v2(
    affinity_mats: dict[str, pd.DataFrame],
    modality_weights: pd.DataFrame,
    combined_feature_matrix: pd.DataFrame,
    interaction_pairs: list[tuple[str, str]] | None,
    lambda_interaction: float,
    k_values: list[int],
    out_tables_dir: Path,
    out_figures_dir: Path,
    prefix: str = "proposed_v2",
    random_state: int = 42,
) -> dict[str, Any]:
    out_tables_dir.mkdir(parents=True, exist_ok=True)
    out_figures_dir.mkdir(parents=True, exist_ok=True)

    weights = modality_weights.set_index("modality")["final_weight"].to_dict()
    base_fused = None
    for m, mat in affinity_mats.items():
        if m not in weights:
            continue
        base_fused = mat * weights[m] if base_fused is None else base_fused + mat * weights[m]
    if base_fused is None:
        raise ValueError("No modality overlap for v2 fusion.")

    if interaction_pairs is None:
        interaction_pairs = [
            ("rna", "methylation"),
            ("rna", "mirna"),
            ("cnv", "rna"),
            ("mutation", "rna"),
        ]

    interaction_rows: list[dict[str, Any]] = []
    inter_fused = pd.DataFrame(0.0, index=base_fused.index, columns=base_fused.columns)
    valid_pairs = []
    scores = []
    inter_mats: dict[tuple[str, str], pd.DataFrame] = {}

    for m, n in interaction_pairs:
        if m not in affinity_mats or n not in affinity_mats:
            continue
        s_mn = _build_interaction_matrix(affinity_mats[m], affinity_mats[n])
        c = _interaction_coupling_score(affinity_mats[m], affinity_mats[n])
        valid_pairs.append((m, n))
        scores.append(c)
        inter_mats[(m, n)] = s_mn

    if scores:
        v = pd.Series(scores, index=[f"{m}__{n}" for m, n in valid_pairs])
        v_norm = _safe_minmax(v)
        if v_norm.sum() <= 1e-12:
            v_norm = pd.Series(np.ones(len(v_norm)) / len(v_norm), index=v_norm.index)
        else:
            v_norm = v_norm / v_norm.sum()

        for (m, n), vn in zip(valid_pairs, v_norm.values):
            inter_fused += inter_mats[(m, n)] * float(vn)
            interaction_rows.append(
                {
                    "interaction": f"{m}-{n}",
                    "I_mn": float(scores[len(interaction_rows)]),
                    "v_mn": float(vn),
                    "notes": "geometric_mean_affinity + upper_triangle_correlation",
                }
            )

    fused = ((base_fused + lambda_interaction * inter_fused) + (base_fused + lambda_interaction * inter_fused).T) / 2.0
    fused = fused.clip(lower=0)

    modality_weights.to_csv(out_tables_dir / f"{prefix}_weight_components.csv", index=False, encoding="utf-8")
    pd.DataFrame(interaction_rows).to_csv(out_tables_dir / f"{prefix}_interaction_components.csv", index=False, encoding="utf-8")

    sns.set_theme(style="white")
    plt.figure(figsize=(7, 6))
    sns.heatmap(fused.iloc[:100, :100], cmap="mako")
    plt.title("Proposed v2 Fused Network (first 100)")
    plt.tight_layout()
    plt.savefig(out_figures_dir / f"{prefix}_fused_network_heatmap.png", dpi=220)
    plt.close()

    plt.figure(figsize=(7, 4))
    tmp = modality_weights[["modality", "final_weight"]].sort_values("final_weight", ascending=False)
    sns.barplot(data=tmp, x="modality", y="final_weight", color="#3a86ff")
    plt.title("Proposed v2 Modality Weights")
    plt.tight_layout()
    plt.savefig(out_figures_dir / f"{prefix}_weight_barplot.png", dpi=220)
    plt.close()

    if interaction_rows:
        i_df = pd.DataFrame(interaction_rows).sort_values("v_mn", ascending=False)
        plt.figure(figsize=(7, 4))
        sns.barplot(data=i_df, x="interaction", y="v_mn", color="#ff006e")
        plt.title("Proposed v2 Interaction Weights")
        plt.tight_layout()
        plt.savefig(out_figures_dir / f"{prefix}_interaction_barplot.png", dpi=220)
        plt.close()

    x = combined_feature_matrix.loc[fused.index].values
    metric_rows: list[dict[str, Any]] = []
    size_rows: list[dict[str, Any]] = []
    labels_map: dict[int, pd.Series] = {}

    for k in k_values:
        labels = SpectralClustering(n_clusters=k, affinity="precomputed", random_state=random_state).fit_predict(fused.values)
        s = pd.Series(labels, index=fused.index, name="cluster")
        labels_map[k] = s
        s.to_frame().to_csv(out_tables_dir / f"{prefix}_labels_k{k}.csv", encoding="utf-8")

        uniq = np.unique(labels)
        if uniq.size >= 2:
            sil = float(silhouette_score(x, labels))
            ch = float(calinski_harabasz_score(x, labels))
            dbi = float(davies_bouldin_score(x, labels))
        else:
            sil, ch, dbi = np.nan, np.nan, np.nan
        counts = s.value_counts().values
        bal = float(counts.min() / counts.max()) if counts.size > 0 else 0.0
        metric_rows.append(
            {
                "method": "proposed_weighted_fusion_v2",
                "variant": "interaction_aware",
                "k": int(k),
                "n_samples": int(len(s)),
                "silhouette": sil,
                "calinski_harabasz": ch,
                "davies_bouldin": dbi,
                "cluster_balance_score": bal,
                "lambda_interaction": float(lambda_interaction),
                "notes": "improvement_v2",
            }
        )
        for cid, csize in s.value_counts().sort_index().items():
            size_rows.append({"k": int(k), "cluster_id": int(cid), "cluster_size": int(csize), "method": "proposed_weighted_fusion_v2"})

    metrics_df = pd.DataFrame(metric_rows)
    sizes_df = pd.DataFrame(size_rows)
    metrics_df.to_csv(out_tables_dir / f"{prefix}_metrics.csv", index=False, encoding="utf-8")
    sizes_df.to_csv(out_tables_dir / f"{prefix}_cluster_sizes.csv", index=False, encoding="utf-8")

    return {
        "fused": fused,
        "metrics": metrics_df,
        "cluster_sizes": sizes_df,
        "labels": labels_map,
        "interaction_components": pd.DataFrame(interaction_rows),
    }
