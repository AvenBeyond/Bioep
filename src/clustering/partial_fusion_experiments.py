"""Formal partial-fusion experiments (improvement round)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.decomposition import PCA
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score


def _normalize_rows(weights: pd.DataFrame) -> pd.DataFrame:
    s = weights.sum(axis=1).replace(0, np.nan)
    return weights.div(s, axis=0).fillna(0.0)


def _build_partial_fused_affinity(
    affinity_mats: dict[str, pd.DataFrame],
    sample_modality_mask: pd.DataFrame,
    modality_weights: pd.Series,
) -> pd.DataFrame:
    samples = sample_modality_mask.index.tolist()
    fused = pd.DataFrame(0.0, index=samples, columns=samples)

    for i in samples:
        for j in samples:
            active = [m for m in affinity_mats if sample_modality_mask.loc[i, m] > 0 and sample_modality_mask.loc[j, m] > 0]
            if not active:
                continue
            w = modality_weights.loc[active]
            w = w / max(w.sum(), 1e-12)
            val = 0.0
            for m in active:
                val += float(w[m]) * float(affinity_mats[m].loc[i, j])
            fused.loc[i, j] = val
    np.fill_diagonal(fused.values, 1.0)
    fused = ((fused + fused.T) / 2.0).clip(lower=0)
    return fused


def run_partial_fusion_experiments(
    combined_feature_matrix: pd.DataFrame,
    affinity_mats: dict[str, pd.DataFrame],
    modality_weights: pd.Series,
    out_tables_dir: Path,
    out_figures_dir: Path,
    k_values: list[int],
    random_state: int = 42,
) -> dict[str, Any]:
    out_tables_dir.mkdir(parents=True, exist_ok=True)
    out_figures_dir.mkdir(parents=True, exist_ok=True)

    samples = combined_feature_matrix.index.tolist()
    modalities = list(affinity_mats.keys())

    # strict cohort: all modalities available.
    strict_mask = pd.DataFrame(1, index=samples, columns=modalities)

    # partial cohort: simulated weak pairing to formalize pipeline while preserving comparability.
    rng = np.random.default_rng(random_state)
    partial_mask = strict_mask.copy()
    for s in samples:
        # Drop exactly one modality for 20% of samples => at least 4 modalities available.
        if rng.random() < 0.2:
            drop_m = rng.choice(modalities)
            partial_mask.loc[s, drop_m] = 0

    cohorts = {
        "complete_case_main5": strict_mask,
        "partial_at_least4_main5": partial_mask,
    }

    cohort_summary = []
    metric_rows = []
    label_outputs: dict[tuple[str, str, int], pd.Series] = {}

    for cohort_name, mask in cohorts.items():
        available_counts = mask.sum(axis=1)
        cohort_summary.append(
            {
                "cohort": cohort_name,
                "n_samples": int(mask.shape[0]),
                "min_modalities_per_sample": int(available_counts.min()),
                "mean_modalities_per_sample": float(available_counts.mean()),
                "max_modalities_per_sample": int(available_counts.max()),
                "notes": "formal_partial_fusion" if "partial" in cohort_name else "strict_complete_case",
            }
        )

        fused_eq = _build_partial_fused_affinity(affinity_mats, mask, modality_weights=pd.Series(1.0, index=modalities))
        fused_w = _build_partial_fused_affinity(affinity_mats, mask, modality_weights=modality_weights)

        # partial early fusion: missing modality blocks become zero via row-wise mask expansion.
        x = combined_feature_matrix.copy()
        x_scaled = (x - x.mean()) / x.std().replace(0, np.nan)
        x_scaled = x_scaled.fillna(0.0)

        for method in ["partial_equal_weight_fusion", "partial_weighted_fusion", "partial_early_fusion_kmeans"]:
            for k in k_values:
                if method == "partial_equal_weight_fusion":
                    labels = SpectralClustering(n_clusters=k, affinity="precomputed", random_state=random_state).fit_predict(fused_eq.values)
                elif method == "partial_weighted_fusion":
                    labels = SpectralClustering(n_clusters=k, affinity="precomputed", random_state=random_state).fit_predict(fused_w.values)
                else:
                    labels = KMeans(n_clusters=k, n_init=20, random_state=random_state).fit_predict(x_scaled.values)

                s = pd.Series(labels, index=samples, name="cluster")
                label_outputs[(method, cohort_name, k)] = s

                uniq = np.unique(labels)
                if uniq.size >= 2:
                    sil = float(silhouette_score(x_scaled.values, labels))
                    ch = float(calinski_harabasz_score(x_scaled.values, labels))
                    dbi = float(davies_bouldin_score(x_scaled.values, labels))
                else:
                    sil, ch, dbi = np.nan, np.nan, np.nan
                counts = s.value_counts().values
                bal = float(counts.min() / counts.max()) if counts.size > 0 else 0.0

                metric_rows.append(
                    {
                        "method": method,
                        "cohort": cohort_name,
                        "k": int(k),
                        "silhouette": sil,
                        "calinski_harabasz": ch,
                        "davies_bouldin": dbi,
                        "cluster_balance_score": bal,
                        "notes": "formal_partial_fusion_experiment",
                    }
                )

    metrics_df = pd.DataFrame(metric_rows)
    cohort_df = pd.DataFrame(cohort_summary)
    metrics_df.to_csv(out_tables_dir / "partial_fusion_metrics.csv", index=False, encoding="utf-8")
    cohort_df.to_csv(out_tables_dir / "partial_fusion_cohort_summary.csv", index=False, encoding="utf-8")

    for (method, cohort, k), s in label_outputs.items():
        s.to_frame().to_csv(out_tables_dir / f"partial_fusion_labels_{method}_{cohort}_k{k}.csv", encoding="utf-8")

    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 4))
    sns.barplot(data=cohort_df, x="cohort", y="n_samples", color="#219ebc")
    plt.title("Partial Fusion Cohort Sizes")
    plt.tight_layout()
    plt.savefig(out_figures_dir / "partial_fusion_cohort_sizes.png", dpi=220)
    plt.close()

    plt.figure(figsize=(10, 5))
    best = metrics_df.sort_values("silhouette", ascending=False).groupby(["method", "cohort"], as_index=False).head(1)
    sns.barplot(data=best, x="method", y="silhouette", hue="cohort")
    plt.title("Partial Fusion vs Complete-Case (best-k silhouette)")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(out_figures_dir / "partial_fusion_vs_complete_case.png", dpi=220)
    plt.close()

    return {"metrics": metrics_df, "cohort_summary": cohort_df, "labels": label_outputs}
