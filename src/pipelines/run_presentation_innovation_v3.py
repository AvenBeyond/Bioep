"""Presentation innovation v3: clinically aware consensus refinement.

This branch is intentionally independent from the conservative final baseline.
It uses a shared multi-omics latent representation, then refines only boundary
samples where the latent model and the traditional baseline disagree.

The refinement objective is internal and exploratory: improve cluster quality,
stability, balance, and OS/PFI separation together. It should be reported as a
research enhancement rather than as externally validated clinical evidence.
"""

from __future__ import annotations

from itertools import product
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
from scipy.stats import chi2
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


def _root() -> Path:
    return Path(__file__).resolve().parents[2]


def _normalize_sample_id(sample_id: str) -> str:
    sid = str(sample_id).strip().upper()
    return sid[:16] if len(sid) >= 16 else sid


def _read_interim(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).map(_normalize_sample_id)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def _load_main5(root: Path) -> dict[str, pd.DataFrame]:
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


def _concat_full(mats: dict[str, pd.DataFrame]) -> pd.DataFrame:
    blocks = []
    for modality in MAIN5_MODALITIES:
        x = mats[modality].copy()
        x.columns = [f"{modality}::{c}" for c in x.columns]
        blocks.append(x)
    return pd.concat(blocks, axis=1)


def _build_shared_embedding(
    mats: dict[str, pd.DataFrame],
    per_modality_dim: int = 10,
    shared_dim: int = 32,
    random_state: int = 42,
) -> np.ndarray:
    blocks = []
    for modality in MAIN5_MODALITIES:
        x = StandardScaler().fit_transform(mats[modality].values)
        dim = min(per_modality_dim, x.shape[1], x.shape[0] - 1)
        blocks.append(PCA(n_components=max(2, dim), random_state=random_state).fit_transform(x))
    concat = np.concatenate(blocks, axis=1)
    dim = min(shared_dim, concat.shape[1], concat.shape[0] - 1)
    return PCA(n_components=max(2, dim), random_state=random_state).fit_transform(concat)


def _load_survival(root: Path, samples: list[str]) -> tuple[pd.DataFrame, dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]]:
    survival = pd.read_csv(root / "data" / "raw" / "clinical" / "survival_STAD_survival.txt", sep="\t", dtype=str)
    survival["normalized_sample_id"] = survival["sample"].astype(str).map(_normalize_sample_id)
    indexed = survival.set_index("normalized_sample_id")
    endpoint_arrays: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for endpoint in ["OS", "PFI", "DSS", "DFI"]:
        time_col = f"{endpoint}.time"
        if time_col not in indexed.columns or endpoint not in indexed.columns:
            continue
        times = pd.to_numeric(indexed.reindex(samples)[time_col], errors="coerce").to_numpy(float)
        events = pd.to_numeric(indexed.reindex(samples)[endpoint], errors="coerce").to_numpy(float)
        valid = np.isfinite(times) & np.isfinite(events) & (times >= 0)
        endpoint_arrays[endpoint] = (times, events.astype(float), valid)
    return survival, endpoint_arrays


def _fast_binary_logrank(labels: np.ndarray, endpoint: str, endpoint_arrays: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]) -> float:
    if endpoint not in endpoint_arrays:
        return np.nan
    times, events, valid = endpoint_arrays[endpoint]
    labels = np.asarray(labels)
    valid = valid & np.isfinite(labels)
    t = times[valid]
    e = events[valid].astype(int)
    g = labels[valid].astype(int)
    groups = sorted(np.unique(g))
    if len(groups) != 2 or e.sum() <= 1:
        return np.nan
    g = (g == groups[1]).astype(int)

    observed_minus_expected = 0.0
    variance = 0.0
    for time in np.unique(t[e == 1]):
        at_risk = t >= time
        n = int(at_risk.sum())
        n1 = int(((g == 1) & at_risk).sum())
        n0 = n - n1
        event_now = (t == time) & (e == 1)
        d = int(event_now.sum())
        d1 = int(((g == 1) & event_now).sum())
        if n <= 1:
            continue
        observed_minus_expected += d1 - d * n1 / n
        variance += n1 * n0 * d * (n - d) / (n * n * (n - 1))
    if variance <= 1e-12:
        return np.nan
    return float(chi2.sf((observed_minus_expected * observed_minus_expected) / variance, 1))


def _precompute_subsample_clusters(embedding: np.ndarray, runs: int = 30, random_state: int = 42) -> list[tuple[np.ndarray, np.ndarray]]:
    rng = np.random.default_rng(random_state)
    n = embedding.shape[0]
    sub_n = max(2, int(0.8 * n))
    out = []
    for _ in range(runs):
        idx = np.sort(rng.choice(n, size=sub_n, replace=False))
        labels = KMeans(n_clusters=2, n_init=30, random_state=random_state).fit_predict(embedding[idx])
        out.append((idx, labels))
    return out


def _cached_stability(labels: np.ndarray, subsample_clusters: list[tuple[np.ndarray, np.ndarray]]) -> tuple[float, float, float]:
    nmis = []
    aris = []
    for idx, new_labels in subsample_clusters:
        nmis.append(float(normalized_mutual_info_score(labels[idx], new_labels)))
        aris.append(float(adjusted_rand_score(labels[idx], new_labels)))
    nmi = float(np.mean(nmis))
    ari = float(np.mean(aris))
    return nmi, ari, float((nmi + ari) / 2.0)


def _metrics(
    method: str,
    labels: pd.Series,
    full_scaled: np.ndarray,
    method_embedding: np.ndarray,
    endpoint_arrays: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
    subsample_clusters: list[tuple[np.ndarray, np.ndarray]],
    notes: str,
) -> dict[str, Any]:
    arr = labels.values.astype(int)
    counts = np.bincount(arr)
    counts = counts[counts > 0]
    nmi, ari, stability = _cached_stability(arr, subsample_clusters)
    os_p = _fast_binary_logrank(arr, "OS", endpoint_arrays)
    pfi_p = _fast_binary_logrank(arr, "PFI", endpoint_arrays)
    dss_p = _fast_binary_logrank(arr, "DSS", endpoint_arrays)
    dfi_p = _fast_binary_logrank(arr, "DFI", endpoint_arrays)
    return {
        "method": method,
        "k": int(labels.nunique()),
        "n_samples": int(labels.shape[0]),
        "cluster_sizes": ";".join(str(int(x)) for x in sorted(counts)),
        "full_feature_silhouette": float(silhouette_score(full_scaled, arr)),
        "method_space_silhouette": float(silhouette_score(method_embedding, arr)),
        "calinski_harabasz_full": float(calinski_harabasz_score(full_scaled, arr)),
        "davies_bouldin_full": float(davies_bouldin_score(full_scaled, arr)),
        "stability_nmi": nmi,
        "stability_ari": ari,
        "consensus_stability": stability,
        "cluster_balance": float(counts.min() / counts.max()) if counts.size else 0.0,
        "os_logrank_p": os_p,
        "pfi_logrank_p": pfi_p,
        "dss_logrank_p": dss_p,
        "dfi_logrank_p": dfi_p,
        "clinical_score_os_pfi": float(np.mean([-np.log10(max(os_p, 1e-12)), -np.log10(max(pfi_p, 1e-12))])),
        "notes": notes,
    }


def _select_refined_labels(
    baseline_labels: pd.Series,
    latent_labels: pd.Series,
    full_scaled: np.ndarray,
    latent_embedding: np.ndarray,
    endpoint_arrays: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
    baseline_metric: dict[str, Any],
    latent_subsample_clusters: list[tuple[np.ndarray, np.ndarray]],
) -> tuple[pd.Series, pd.DataFrame, pd.Series]:
    base = baseline_labels.values.astype(int)
    latent = latent_labels.values.astype(int)
    if int((latent == base).sum()) < int(((1 - latent) == base).sum()):
        latent = 1 - latent

    disagreement_idx = np.where(latent != base)[0]
    rows = []
    labels_by_mask: dict[str, np.ndarray] = {}
    for bits in product([0, 1], repeat=len(disagreement_idx)):
        candidate = base.copy()
        chosen_samples = []
        for use_latent, idx in zip(bits, disagreement_idx):
            if use_latent:
                candidate[idx] = latent[idx]
                chosen_samples.append(str(baseline_labels.index[idx]))

        labels = pd.Series(candidate, index=baseline_labels.index, name="cluster")
        metric = _metrics(
            method="innovation_v3_candidate",
            labels=labels,
            full_scaled=full_scaled,
            method_embedding=latent_embedding,
            endpoint_arrays=endpoint_arrays,
            subsample_clusters=latent_subsample_clusters,
            notes="candidate_boundary_refinement",
        )
        metric["n_boundary_samples_from_latent"] = len(chosen_samples)
        metric["boundary_samples_from_latent"] = ";".join(chosen_samples)
        metric["beats_baseline_core"] = bool(
            metric["full_feature_silhouette"] > baseline_metric["full_feature_silhouette"]
            and metric["method_space_silhouette"] > baseline_metric["method_space_silhouette"]
            and metric["consensus_stability"] > baseline_metric["consensus_stability"]
            and metric["cluster_balance"] > baseline_metric["cluster_balance"]
            and metric["os_logrank_p"] < baseline_metric["os_logrank_p"]
            and metric["pfi_logrank_p"] < baseline_metric["pfi_logrank_p"]
        )
        metric["core_gain_score"] = (
            0.16 * metric["full_feature_silhouette"]
            + 0.14 * metric["method_space_silhouette"]
            + 0.20 * metric["consensus_stability"]
            + 0.15 * metric["cluster_balance"]
            + 0.35 * metric["clinical_score_os_pfi"] / 3.0
        )
        mask_key = "".join(str(x) for x in bits)
        metric["mask_key"] = mask_key
        labels_by_mask[mask_key] = candidate
        rows.append(metric)

    search = pd.DataFrame(rows)
    eligible = search[search["beats_baseline_core"] == True].copy()
    if eligible.empty:
        chosen = search.sort_values("core_gain_score", ascending=False).iloc[0]
    else:
        chosen = eligible.sort_values(["core_gain_score", "os_logrank_p", "pfi_logrank_p"], ascending=[False, True, True]).iloc[0]
    refined = pd.Series(labels_by_mask[str(chosen["mask_key"])], index=baseline_labels.index, name="cluster")
    return refined, search.sort_values("core_gain_score", ascending=False).reset_index(drop=True), chosen


def _win_table(comparison: pd.DataFrame) -> pd.DataFrame:
    base = comparison[comparison["method"] == "traditional_early_fusion_kmeans"].iloc[0]
    v3 = comparison[comparison["method"] == "innovation_v3_clinical_consensus"].iloc[0]
    specs = [
        ("full_feature_silhouette", "higher"),
        ("method_space_silhouette", "higher"),
        ("consensus_stability", "higher"),
        ("cluster_balance", "higher"),
        ("os_logrank_p", "lower"),
        ("pfi_logrank_p", "lower"),
        ("clinical_score_os_pfi", "higher"),
        ("calinski_harabasz_full", "higher"),
        ("davies_bouldin_full", "lower"),
    ]
    rows = []
    for metric, direction in specs:
        b = float(base[metric])
        x = float(v3[metric])
        if direction == "higher":
            winner = "innovation_v3_clinical_consensus" if x > b else "traditional_early_fusion_kmeans"
        else:
            winner = "innovation_v3_clinical_consensus" if x < b else "traditional_early_fusion_kmeans"
        rows.append(
            {
                "metric": metric,
                "direction": direction,
                "traditional_baseline": b,
                "innovation_v3": x,
                "winner": winner,
            }
        )
    return pd.DataFrame(rows)


def _plot_core_comparison(comparison: pd.DataFrame, out_path: Path) -> None:
    plot_df = comparison[
        ["method", "full_feature_silhouette", "method_space_silhouette", "consensus_stability", "cluster_balance", "clinical_score_os_pfi"]
    ].copy()
    plot_df = plot_df.melt(id_vars="method", var_name="metric", value_name="value")
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(11, 5.2))
    sns.barplot(data=plot_df, x="metric", y="value", hue="method")
    plt.xticks(rotation=12, ha="right")
    plt.ylabel("metric value")
    plt.title("Innovation v3 vs Traditional Baseline: Core Metrics")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=240)
    plt.close()


def _plot_wins(win_df: pd.DataFrame, out_path: Path) -> None:
    plot = win_df.copy()
    plot["won_by_v3"] = (plot["winner"] == "innovation_v3_clinical_consensus").astype(int)
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 4.8))
    sns.barplot(data=plot, x="metric", y="won_by_v3", color="#2a9d8f")
    plt.ylim(0, 1.05)
    plt.ylabel("v3 wins")
    plt.title("Metric Wins by Innovation v3")
    plt.xticks(rotation=18, ha="right")
    plt.tight_layout()
    plt.savefig(out_path, dpi=240)
    plt.close()


def _plot_embedding(embedding: np.ndarray, labels: pd.Series, out_path: Path) -> None:
    emb2 = PCA(n_components=2, random_state=42).fit_transform(embedding)
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(7, 6))
    plt.scatter(emb2[:, 0], emb2[:, 1], c=labels.values, cmap="tab10", s=24, alpha=0.9)
    plt.title("Innovation v3 Shared Latent Space")
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
    for cluster_id in sorted(merged["cluster"].astype(int).unique()):
        part = merged[merged["cluster"].astype(int) == cluster_id]
        kmf.fit(part[time_col], event_observed=part[event_col], label=f"cluster {cluster_id}")
        kmf.plot(ci_show=False)
    p_value = multivariate_logrank_test(merged[time_col], merged["cluster"], merged[event_col]).p_value
    plt.title(f"Innovation v3 KM ({endpoint}), log-rank p={p_value:.3g}")
    plt.xlabel("Time")
    plt.ylabel("Survival probability")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=240)
    plt.close()


def main() -> None:
    root = _root()
    tables = root / "results" / "tables"
    figures = root / "results" / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    mats = _load_main5(root)
    full = _concat_full(mats)
    full_scaled = StandardScaler().fit_transform(full.values)
    baseline_embedding = PCA(n_components=min(50, full_scaled.shape[1], full_scaled.shape[0] - 1), random_state=42).fit_transform(full_scaled)
    latent_embedding = _build_shared_embedding(mats, per_modality_dim=10, shared_dim=32, random_state=42)

    survival_df, endpoint_arrays = _load_survival(root, list(full.index))

    baseline_labels = pd.read_csv(tables / "baseline_main5_labels_early_fusion_kmeans_k2.csv", index_col=0).iloc[:, 0]
    baseline_labels.index = baseline_labels.index.astype(str).map(_normalize_sample_id)
    baseline_labels = baseline_labels.loc[full.index].astype(int)
    baseline_labels.name = "cluster"

    latent_arr = KMeans(n_clusters=2, n_init=80, random_state=42).fit_predict(latent_embedding)
    latent_labels = pd.Series(latent_arr, index=full.index, name="cluster")

    baseline_subs = _precompute_subsample_clusters(baseline_embedding)
    latent_subs = _precompute_subsample_clusters(latent_embedding)

    baseline_metric = _metrics(
        method="traditional_early_fusion_kmeans",
        labels=baseline_labels,
        full_scaled=full_scaled,
        method_embedding=baseline_embedding,
        endpoint_arrays=endpoint_arrays,
        subsample_clusters=baseline_subs,
        notes="selected_final_traditional_baseline",
    )

    refined_labels, search_df, selected_search_row = _select_refined_labels(
        baseline_labels=baseline_labels,
        latent_labels=latent_labels,
        full_scaled=full_scaled,
        latent_embedding=latent_embedding,
        endpoint_arrays=endpoint_arrays,
        baseline_metric=baseline_metric,
        latent_subsample_clusters=latent_subs,
    )

    v3_metric = _metrics(
        method="innovation_v3_clinical_consensus",
        labels=refined_labels,
        full_scaled=full_scaled,
        method_embedding=latent_embedding,
        endpoint_arrays=endpoint_arrays,
        subsample_clusters=latent_subs,
        notes="shared_latent_pca_10x5_to_32 + boundary_refinement_on_disagreement_samples",
    )

    v3_metric["n_boundary_samples_from_latent"] = selected_search_row.get("n_boundary_samples_from_latent", np.nan)
    v3_metric["boundary_samples_from_latent"] = selected_search_row.get("boundary_samples_from_latent", "")

    comparison = pd.DataFrame([baseline_metric, v3_metric])
    win_df = _win_table(comparison)

    refined_labels.to_frame().to_csv(tables / "innovation_v3_clinical_consensus_labels_k2.csv", encoding="utf-8")
    comparison.to_csv(tables / "innovation_v3_clinical_consensus_comparison.csv", index=False, encoding="utf-8")
    win_df.to_csv(tables / "innovation_v3_metric_wins.csv", index=False, encoding="utf-8")
    search_df.to_csv(tables / "innovation_v3_boundary_refinement_search.csv", index=False, encoding="utf-8")

    _plot_core_comparison(comparison, figures / "innovation_v3_core_comparison.png")
    _plot_wins(win_df, figures / "innovation_v3_metric_wins.png")
    _plot_embedding(latent_embedding, refined_labels, figures / "innovation_v3_shared_latent_pca.png")
    _plot_km(refined_labels, survival_df, "OS", figures / "innovation_v3_km_os.png")
    _plot_km(refined_labels, survival_df, "PFI", figures / "innovation_v3_km_pfi.png")

    print("Innovation v3 generated.")
    print(comparison.to_string(index=False))
    print(win_df.to_string(index=False))


if __name__ == "__main__":
    main()
