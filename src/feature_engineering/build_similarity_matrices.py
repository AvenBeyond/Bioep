"""Build sample-sample similarity matrices for round-1 modeling.

Input:
- interim modality matrices (sample x feature)

Output:
- aligned matrices in data/interim/aligned_matrices
- similarity/affinity matrices in results/tables/similarity_matrices
- logs/tables and per-modality heatmaps
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances

from src.utils.config_utils import load_project_config
from src.utils.io_utils import ensure_dir


ROUND1_MODALITIES = ["mutation", "cnv", "rna", "mirna"]


def pearson_similarity(x: pd.DataFrame) -> pd.DataFrame:
    sim = np.corrcoef(x.values)
    return pd.DataFrame(sim, index=x.index, columns=x.index)


def cosine_affinity(x: pd.DataFrame) -> pd.DataFrame:
    sim = cosine_similarity(x.values)
    return pd.DataFrame(sim, index=x.index, columns=x.index)


def euclidean_affinity(x: pd.DataFrame, sigma: float = 1.0) -> pd.DataFrame:
    dist = euclidean_distances(x.values)
    sim = np.exp(-(dist ** 2) / (2 * sigma**2))
    return pd.DataFrame(sim, index=x.index, columns=x.index)


def _is_already_normalized(x: pd.DataFrame, atol: float = 0.2) -> bool:
    mean_ok = bool(np.all(np.abs(x.mean(axis=0).values) < atol))
    stds = x.std(axis=0).replace(0, np.nan).dropna().values
    std_ok = bool(np.all(np.abs(stds - 1.0) < atol)) if stds.size > 0 else False
    return mean_ok and std_ok


def _zscore_if_needed(x: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    matrix = x.copy()
    if _is_already_normalized(matrix):
        return matrix, "already_normalized"
    mean = matrix.mean(axis=0)
    std = matrix.std(axis=0).replace(0, np.nan)
    matrix = ((matrix - mean) / std).fillna(0.0)
    return matrix, "zscore_applied"


def _load_interim_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).str.upper().str[:16]
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def build_similarity_round1(project_root: Path, modalities: list[str] | None = None) -> dict[str, Any]:
    paths_cfg, _ = load_project_config(project_root)
    tables_dir = Path(paths_cfg["results"]["tables"])
    figures_dir = Path(paths_cfg["results"]["figures"])
    interim_dir = Path(paths_cfg["data"]["interim_root"])

    aligned_dir = interim_dir / "aligned_matrices"
    sim_dir = tables_dir / "similarity_matrices"
    ensure_dir(aligned_dir)
    ensure_dir(sim_dir)
    ensure_dir(figures_dir)

    use_modalities = modalities or ROUND1_MODALITIES
    interim_files = {m: interim_dir / f"{m}_round1.csv" for m in use_modalities}

    matrices: dict[str, pd.DataFrame] = {}
    for m, p in interim_files.items():
        if not p.exists():
            continue
        matrices[m] = _load_interim_matrix(p)

    if not matrices:
        raise FileNotFoundError("No round1 interim matrices found for similarity building.")

    cohort_samples = sorted(set.intersection(*[set(df.index) for df in matrices.values()]))
    if len(cohort_samples) == 0:
        raise ValueError("No intersected samples across selected modalities.")

    shape_rows: list[dict[str, Any]] = []
    build_log_rows: list[dict[str, Any]] = []
    similarities: dict[str, pd.DataFrame] = {}
    affinities: dict[str, pd.DataFrame] = {}

    sns.set_theme(style="white")

    for modality, mat in matrices.items():
        aligned = mat.loc[cohort_samples].copy()
        aligned, norm_status = _zscore_if_needed(aligned)
        aligned.to_csv(aligned_dir / f"{modality}_aligned.csv", encoding="utf-8")

        pear = pearson_similarity(aligned).clip(-1, 1)
        cos = cosine_affinity(aligned)
        aff = euclidean_affinity(aligned, sigma=1.0)

        pear.to_csv(sim_dir / f"{modality}_pearson_similarity.csv", encoding="utf-8")
        cos.to_csv(sim_dir / f"{modality}_cosine_similarity.csv", encoding="utf-8")
        aff.to_csv(sim_dir / f"{modality}_euclidean_affinity.csv", encoding="utf-8")

        similarities[modality] = pear
        affinities[modality] = aff

        shape_rows.append(
            {
                "modality": modality,
                "input_samples": mat.shape[0],
                "input_features": mat.shape[1],
                "aligned_samples": aligned.shape[0],
                "aligned_features": aligned.shape[1],
                "notes": norm_status,
            }
        )

        for method, out_mat in {
            "pearson": pear,
            "cosine": cos,
            "euclidean_affinity": aff,
        }.items():
            build_log_rows.append(
                {
                    "modality": modality,
                    "input_samples": mat.shape[0],
                    "input_features": mat.shape[1],
                    "aligned_samples": aligned.shape[0],
                    "aligned_features": aligned.shape[1],
                    "similarity_method": method,
                    "output_matrix_shape": f"{out_mat.shape[0]}x{out_mat.shape[1]}",
                    "notes": norm_status,
                }
            )

        plt.figure(figsize=(7, 6))
        sns.heatmap(pear.iloc[:100, :100], cmap="vlag", center=0, cbar=True)
        plt.title(f"Sample Similarity Heatmap ({modality}, first 100 samples)")
        plt.xlabel("Samples")
        plt.ylabel("Samples")
        plt.tight_layout()
        plt.savefig(figures_dir / f"sample_similarity_heatmap_{modality}.png", dpi=220)
        plt.close()

    shape_df = pd.DataFrame(shape_rows)
    log_df = pd.DataFrame(build_log_rows)
    shape_df.to_csv(tables_dir / "modality_modeling_shapes.csv", index=False, encoding="utf-8")
    log_df.to_csv(tables_dir / "similarity_build_log.csv", index=False, encoding="utf-8")

    return {
        "cohort_samples": cohort_samples,
        "modalities": sorted(matrices.keys()),
        "similarities": similarities,
        "affinities": affinities,
        "aligned_dir": aligned_dir,
        "similarity_dir": sim_dir,
    }


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[2]
    output = build_similarity_round1(root)
    print(f"Built similarity for modalities: {output['modalities']}")
    print(f"Cohort size: {len(output['cohort_samples'])}")
