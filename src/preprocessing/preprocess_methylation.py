"""Methylation preprocessing utilities.

Input:
- data/raw/methylation/HumanMethylation450 (probe x sample)
- optional probe map file

Output:
- round1 preview outputs
- round2 modeling matrix and processing logs

Purpose:
- execute chunked scan for large matrix
- compute missingness and dimension summary
- provide round2 modeling-ready matrix (sample x feature)
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from src.preprocessing.preprocess_logging import append_preprocessing_dimension_change
from src.utils.io_utils import ensure_dir


def _normalize_sample_id(sample_id: str) -> str:
    sid = str(sample_id).strip().upper()
    return sid[:16] if len(sid) >= 16 else sid


def run_methylation_round1(
    input_path: Path,
    summary_output_path: Path,
    preview_output_path: Path,
    selected_samples: Iterable[str] | None = None,
    probe_map_path: Path | None = None,
    top_var_preview: int = 2000,
    preview_scan_rows: int = 50000,
    chunk_size: int = 5000,
    log_path: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    ensure_dir(summary_output_path.parent)

    header = pd.read_csv(input_path, sep="\t", nrows=0)
    sample_cols = header.columns[1:]
    normalized_cols = [_normalize_sample_id(c) for c in sample_cols]

    if selected_samples is not None:
        selected_set = set(selected_samples)
        keep_idx = [i for i, sid in enumerate(normalized_cols) if sid in selected_set]
    else:
        keep_idx = list(range(len(normalized_cols)))

    usecols = [0] + [i + 1 for i in keep_idx]
    kept_sample_ids = [normalized_cols[i] for i in keep_idx]

    total_rows = 0
    total_missing = 0
    total_values = 0
    preview_parts: list[pd.DataFrame] = []

    for chunk in pd.read_csv(input_path, sep="\t", usecols=usecols, chunksize=chunk_size, dtype=str):
        chunk.columns = ["probe_id"] + kept_sample_ids
        value_df = chunk.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")

        total_rows += chunk.shape[0]
        total_missing += int(value_df.isna().sum().sum())
        total_values += int(value_df.size)

        if sum(x.shape[0] for x in preview_parts) < preview_scan_rows:
            tmp = pd.concat([chunk.iloc[:, :1], value_df], axis=1)
            preview_parts.append(tmp)

    miss_rate = (total_missing / total_values) if total_values > 0 else 0.0

    preview_df = pd.concat(preview_parts, axis=0, ignore_index=True)
    preview_df = preview_df.iloc[:preview_scan_rows, :]
    values = preview_df.iloc[:, 1:]
    variances = values.var(axis=1, skipna=True)
    top_idx = variances.sort_values(ascending=False).head(min(top_var_preview, preview_df.shape[0])).index
    preview_top = preview_df.loc[top_idx].reset_index(drop=True)

    # Optional probe map join hook: record available mapping coverage for preview probes.
    probe_map_join_note = "probe_map_not_used"
    if probe_map_path is not None and probe_map_path.exists():
        probe_map_df = pd.read_csv(probe_map_path, sep="\t", dtype=str)
        map_key = probe_map_df.columns[0]
        mapped = preview_top[["probe_id"]].merge(probe_map_df[[map_key]], left_on="probe_id", right_on=map_key, how="left")
        matched = int(mapped[map_key].notna().sum())
        probe_map_join_note = f"probe_map_matched_in_preview={matched}/{preview_top.shape[0]}"

    summary_df = pd.DataFrame(
        [
            {
                "modality": "methylation",
                "read_mode": "chunked_read",
                "chunk_size": chunk_size,
                "sample_count_kept": len(kept_sample_ids),
                "probe_count_scanned": total_rows,
                "missing_cells": total_missing,
                "total_cells": total_values,
                "missing_rate": miss_rate,
                "preview_rows_for_variance": preview_df.shape[0],
                "top_var_preview_features": preview_top.shape[0],
                "notes": probe_map_join_note,
            }
        ]
    )

    summary_df.to_csv(summary_output_path, index=False, encoding="utf-8")
    preview_top.to_csv(preview_output_path, index=False, encoding="utf-8")

    if log_path is not None:
        append_preprocessing_dimension_change(
            log_path=log_path,
            modality="methylation",
            input_shape=(len(normalized_cols), total_rows),
            output_shape=(len(kept_sample_ids), preview_top.shape[0]),
            filtering_steps=f"chunked_scan;missing_stats;top_var_preview<={top_var_preview}",
            read_mode="chunked_read",
            notes=f"full_matrix_scanned_for_stats;{probe_map_join_note}",
        )

    return summary_df, preview_top


def _load_probe_map(probe_map_path: Path | None) -> pd.DataFrame:
    if probe_map_path is None or not probe_map_path.exists():
        return pd.DataFrame(columns=["probe_id", "gene", "chrom", "chromStart", "chromEnd", "strand"])

    probe_map_df = pd.read_csv(probe_map_path, sep="\t", dtype=str)
    probe_map_df = probe_map_df.rename(columns={probe_map_df.columns[0]: "probe_id"})
    return probe_map_df


def run_methylation_round2_modeling(
    input_path: Path,
    output_matrix_path: Path,
    summary_output_path: Path,
    filter_log_output_path: Path,
    selected_samples: Iterable[str] | None = None,
    probe_map_path: Path | None = None,
    missingness_filter_threshold: float = 0.30,
    variance_filter_quantile: float = 0.75,
    max_features: int = 5000,
    chunk_size: int = 2000,
    figures_dir: Path | None = None,
    log_path: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build round2 methylation modeling matrix (sample x selected probes).

    Current round uses probe-level selection. Probe-to-gene/promoter aggregation is left
    as an explicit next step and tracked in notes.
    """
    ensure_dir(output_matrix_path.parent)
    ensure_dir(summary_output_path.parent)
    ensure_dir(filter_log_output_path.parent)
    if figures_dir is not None:
        ensure_dir(figures_dir)

    header = pd.read_csv(input_path, sep="\t", nrows=0)
    sample_cols = header.columns[1:]
    normalized_cols = [_normalize_sample_id(c) for c in sample_cols]

    if selected_samples is not None:
        selected_set = set(selected_samples)
        keep_idx = [i for i, sid in enumerate(normalized_cols) if sid in selected_set]
    else:
        keep_idx = list(range(len(normalized_cols)))

    kept_sample_ids = [normalized_cols[i] for i in keep_idx]
    usecols = [0] + [i + 1 for i in keep_idx]

    # Pass 1: compute per-probe missingness and variance statistics.
    stat_rows: list[pd.DataFrame] = []
    total_probe_count = 0
    for chunk in pd.read_csv(input_path, sep="\t", usecols=usecols, chunksize=chunk_size, dtype=str):
        chunk.columns = ["probe_id"] + kept_sample_ids
        value_df = chunk.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")

        probe_missing = value_df.isna().mean(axis=1)
        probe_var = value_df.var(axis=1, skipna=True).fillna(0.0)
        stats_df = pd.DataFrame(
            {
                "probe_id": chunk["probe_id"].astype(str).values,
                "probe_missing_rate": probe_missing.values,
                "probe_variance": probe_var.values,
            }
        )
        stat_rows.append(stats_df)
        total_probe_count += chunk.shape[0]

    probe_stats = pd.concat(stat_rows, axis=0, ignore_index=True)

    miss_pass = probe_stats[probe_stats["probe_missing_rate"] <= missingness_filter_threshold].copy()
    if miss_pass.empty:
        raise ValueError("No methylation probes left after missingness filtering.")

    var_cut = float(miss_pass["probe_variance"].quantile(variance_filter_quantile))
    var_pass = miss_pass[miss_pass["probe_variance"] >= var_cut].copy()
    var_pass = var_pass.sort_values("probe_variance", ascending=False)
    if var_pass.shape[0] > max_features:
        var_pass = var_pass.head(max_features)

    selected_probe_ids = set(var_pass["probe_id"].astype(str).tolist())

    # Pass 2: materialize selected probes into matrix.
    kept_chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(input_path, sep="\t", usecols=usecols, chunksize=chunk_size, dtype=str):
        chunk.columns = ["probe_id"] + kept_sample_ids
        sub = chunk[chunk["probe_id"].astype(str).isin(selected_probe_ids)].copy()
        if sub.empty:
            continue
        values = sub.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
        # Probe-wise median imputation keeps beta-value scale for missing entries.
        values = values.T
        med = values.median(axis=0)
        values = values.fillna(med).fillna(0.0)
        values = values.T
        sub.iloc[:, 1:] = values.values
        kept_chunks.append(sub)

    if not kept_chunks:
        raise ValueError("No selected methylation probes were materialized in pass2.")

    selected_df = pd.concat(kept_chunks, axis=0, ignore_index=True)
    selected_df = selected_df.drop_duplicates(subset=["probe_id"], keep="first")
    selected_df = selected_df.set_index("probe_id")

    matrix = selected_df.T
    matrix.index = [_normalize_sample_id(x) for x in matrix.index]
    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)

    mean = matrix.mean(axis=0)
    std = matrix.std(axis=0).replace(0, np.nan)
    matrix = ((matrix - mean) / std).fillna(0.0)

    matrix.to_csv(output_matrix_path, encoding="utf-8")

    probe_map_df = _load_probe_map(probe_map_path)
    mapped_gene_fraction = np.nan
    if not probe_map_df.empty:
        merged = var_pass.merge(probe_map_df[["probe_id", "gene"]], on="probe_id", how="left")
        mapped_gene_fraction = float(merged["gene"].notna().mean()) if merged.shape[0] > 0 else np.nan

    summary_df = pd.DataFrame(
        [
            {
                "modality": "methylation",
                "input_probe_count": int(total_probe_count),
                "output_feature_count": int(matrix.shape[1]),
                "input_sample_count": int(len(kept_sample_ids)),
                "output_sample_count": int(matrix.shape[0]),
                "missingness_filter_threshold": float(missingness_filter_threshold),
                "variance_filter_strategy": f"quantile>={variance_filter_quantile}, top<={max_features}",
                "final_matrix_ready": bool(matrix.shape[0] > 0 and matrix.shape[1] > 0),
                "probe_map_gene_coverage": mapped_gene_fraction,
                "notes": "probe_level_modeling_matrix;probe_to_gene_aggregation_interface_reserved",
            }
        ]
    )
    summary_df.to_csv(summary_output_path, index=False, encoding="utf-8")

    filter_log_df = var_pass.copy()
    filter_log_df["missingness_filter_threshold"] = float(missingness_filter_threshold)
    filter_log_df["variance_filter_quantile"] = float(variance_filter_quantile)
    filter_log_df["selected_for_modeling"] = True
    filter_log_df.to_csv(filter_log_output_path, index=False, encoding="utf-8")

    if figures_dir is not None:
        sns.set_theme(style="whitegrid")

        plt.figure(figsize=(8, 5))
        sns.histplot(probe_stats["probe_missing_rate"], bins=60, color="#2a6f97")
        plt.axvline(missingness_filter_threshold, color="#c1121f", linestyle="--", linewidth=2)
        plt.title("Methylation Probe Missingness Distribution")
        plt.xlabel("Probe missing rate")
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(figures_dir / "methylation_missingness_summary.png", dpi=220)
        plt.close()

        plt.figure(figsize=(8, 5))
        sns.histplot(np.log10(miss_pass["probe_variance"].clip(lower=1e-12)), bins=60, color="#588157")
        plt.axvline(np.log10(max(var_cut, 1e-12)), color="#bc6c25", linestyle="--", linewidth=2)
        plt.title("Methylation Probe Variance Distribution (log10)")
        plt.xlabel("log10(variance)")
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(figures_dir / "methylation_variance_distribution.png", dpi=220)
        plt.close()

    if log_path is not None:
        append_preprocessing_dimension_change(
            log_path=log_path,
            modality="methylation_round2",
            input_shape=(len(kept_sample_ids), total_probe_count),
            output_shape=matrix.shape,
            filtering_steps=f"missing_rate<={missingness_filter_threshold};variance_quantile>={variance_filter_quantile};top<={max_features};zscore",
            read_mode="chunked_two_pass",
            notes="round2_modeling_matrix",
        )

    return matrix, summary_df, filter_log_df
