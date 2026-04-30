"""Ablation study utilities for final model selection.

Input: modular experiment settings
Output: comparison table for ablation variants
Purpose: test innovation contributions and robustness.
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


def compute_ablation_rank(ablation_df: pd.DataFrame) -> pd.DataFrame:
    df = ablation_df.copy()
    # Lower p-values are better, so convert to separation scores.
    df["os_sep"] = -np.log10(df["logrank_os_p"].clip(lower=1e-12))
    df["pfi_sep"] = -np.log10(df["logrank_pfi_p"].clip(lower=1e-12))

    score = (
        0.30 * df["silhouette"].fillna(0.0)
        + 0.25 * df["stability_score"].fillna(0.0)
        + 0.20 * df["cluster_balance_score"].fillna(0.0)
        + 0.15 * df["os_sep"].fillna(0.0)
        + 0.10 * df["pfi_sep"].fillna(0.0)
    )
    df["overall_score"] = score
    df = df.sort_values("overall_score", ascending=False).reset_index(drop=True)
    df["overall_rank"] = df.index + 1
    return df


def save_ablation_plots(ablation_df: pd.DataFrame, figures_dir: Path) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")

    plt.figure(figsize=(12, 5))
    top = ablation_df.sort_values("overall_rank").head(15)
    sns.barplot(data=top, x="experiment_name", y="silhouette", color="#264653")
    plt.xticks(rotation=60, ha="right")
    plt.title("Ablation Comparison (Top by Rank)")
    plt.tight_layout()
    plt.savefig(figures_dir / "ablation_comparison.png", dpi=220)
    plt.close()

    weight_df = ablation_df[ablation_df["experiment_name"].str.contains("weight_")]
    if not weight_df.empty:
        plt.figure(figsize=(10, 4))
        sns.barplot(data=weight_df, x="experiment_name", y="silhouette", color="#2a9d8f")
        plt.xticks(rotation=45, ha="right")
        plt.title("Weight-Component Ablation")
        plt.tight_layout()
        plt.savefig(figures_dir / "weight_component_ablation.png", dpi=220)
        plt.close()

    mod_df = ablation_df[ablation_df["experiment_name"].str.contains("drop_")]
    if not mod_df.empty:
        plt.figure(figsize=(10, 4))
        sns.barplot(data=mod_df, x="experiment_name", y="silhouette", color="#e76f51")
        plt.xticks(rotation=45, ha="right")
        plt.title("Modality Ablation")
        plt.tight_layout()
        plt.savefig(figures_dir / "modality_ablation.png", dpi=220)
        plt.close()


def write_ablation_outputs(ablation_df: pd.DataFrame, tables_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    tables_dir.mkdir(parents=True, exist_ok=True)
    ranked = compute_ablation_rank(ablation_df)
    ranked.to_csv(tables_dir / "ablation_results.csv", index=False, encoding="utf-8")
    rank_summary = ranked[["experiment_name", "modality_set", "k", "overall_rank", "overall_score", "notes"]].copy()
    rank_summary.to_csv(tables_dir / "ablation_rank_summary.csv", index=False, encoding="utf-8")
    return ranked, rank_summary
