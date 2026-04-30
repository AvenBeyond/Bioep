"""CNV preprocessing (round 1 executable).

Input:
- data/raw/cnv/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes (gene x sample)

Output:
- data/interim/cnv_round1.csv (sample x feature)

Purpose:
- normalize sample IDs
- enforce GISTIC-like coding and low-variance filtering
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from src.preprocessing.preprocess_logging import append_preprocessing_dimension_change
from src.utils.io_utils import ensure_dir


def _normalize_sample_id(sample_id: str) -> str:
    sid = str(sample_id).strip().upper()
    return sid[:16] if len(sid) >= 16 else sid


def preprocess_cnv_dataframe(
    df_gene_by_sample: pd.DataFrame,
    min_variance: float = 0.01,
    selected_samples: Iterable[str] | None = None,
) -> pd.DataFrame:
    matrix = df_gene_by_sample.copy().T
    matrix.index = [_normalize_sample_id(x) for x in matrix.index]

    matrix = matrix.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    matrix = matrix.clip(lower=-2, upper=2)
    matrix = matrix.round(0).astype(np.int8)

    if selected_samples is not None:
        selected_set = set(selected_samples)
        matrix = matrix.loc[matrix.index.isin(selected_set)]

    variances = matrix.var(axis=0)
    kept_features = variances[variances >= min_variance].index
    matrix = matrix.loc[:, kept_features]
    return matrix


def run_cnv_round1(
    input_path: Path,
    output_path: Path,
    min_variance: float = 0.01,
    selected_samples: Iterable[str] | None = None,
    log_path: Path | None = None,
) -> pd.DataFrame:
    raw_df = pd.read_csv(input_path, sep="\t", index_col=0, dtype=str)
    input_shape = (raw_df.shape[1], raw_df.shape[0])

    out_df = preprocess_cnv_dataframe(raw_df, min_variance=min_variance, selected_samples=selected_samples)
    ensure_dir(output_path.parent)
    out_df.to_csv(output_path, encoding="utf-8")

    if log_path is not None:
        append_preprocessing_dimension_change(
            log_path=log_path,
            modality="cnv",
            input_shape=input_shape,
            output_shape=out_df.shape,
            filtering_steps=f"clip_to_-2_2;round_to_int;variance>={min_variance}",
            read_mode="full_read",
            notes="round1_real_preprocessing",
        )

    return out_df
