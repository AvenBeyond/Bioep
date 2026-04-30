"""RPPA preprocessing (round 1 executable).

Input:
- data/raw/rppa/RPPA (feature x sample)

Output:
- data/interim/rppa_round1.csv (sample x feature)

Purpose:
- treat RPPA as supplementary modality with basic filtering and scaling
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


def preprocess_rppa_dataframe(
    df_feature_by_sample: pd.DataFrame,
    min_nonzero_rate: float = 0.1,
    selected_samples: Iterable[str] | None = None,
) -> pd.DataFrame:
    matrix = df_feature_by_sample.copy().T
    matrix.index = [_normalize_sample_id(x) for x in matrix.index]
    matrix = matrix.apply(pd.to_numeric, errors="coerce")

    # RPPA may contain sparse NA blocks; fill with feature median.
    matrix = matrix.apply(lambda col: col.fillna(col.median()), axis=0)
    matrix = matrix.fillna(0.0)

    if selected_samples is not None:
        selected_set = set(selected_samples)
        matrix = matrix.loc[matrix.index.isin(selected_set)]

    nonzero_rate = (matrix != 0).mean(axis=0)
    matrix = matrix.loc[:, nonzero_rate >= min_nonzero_rate]

    mean = matrix.mean(axis=0)
    std = matrix.std(axis=0).replace(0, np.nan)
    matrix = ((matrix - mean) / std).fillna(0.0)
    return matrix


def run_rppa_round1(
    input_path: Path,
    output_path: Path,
    min_nonzero_rate: float = 0.1,
    selected_samples: Iterable[str] | None = None,
    log_path: Path | None = None,
) -> pd.DataFrame:
    raw_df = pd.read_csv(input_path, sep="\t", index_col=0, dtype=str)
    input_shape = (raw_df.shape[1], raw_df.shape[0])

    out_df = preprocess_rppa_dataframe(raw_df, min_nonzero_rate=min_nonzero_rate, selected_samples=selected_samples)
    ensure_dir(output_path.parent)
    out_df.to_csv(output_path, encoding="utf-8")

    if log_path is not None:
        append_preprocessing_dimension_change(
            log_path=log_path,
            modality="rppa",
            input_shape=input_shape,
            output_shape=out_df.shape,
            filtering_steps=f"median_impute;nonzero_rate>={min_nonzero_rate};zscore",
            read_mode="full_read",
            notes="round1_real_preprocessing_supplementary",
        )

    return out_df
