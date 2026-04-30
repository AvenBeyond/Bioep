"""Preprocessing log helpers.

Input:
- per-modality preprocessing run metadata

Output:
- append rows to results/logs/preprocessing_dimension_changes.csv
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import csv

from src.utils.io_utils import ensure_dir


def append_preprocessing_dimension_change(
    log_path: Path,
    modality: str,
    input_shape: tuple[int, int],
    output_shape: tuple[int, int],
    filtering_steps: str,
    read_mode: str,
    notes: str,
) -> None:
    ensure_dir(log_path.parent)
    file_exists = log_path.exists()

    fieldnames = [
        "timestamp",
        "modality",
        "input_shape",
        "output_shape",
        "sample_count_before",
        "sample_count_after",
        "feature_count_before",
        "feature_count_after",
        "filtering_steps",
        "read_mode",
        "notes",
    ]

    with log_path.open("a", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()
        writer.writerow(
            {
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "modality": modality,
                "input_shape": f"{input_shape[0]}x{input_shape[1]}",
                "output_shape": f"{output_shape[0]}x{output_shape[1]}",
                "sample_count_before": input_shape[0],
                "sample_count_after": output_shape[0],
                "feature_count_before": input_shape[1],
                "feature_count_after": output_shape[1],
                "filtering_steps": filtering_steps,
                "read_mode": read_mode,
                "notes": notes,
            }
        )
