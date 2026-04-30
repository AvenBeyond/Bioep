"""Load all modalities and build a sample inventory.

Input:
- paths.yaml resolved raw file paths

Output:
- lightweight per-file inventory: existence, shape estimate, head preview
- sample inventory table from real header parsing

Purpose:
- First-round data entry point before full cleaning/alignment.

TODO:
- Extend to row-oriented sample tables if future datasets differ.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from src.utils.io_utils import file_exists, read_head


OMICS_KEYS = [
    "mutation",
    "cnv",
    "rna",
    "methylation",
    "methylation_probe_map",
    "mirna",
    "rppa",
    "clinical_matrix",
    "survival",
]


def estimate_shape_tsv(path: Path) -> tuple[int, int]:
    line_count = 0
    col_count = 0
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            if i == 0:
                col_count = len(line.rstrip("\n").split("\t"))
            line_count += 1
    return line_count, col_count


def normalize_sample_id(sample_id: str) -> str:
    # Candidate strategy: uppercase + retain first 15 chars when possible.
    sid = str(sample_id).strip().upper()
    return sid[:15] if len(sid) >= 15 else sid


def build_file_inventory(paths_cfg: dict[str, Any]) -> pd.DataFrame:
    rows = []
    for key in OMICS_KEYS:
        path = Path(paths_cfg["data"][key])
        exists = file_exists(path)
        line_count, col_count = (estimate_shape_tsv(path) if exists else (0, 0))
        head = read_head(path, n_lines=2) if exists else []
        rows.append(
            {
                "key": key,
                "path": str(path),
                "exists": exists,
                "line_count": line_count,
                "column_count": col_count,
                "head_line_1": head[0] if len(head) > 0 else "",
                "head_line_2": head[1] if len(head) > 1 else "",
            }
        )
    return pd.DataFrame(rows)


def build_sample_inventory(paths_cfg: dict[str, Any]) -> pd.DataFrame:
    # Parse sample IDs from header columns for current TCGA matrix layout.
    records = []
    for key in OMICS_KEYS:
        path = Path(paths_cfg["data"][key])
        if not file_exists(path):
            continue
        with path.open("r", encoding="utf-8", errors="replace") as f:
            header = f.readline().rstrip("\n").split("\t")
        sample_ids = header[1:] if len(header) > 1 else []
        for sid in sample_ids:
            records.append(
                {
                    "modality": key,
                    "raw_sample_id": sid,
                    "normalized_sample_id": normalize_sample_id(sid),
                }
            )
    if not records:
        return pd.DataFrame(columns=["modality", "raw_sample_id", "normalized_sample_id"])
    return pd.DataFrame(records)
