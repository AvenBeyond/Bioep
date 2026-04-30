from pathlib import Path

import pandas as pd

from src.data_loading.load_all_modalities import build_file_inventory
from src.utils.config_utils import load_project_config


def test_file_inventory_has_all_expected_rows():
    root = Path(__file__).resolve().parents[1]
    paths_cfg, _ = load_project_config(root)
    inventory = build_file_inventory(paths_cfg)
    assert inventory.shape[0] == 9
    assert inventory["exists"].all()


def test_raw_integrity_log_exists_and_complete():
    root = Path(__file__).resolve().parents[1]
    inventory_path = root / "results" / "logs" / "raw_file_inventory.csv"
    assert inventory_path.exists(), "Missing raw_file_inventory.csv; run `python main.py --integrity` first."

    df = pd.read_csv(inventory_path)
    required_columns = {
        "file_path",
        "original_source_path",
        "file_exists",
        "file_size_bytes",
        "file_size_mb",
        "delimiter_guess",
        "encoding_guess",
        "total_lines",
        "header_column_count",
        "parsed_column_count",
        "read_mode",
        "full_scan_completed",
        "preview_rows_captured",
        "preview_columns_captured",
        "notes",
    }
    assert required_columns.issubset(set(df.columns))
    assert int(df.shape[0]) == 9
    assert df["file_exists"].all()
    assert df["full_scan_completed"].all()
