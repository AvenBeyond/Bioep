"""Project main entrypoint.

Input:
- configs/paths.yaml
- configs/default.yaml

Output:
- status summary in terminal
- inventory CSV in results/tables

Purpose:
- first-round executable check: config loading, path validation, raw data presence.

TODO:
- wire full preprocessing, clustering, and evaluation pipeline.
"""

from __future__ import annotations

from pathlib import Path
import sys

from src.data_loading.load_all_modalities import build_file_inventory, build_sample_inventory
from src.data_loading.raw_integrity_checks import export_integrity_markdown, run_raw_integrity_checks
from src.utils.config_utils import load_project_config
from src.utils.io_utils import ensure_dir
from src.utils.logging_utils import get_logger


def run_status_check(project_root: Path) -> int:
    paths_cfg, default_cfg = load_project_config(project_root)
    logger = get_logger("main", Path(paths_cfg["results"]["logs"]))

    ensure_dir(Path(paths_cfg["results"]["tables"]))
    inventory_df = build_file_inventory(paths_cfg)
    sample_inv_df = build_sample_inventory(paths_cfg)

    inventory_path = Path(paths_cfg["results"]["tables"]) / "raw_data_inventory.csv"
    sample_inv_path = Path(paths_cfg["results"]["tables"]) / "sample_inventory.csv"
    inventory_df.to_csv(inventory_path, index=False, encoding="utf-8")
    sample_inv_df.to_csv(sample_inv_path, index=False, encoding="utf-8")

    existing = int(inventory_df["exists"].sum())
    total = int(inventory_df.shape[0])
    logger.info("Project: %s", default_cfg.get("project_name", "unknown"))
    logger.info("Raw data files found: %s/%s", existing, total)
    logger.info("Inventory saved to: %s", inventory_path)
    logger.info("Sample inventory saved to: %s", sample_inv_path)

    print("=== TCGA-STAD Project Status ===")
    print(f"Project root: {project_root}")
    print(f"Raw data files found: {existing}/{total}")
    print(f"Inventory: {inventory_path}")
    print(f"Sample inventory: {sample_inv_path}")

    return 0 if existing == total else 1


def run_integrity_audit(project_root: Path) -> int:
    paths_cfg, default_cfg = load_project_config(project_root)
    logger = get_logger("integrity_audit", Path(paths_cfg["results"]["logs"]))

    inventory_df, verification_df, _ = run_raw_integrity_checks(paths_cfg, default_cfg)
    export_integrity_markdown(inventory_df, verification_df, project_root / "docs" / "11_data_reading_integrity_check.md")

    full_scanned = int(inventory_df["full_scan_completed"].sum())
    total = int(inventory_df.shape[0])
    copy_ok = int(verification_df["hash_match"].sum())

    logger.info("Raw integrity scan completed: %s/%s files fully scanned", full_scanned, total)
    logger.info("Raw copy verification passed: %s/%s files", copy_ok, total)
    print("=== Raw Data Integrity Audit ===")
    print(f"Fully scanned files: {full_scanned}/{total}")
    print(f"Copy verification hash match: {copy_ok}/{total}")
    print(f"Inventory CSV: {paths_cfg['results']['raw_file_inventory_csv']}")
    print(f"Inventory JSON: {paths_cfg['results']['raw_file_inventory_json']}")
    print(f"Copy verification CSV: {paths_cfg['results']['raw_copy_verification_csv']}")
    print("Integrity doc: docs/11_data_reading_integrity_check.md")

    return 0 if full_scanned == total and copy_ok == total else 1


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    if len(sys.argv) > 1 and sys.argv[1] == "--integrity":
        sys.exit(run_integrity_audit(root))
    sys.exit(run_status_check(root))
