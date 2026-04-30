"""Run round-1 preprocessing for all modalities.

Input:
- raw modality files
- sample_inventory / overlap analysis outputs for cohort selection

Output:
- interim modality files
- results/logs/preprocessing_dimension_changes.csv
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.preprocessing.preprocess_cnv import run_cnv_round1
from src.preprocessing.preprocess_methylation import run_methylation_round1
from src.preprocessing.preprocess_mirna import run_mirna_round1
from src.preprocessing.preprocess_mutation import run_mutation_round1
from src.preprocessing.preprocess_rna import run_rna_round1
from src.preprocessing.preprocess_rppa import run_rppa_round1
from src.utils.config_utils import load_project_config
from src.utils.io_utils import ensure_dir


MAIN_MODALITIES = ["mutation", "cnv", "methylation", "rna", "mirna"]


def load_main_intersection_samples(project_root: Path) -> set[str]:
    sample_inventory_path = project_root / "results" / "tables" / "sample_inventory.csv"
    if not sample_inventory_path.exists():
        raise FileNotFoundError("sample_inventory.csv not found. Run Phase 7 sample_alignment_analysis first.")

    inv = pd.read_csv(sample_inventory_path)
    modality_sets = {
        m: set(inv.loc[inv["modality"] == m, "normalized_sample_id"].dropna().astype(str).unique())
        for m in MAIN_MODALITIES
    }
    return set.intersection(*modality_sets.values())


def run_round1(project_root: Path) -> None:
    paths_cfg, _ = load_project_config(project_root)
    data_cfg = paths_cfg["data"]

    interim_root = Path(data_cfg["interim_root"])
    ensure_dir(interim_root)

    logs_dir = Path(paths_cfg["results"]["logs"])
    ensure_dir(logs_dir)
    preprocess_log = logs_dir / "preprocessing_dimension_changes.csv"

    main_samples = load_main_intersection_samples(project_root)

    run_mutation_round1(
        input_path=Path(data_cfg["mutation"]),
        output_path=interim_root / "mutation_round1.csv",
        min_event_rate=0.02,
        selected_samples=main_samples,
        log_path=preprocess_log,
    )
    run_cnv_round1(
        input_path=Path(data_cfg["cnv"]),
        output_path=interim_root / "cnv_round1.csv",
        min_variance=0.01,
        selected_samples=main_samples,
        log_path=preprocess_log,
    )
    run_rna_round1(
        input_path=Path(data_cfg["rna"]),
        output_path=interim_root / "rna_round1.csv",
        min_nonzero_rate=0.1,
        top_var_features=5000,
        selected_samples=main_samples,
        log_path=preprocess_log,
    )
    run_mirna_round1(
        input_path=Path(data_cfg["mirna"]),
        output_path=interim_root / "mirna_round1.csv",
        min_nonzero_rate=0.1,
        top_var_features=1000,
        selected_samples=main_samples,
        log_path=preprocess_log,
    )
    run_methylation_round1(
        input_path=Path(data_cfg["methylation"]),
        summary_output_path=interim_root / "methylation_round1_summary.csv",
        preview_output_path=interim_root / "methylation_round1_preview_topvar.csv",
        selected_samples=main_samples,
        probe_map_path=Path(data_cfg["methylation_probe_map"]),
        top_var_preview=2000,
        preview_scan_rows=50000,
        chunk_size=5000,
        log_path=preprocess_log,
    )

    # RPPA kept as supplementary track, uses its own sample set.
    run_rppa_round1(
        input_path=Path(data_cfg["rppa"]),
        output_path=interim_root / "rppa_round1.csv",
        min_nonzero_rate=0.1,
        selected_samples=None,
        log_path=preprocess_log,
    )


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[2]
    run_round1(root)
    print("Round-1 preprocessing completed.")
