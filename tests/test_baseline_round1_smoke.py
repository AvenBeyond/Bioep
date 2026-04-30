from pathlib import Path

import pandas as pd

from src.pipelines.run_phase11_15 import run_phase11_15


def test_baseline_round1_pipeline_smoke():
    root = Path(__file__).resolve().parents[1]
    run_phase11_15(root)

    metrics_path = root / "results" / "tables" / "baseline_round1_metrics.csv"
    sizes_path = root / "results" / "tables" / "baseline_round1_cluster_sizes.csv"
    stability_path = root / "results" / "tables" / "cluster_stability_round1.csv"

    assert metrics_path.exists()
    assert sizes_path.exists()
    assert stability_path.exists()

    metrics = pd.read_csv(metrics_path)
    assert not metrics.empty
    assert metrics["k"].isin([2, 3, 4, 5, 6]).all()
