from pathlib import Path

import pandas as pd

from src.evaluation.ablation_study import save_ablation_plots, write_ablation_outputs


def test_ablation_smoke(tmp_path: Path):
    df = pd.DataFrame(
        [
            {
                "experiment_name": "early_fusion_kmeans",
                "modality_set": "main5",
                "k": 2,
                "silhouette": 0.11,
                "stability_score": 0.70,
                "logrank_os_p": 0.04,
                "logrank_pfi_p": 0.10,
                "cluster_balance_score": 0.9,
                "notes": "smoke",
            },
            {
                "experiment_name": "weight_full",
                "modality_set": "main5",
                "k": 2,
                "silhouette": 0.12,
                "stability_score": 0.65,
                "logrank_os_p": 0.05,
                "logrank_pfi_p": 0.20,
                "cluster_balance_score": 0.8,
                "notes": "smoke",
            },
        ]
    )

    ranked, summary = write_ablation_outputs(df, tmp_path)
    save_ablation_plots(ranked, tmp_path)

    assert not ranked.empty
    assert not summary.empty
    assert (tmp_path / "ablation_results.csv").exists()
    assert (tmp_path / "ablation_comparison.png").exists()
