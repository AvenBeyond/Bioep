from pathlib import Path

import pandas as pd

from src.feature_engineering.build_similarity_matrices import build_similarity_round1


def test_similarity_build_smoke():
    root = Path(__file__).resolve().parents[1]
    out = build_similarity_round1(root, modalities=["mutation", "cnv", "rna", "mirna"])

    assert len(out["cohort_samples"]) > 0
    assert set(out["modalities"]) == {"mutation", "cnv", "rna", "mirna"}

    shape_path = root / "results" / "tables" / "modality_modeling_shapes.csv"
    log_path = root / "results" / "tables" / "similarity_build_log.csv"
    assert shape_path.exists()
    assert log_path.exists()

    shape_df = pd.read_csv(shape_path)
    assert not shape_df.empty
    assert (shape_df["aligned_samples"] > 0).all()
