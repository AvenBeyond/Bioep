from pathlib import Path

import numpy as np
import pandas as pd

from src.clustering.proposed_weighted_fusion import (
    build_proposed_weight_components,
    estimate_final_weights,
    run_weighted_fusion_clustering,
)


def test_weighted_fusion_smoke(tmp_path: Path):
    samples = [f"S{i}" for i in range(12)]
    rng = np.random.default_rng(42)

    def mk_aff():
        a = rng.random((12, 12))
        a = (a + a.T) / 2
        np.fill_diagonal(a, 1.0)
        return pd.DataFrame(a, index=samples, columns=samples)

    affs = {m: mk_aff() for m in ["mutation", "cnv", "methylation", "rna", "mirna"]}
    x = pd.DataFrame(rng.normal(size=(12, 30)), index=samples)

    comp = build_proposed_weight_components(
        modalities=list(affs.keys()),
        q_scores={m: 0.5 for m in affs},
        c_scores={m: 0.4 for m in affs},
        s_scores={m: 1.0 for m in affs},
        m_scores={m: 0.1 for m in affs},
    )
    w = estimate_final_weights(comp, alpha=1.0, beta=1.0, delta=1.0, gamma=1.0)

    out = run_weighted_fusion_clustering(
        affinity_mats=affs,
        weight_table=w,
        combined_feature_matrix=x,
        k_values=[2, 3],
        modality_set_name="main5",
        out_tables_dir=tmp_path,
        out_figures_dir=tmp_path,
        label_prefix="smoke_proposed",
        random_state=42,
    )

    assert not out["metrics"].empty
    assert (tmp_path / "smoke_proposed_labels_k2.csv").exists()
    assert (tmp_path / "proposed_weight_barplot.png").exists()
