from pathlib import Path

import numpy as np
import pandas as pd

from src.evaluation.subtype_classifier import run_subtype_classifier_internal_validation, save_classifier_plots


def test_classifier_smoke(tmp_path: Path):
    rng = np.random.default_rng(42)
    n = 60
    x = pd.DataFrame(rng.normal(size=(n, 20)), index=[f"S{i}" for i in range(n)], columns=[f"f{j}" for j in range(20)])
    y = pd.Series([0] * 20 + [1] * 20 + [2] * 20, index=x.index)

    out = run_subtype_classifier_internal_validation(
        features=x,
        labels=y,
        cv_folds=3,
        cv_repeats=1,
        random_state=42,
    )
    save_classifier_plots(out["cv_results"], out["confusion"], out["feature_importance"], tmp_path)

    assert not out["cv_results"].empty
    assert (tmp_path / "subtype_classifier_performance.png").exists()
