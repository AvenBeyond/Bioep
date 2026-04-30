from pathlib import Path

from src.utils.config_utils import load_project_config


def test_config_loading_success():
    root = Path(__file__).resolve().parents[1]
    paths_cfg, default_cfg = load_project_config(root)
    assert "data" in paths_cfg
    assert "results" in paths_cfg
    assert default_cfg["random_seed"] == 42
    assert default_cfg["k_range"] == [2, 3, 4, 5, 6]
