"""Configuration utilities.

Input:
- YAML files in configs/

Output:
- Parsed Python dictionaries with absolute pathlib paths

Purpose:
- Centralize configuration loading and path resolution.

TODO:
- Add schema validation and richer runtime checks.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def resolve_paths(paths_cfg: dict[str, Any], project_root: Path) -> dict[str, Any]:
    resolved = {}
    for key, value in paths_cfg.items():
        if isinstance(value, dict):
            resolved[key] = resolve_paths(value, project_root)
        elif isinstance(value, str):
            resolved[key] = (project_root / value).resolve()
        else:
            resolved[key] = value
    return resolved


def load_project_config(project_root: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    paths_cfg = load_yaml(project_root / "configs" / "paths.yaml")
    default_cfg = load_yaml(project_root / "configs" / "default.yaml")
    resolved_paths = resolve_paths(paths_cfg, project_root)
    return resolved_paths, default_cfg
