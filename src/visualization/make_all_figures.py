"""Visualization skeleton.

Input: experiment outputs in results/tables
Output: publication/defense-ready figures in results/figures
Purpose: central place for all figure generation hooks.
TODO: implement each figure function and style templates.
"""

from __future__ import annotations

from pathlib import Path


def make_all_figures(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
