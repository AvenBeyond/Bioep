from pathlib import Path

from main import run_status_check


def test_main_smoke_runs():
    root = Path(__file__).resolve().parents[1]
    exit_code = run_status_check(root)
    assert exit_code == 0
