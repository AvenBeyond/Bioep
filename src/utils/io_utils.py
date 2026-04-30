"""I/O utilities.

Input:
- Path objects for files/directories

Output:
- File existence checks and lightweight previews

Purpose:
- Keep file system helper logic in one place.

TODO:
- Add chunked readers for very large omics files.
"""

from __future__ import annotations

from pathlib import Path


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def file_exists(path: Path) -> bool:
    return path.exists() and path.is_file()


def read_head(path: Path, n_lines: int = 2) -> list[str]:
    lines: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for _ in range(n_lines):
            line = f.readline()
            if not line:
                break
            lines.append(line.rstrip("\n"))
    return lines
