"""Raw file integrity checks and inventory logging.

Input:
- paths/default configs
- raw files under data/raw and source files under experiment_data

Output:
- results/logs/raw_file_inventory.csv
- results/logs/raw_file_inventory.json
- results/logs/raw_copy_verification.csv

Purpose:
- Provide reproducible evidence that dimensions are based on full streamed scans.

TODO:
- Add optional parquet conversion suggestions based on scan metrics.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Any
import csv
import hashlib
import json

import pandas as pd

from src.utils.io_utils import ensure_dir


OMICS_KEYS = [
    "mutation",
    "cnv",
    "rna",
    "methylation",
    "methylation_probe_map",
    "mirna",
    "rppa",
    "clinical_matrix",
    "survival",
]


@dataclass
class RawFileInventoryRecord:
    key: str
    file_path: str
    original_source_path: str
    file_exists: bool
    file_size_bytes: int
    file_size_mb: float
    delimiter_guess: str
    encoding_guess: str
    total_lines: int
    header_column_count: int
    parsed_column_count: int
    read_mode: str
    full_scan_completed: bool
    preview_rows_captured: str
    preview_columns_captured: int
    notes: str


@dataclass
class RawCopyVerificationRecord:
    key: str
    source_path: str
    target_path: str
    source_exists: bool
    target_exists: bool
    source_size_bytes: int
    target_size_bytes: int
    size_match: bool
    source_hash: str
    target_hash: str
    hash_algorithm: str
    hash_match: bool
    notes: str


def _guess_encoding(path: Path, sample_bytes: int) -> str:
    with path.open("rb") as f:
        chunk = f.read(sample_bytes)
    candidates = ["utf-8-sig", "utf-8", "gb18030", "latin-1"]
    for enc in candidates:
        try:
            chunk.decode(enc)
            return enc
        except UnicodeDecodeError:
            continue
    return "unknown"


def _guess_delimiter(header_line: str) -> str:
    candidates = ["\t", ",", ";", "|"]
    counts = {d: header_line.count(d) for d in candidates}
    best = max(counts, key=counts.get)
    if counts[best] == 0:
        return "unknown"
    return "tab" if best == "\t" else best


def _split_line(line: str, delimiter_guess: str) -> list[str]:
    if delimiter_guess == "tab":
        return line.rstrip("\n\r").split("\t")
    if delimiter_guess in {",", ";", "|"}:
        return line.rstrip("\n\r").split(delimiter_guess)
    return line.rstrip("\n\r").split()


def _compute_hash(path: Path, hash_algorithm: str, chunk_bytes: int) -> str:
    h = hashlib.new(hash_algorithm)
    with path.open("rb") as f:
        while True:
            block = f.read(chunk_bytes)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def _scan_single_file(
    key: str,
    target_path: Path,
    source_path: Path,
    preview_rows: int,
    preview_columns: int,
    scan_encoding_bytes: int,
) -> RawFileInventoryRecord:
    exists = target_path.exists() and target_path.is_file()
    if not exists:
        return RawFileInventoryRecord(
            key=key,
            file_path=str(target_path),
            original_source_path=str(source_path),
            file_exists=False,
            file_size_bytes=0,
            file_size_mb=0.0,
            delimiter_guess="unknown",
            encoding_guess="unknown",
            total_lines=0,
            header_column_count=0,
            parsed_column_count=0,
            read_mode="header_preview",
            full_scan_completed=False,
            preview_rows_captured="[]",
            preview_columns_captured=0,
            notes="target file missing",
        )

    size_bytes = target_path.stat().st_size
    encoding_guess = _guess_encoding(target_path, scan_encoding_bytes)

    delimiter_guess = "unknown"
    total_lines = 0
    header_column_count = 0
    parsed_column_count = 0
    preview_cache: list[list[str]] = []
    col_counter: Counter[int] = Counter()
    mismatch_lines = 0
    full_scan_completed = False

    try:
        with target_path.open("r", encoding=encoding_guess if encoding_guess != "unknown" else "utf-8", errors="replace") as f:
            header = f.readline()
            if header:
                delimiter_guess = _guess_delimiter(header)
                header_columns = _split_line(header, delimiter_guess)
                header_column_count = len(header_columns)
                col_counter[header_column_count] += 1
                preview_cache.append(header_columns[:preview_columns])
                total_lines += 1

            for idx, line in enumerate(f, start=1):
                cols = _split_line(line, delimiter_guess)
                col_count = len(cols)
                col_counter[col_count] += 1
                if header_column_count > 0 and col_count != header_column_count:
                    mismatch_lines += 1
                if idx <= preview_rows:
                    preview_cache.append(cols[:preview_columns])
                total_lines += 1

        parsed_column_count = col_counter.most_common(1)[0][0] if col_counter else 0
        full_scan_completed = True
    except OSError as exc:
        return RawFileInventoryRecord(
            key=key,
            file_path=str(target_path),
            original_source_path=str(source_path),
            file_exists=True,
            file_size_bytes=size_bytes,
            file_size_mb=round(size_bytes / (1024 * 1024), 3),
            delimiter_guess=delimiter_guess,
            encoding_guess=encoding_guess,
            total_lines=total_lines,
            header_column_count=header_column_count,
            parsed_column_count=parsed_column_count,
            read_mode="streamed_scan",
            full_scan_completed=False,
            preview_rows_captured=json.dumps(preview_cache, ensure_ascii=False),
            preview_columns_captured=preview_columns,
            notes=f"scan failed: {exc}",
        )

    notes: list[str] = []
    if mismatch_lines > 0:
        notes.append(f"column_count_mismatch_lines={mismatch_lines}")
    if size_bytes > 200 * 1024 * 1024:
        notes.append("large_file_stream_scan_recommended")
    if not notes:
        notes.append("streamed full scan completed")

    return RawFileInventoryRecord(
        key=key,
        file_path=str(target_path),
        original_source_path=str(source_path),
        file_exists=True,
        file_size_bytes=size_bytes,
        file_size_mb=round(size_bytes / (1024 * 1024), 3),
        delimiter_guess=delimiter_guess,
        encoding_guess=encoding_guess,
        total_lines=total_lines,
        header_column_count=header_column_count,
        parsed_column_count=parsed_column_count,
        read_mode="streamed_scan",
        full_scan_completed=full_scan_completed,
        preview_rows_captured=json.dumps(preview_cache, ensure_ascii=False),
        preview_columns_captured=preview_columns,
        notes="; ".join(notes),
    )


def run_raw_integrity_checks(paths_cfg: dict[str, Any], default_cfg: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    integrity_cfg = default_cfg.get("data_integrity", {})
    preview_rows = int(integrity_cfg.get("preview_rows", 2))
    preview_columns = int(integrity_cfg.get("preview_columns", 12))
    scan_encoding_bytes = int(integrity_cfg.get("scan_encoding_bytes", 65536))
    hash_algorithm = str(integrity_cfg.get("hash_algorithm", "sha256"))
    hash_chunk_bytes = int(integrity_cfg.get("hash_chunk_bytes", 4 * 1024 * 1024))

    logs_dir = Path(paths_cfg["results"]["logs"])
    ensure_dir(logs_dir)

    inventory_records: list[RawFileInventoryRecord] = []
    verify_records: list[RawCopyVerificationRecord] = []

    for key in OMICS_KEYS:
        target = Path(paths_cfg["data"][key])
        source = Path(paths_cfg.get("source_data", {}).get(key, ""))

        inventory_records.append(
            _scan_single_file(
                key=key,
                target_path=target,
                source_path=source,
                preview_rows=preview_rows,
                preview_columns=preview_columns,
                scan_encoding_bytes=scan_encoding_bytes,
            )
        )

        source_exists = source.exists() and source.is_file()
        target_exists = target.exists() and target.is_file()
        source_size = source.stat().st_size if source_exists else 0
        target_size = target.stat().st_size if target_exists else 0
        size_match = source_size == target_size if source_exists and target_exists else False

        source_hash = _compute_hash(source, hash_algorithm, hash_chunk_bytes) if source_exists else ""
        target_hash = _compute_hash(target, hash_algorithm, hash_chunk_bytes) if target_exists else ""
        hash_match = source_hash == target_hash if source_exists and target_exists else False

        notes = []
        if not source_exists:
            notes.append("source_missing")
        if not target_exists:
            notes.append("target_missing")
        if source_exists and target_exists and not size_match:
            notes.append("size_mismatch")
        if source_exists and target_exists and not hash_match:
            notes.append("hash_mismatch")
        if not notes:
            notes.append("copy_verified")

        verify_records.append(
            RawCopyVerificationRecord(
                key=key,
                source_path=str(source),
                target_path=str(target),
                source_exists=source_exists,
                target_exists=target_exists,
                source_size_bytes=source_size,
                target_size_bytes=target_size,
                size_match=size_match,
                source_hash=source_hash,
                target_hash=target_hash,
                hash_algorithm=hash_algorithm,
                hash_match=hash_match,
                notes="; ".join(notes),
            )
        )

    inventory_df = pd.DataFrame([asdict(x) for x in inventory_records])
    verification_df = pd.DataFrame([asdict(x) for x in verify_records])

    inventory_csv = Path(paths_cfg["results"].get("raw_file_inventory_csv", logs_dir / "raw_file_inventory.csv"))
    inventory_json = Path(paths_cfg["results"].get("raw_file_inventory_json", logs_dir / "raw_file_inventory.json"))
    verify_csv = Path(paths_cfg["results"].get("raw_copy_verification_csv", logs_dir / "raw_copy_verification.csv"))

    inventory_df.to_csv(inventory_csv, index=False, encoding="utf-8")
    verification_df.to_csv(verify_csv, index=False, encoding="utf-8")

    payload = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "inventory": [asdict(x) for x in inventory_records],
        "copy_verification": [asdict(x) for x in verify_records],
    }
    with inventory_json.open("w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

    return inventory_df, verification_df, payload


def export_integrity_markdown(
    inventory_df: pd.DataFrame,
    verification_df: pd.DataFrame,
    out_path: Path,
) -> None:
    ensure_dir(out_path.parent)

    lines: list[str] = []
    lines.append("# 11 数据读取完整性校验")
    lines.append("")
    lines.append("本文件基于 results/logs 中的结构化日志自动生成。")
    lines.append("")
    lines.append("- 说明：\"头部预览\"仅用于识别文件结构，不代表已完成全量读取。")
    lines.append("- 说明：\"已完成全量扫描\"表示该文件已通过流式方式完整遍历，维度统计可作为后续文档依据。")
    lines.append("- 说明：若某文件尚未完成完整解析，当前维度仅为初步统计，待后续预处理阶段进一步确认。")
    lines.append("")

    lines.append("## 1. 原始文件读取完整性总表")
    lines.append("")
    lines.append("| key | file_size_mb | delimiter_guess | encoding_guess | total_lines | header_column_count | parsed_column_count | read_mode | full_scan_completed | notes |")
    lines.append("|---|---:|---|---|---:|---:|---:|---|---|---|")
    for _, row in inventory_df.iterrows():
        lines.append(
            f"| {row['key']} | {row['file_size_mb']} | {row['delimiter_guess']} | {row['encoding_guess']} | "
            f"{row['total_lines']} | {row['header_column_count']} | {row['parsed_column_count']} | "
            f"{row['read_mode']} | {row['full_scan_completed']} | {row['notes']} |"
        )
    lines.append("")

    lines.append("## 2. 文件复制一致性校验")
    lines.append("")
    lines.append("| key | size_match | hash_match | hash_algorithm | notes |")
    lines.append("|---|---|---|---|---|")
    for _, row in verification_df.iterrows():
        lines.append(
            f"| {row['key']} | {row['size_match']} | {row['hash_match']} | {row['hash_algorithm']} | {row['notes']} |"
        )
    lines.append("")

    lines.append("## 3. 风险点与后续建议")
    lines.append("")
    lines.append("- 对超大矩阵（如 methylation）避免默认一次性全量 read_csv，优先使用 chunked_read 或 streamed_scan。")
    lines.append("- 后续预处理脚本需记录输入维度、读取模式（full/chunked）、输出维度、过滤数量与原因。")
    lines.append("- 若出现列数不一致，需在预处理阶段先做格式修复再进入建模。")
    lines.append("- 可评估将超大中间结果转存为 parquet 以降低重复 I/O 成本。")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_preprocess_dimension_change(
    log_path: Path,
    modality: str,
    stage: str,
    read_mode: str,
    input_rows: int,
    input_cols: int,
    output_rows: int,
    output_cols: int,
    filtered_rows: int,
    filtered_cols: int,
    reason: str,
) -> None:
    ensure_dir(log_path.parent)
    fieldnames = [
        "timestamp",
        "modality",
        "stage",
        "read_mode",
        "input_rows",
        "input_cols",
        "output_rows",
        "output_cols",
        "filtered_rows",
        "filtered_cols",
        "reason",
    ]
    file_exists = log_path.exists()
    with log_path.open("a", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()
        writer.writerow(
            {
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "modality": modality,
                "stage": stage,
                "read_mode": read_mode,
                "input_rows": input_rows,
                "input_cols": input_cols,
                "output_rows": output_rows,
                "output_cols": output_cols,
                "filtered_rows": filtered_rows,
                "filtered_cols": filtered_cols,
                "reason": reason,
            }
        )
