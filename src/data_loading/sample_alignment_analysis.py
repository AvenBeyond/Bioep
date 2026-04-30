"""Sample alignment and overlap analysis for multi-omics STAD data.

Input:
- data/raw modality files resolved from configs/paths.yaml

Output:
- results/tables/modality_dimension_summary.csv
- results/tables/sample_inventory.csv
- results/tables/sample_overlap_matrix.csv
- results/tables/multiomics_intersection_summary.csv
- results/tables/patient_alignment_rules.csv
- results/figures/sample_overlap_heatmap.png
- results/figures/modality_sample_counts_barplot.png
- results/figures/main_vs_with_rppa_intersection.png

Purpose:
- Phase 7 data truth: parse sample IDs, normalize TCGA IDs, compute overlap statistics.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from src.utils.config_utils import load_project_config
from src.utils.io_utils import ensure_dir


OVERLAP_MODALITIES = ["mutation", "cnv", "methylation", "rna", "mirna", "rppa"]
MAIN_MODALITIES = ["mutation", "cnv", "methylation", "rna", "mirna"]


@dataclass
class ModalityMeta:
    modality: str
    raw_file_name: str
    raw_file_path: str
    sample_axis_orientation: str
    raw_sample_count: int
    raw_feature_count: int
    parsed_sample_id_example: str
    parsed_feature_id_example: str
    value_type: str
    delimiter: str
    missingness_summary: str
    notes: str


def normalize_sample_id(sample_id: str) -> str:
    sid = str(sample_id).strip().upper()
    if sid == "" or sid in {"NA", "NAN"}:
        return ""
    return sid[:16] if len(sid) >= 16 else sid


def normalize_patient_id(sample_id: str) -> str:
    sid = normalize_sample_id(sample_id)
    return sid[:12] if len(sid) >= 12 else sid


def parse_sample_type_candidate(sample_id: str) -> str:
    sid = normalize_sample_id(sample_id)
    if len(sid) < 15:
        return "Unknown"
    sample_type_code = sid[13:15]
    if sample_type_code == "01":
        return "Primary Tumor"
    if sample_type_code == "11":
        return "Solid Tissue Normal"
    if sample_type_code == "06":
        return "Metastatic"
    if sample_type_code == "02":
        return "Recurrent Tumor"
    return f"Other({sample_type_code})"


def _read_header_and_examples(path: Path) -> tuple[list[str], list[str], str]:
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header_line = f.readline().rstrip("\n\r")
        first_data_line = f.readline().rstrip("\n\r")
    delimiter = "\t" if "\t" in header_line else ","
    header = header_line.split(delimiter)
    first_data = first_data_line.split(delimiter) if first_data_line else []
    return header, first_data, delimiter


def _estimate_missingness(path: Path, sample_rows: int = 20000) -> str:
    # Fast estimate for large files: evaluate first N data rows only.
    try:
        df = pd.read_csv(path, sep="\t", nrows=sample_rows, dtype=str)
        if df.shape[1] <= 1:
            return "n/a"
        value_df = df.iloc[:, 1:]
        miss = value_df.isna().mean().mean() * 100
        return f"estimated_missing={miss:.3f}%_from_first_{min(sample_rows, df.shape[0])}_rows"
    except (OSError, UnicodeDecodeError, ValueError, pd.errors.ParserError) as exc:
        return f"missingness_estimation_failed:{type(exc).__name__}"


def _value_type_for_modality(modality: str) -> str:
    mapping = {
        "mutation": "binary_or_event_indicator",
        "cnv": "discrete_gistic_-2_to_2",
        "methylation": "continuous_beta_0_to_1",
        "rna": "continuous_expression",
        "mirna": "continuous_expression_with_na",
        "rppa": "continuous_protein_expression",
    }
    return mapping.get(modality, "unknown")


def build_phase7_outputs(project_root: Path) -> dict[str, Path]:
    paths_cfg, _ = load_project_config(project_root)

    tables_dir = Path(paths_cfg["results"]["tables"])
    figures_dir = Path(paths_cfg["results"]["figures"])
    ensure_dir(tables_dir)
    ensure_dir(figures_dir)

    modality_rows: list[ModalityMeta] = []
    sample_inventory_rows: list[dict[str, Any]] = []
    modality_to_samples: dict[str, set[str]] = {}

    for modality in OVERLAP_MODALITIES:
        path = Path(paths_cfg["data"][modality])
        header, first_data, delimiter = _read_header_and_examples(path)
        delimiter_name = "tab" if delimiter == "\t" else "comma"

        sample_ids = [x for x in header[1:] if str(x).strip()]
        normalized_samples = [normalize_sample_id(x) for x in sample_ids if normalize_sample_id(x)]
        modality_to_samples[modality] = set(normalized_samples)

        raw_feature_count = max(0, sum(1 for _ in path.open("r", encoding="utf-8", errors="replace")) - 1)
        sample_example = sample_ids[0] if sample_ids else ""
        feature_example = first_data[0] if first_data else ""

        modality_rows.append(
            ModalityMeta(
                modality=modality,
                raw_file_name=path.name,
                raw_file_path=str(path),
                sample_axis_orientation="column",
                raw_sample_count=len(sample_ids),
                raw_feature_count=raw_feature_count,
                parsed_sample_id_example=sample_example,
                parsed_feature_id_example=feature_example,
                value_type=_value_type_for_modality(modality),
                delimiter=delimiter_name,
                missingness_summary=_estimate_missingness(path),
                notes="sample_ids parsed from header columns",
            )
        )

        for raw_sid in sample_ids:
            nsid = normalize_sample_id(raw_sid)
            npid = normalize_patient_id(raw_sid)
            stype = parse_sample_type_candidate(raw_sid)
            sample_inventory_rows.append(
                {
                    "modality": modality,
                    "raw_sample_id": raw_sid,
                    "normalized_sample_id": nsid,
                    "normalized_patient_id": npid,
                    "sample_type_candidate": stype,
                    "included_in_main_analysis_candidate": bool(modality in MAIN_MODALITIES and stype == "Primary Tumor"),
                    "notes": "normalized_by_tcga_barcode_rules",
                }
            )

    dimension_df = pd.DataFrame([m.__dict__ for m in modality_rows])
    sample_inventory_df = pd.DataFrame(sample_inventory_rows)

    overlap_df = pd.DataFrame(index=OVERLAP_MODALITIES, columns=OVERLAP_MODALITIES, dtype=int)
    for m1 in OVERLAP_MODALITIES:
        for m2 in OVERLAP_MODALITIES:
            overlap_df.loc[m1, m2] = len(modality_to_samples[m1].intersection(modality_to_samples[m2]))

    all_samples = set().union(*[modality_to_samples[m] for m in OVERLAP_MODALITIES])
    sample_presence: dict[str, list[str]] = {}
    for sid in all_samples:
        present = [m for m in OVERLAP_MODALITIES if sid in modality_to_samples[m]]
        sample_presence[sid] = present

    only_counts = {
        f"{m}_only_sample_count": sum(1 for mods in sample_presence.values() if mods == [m])
        for m in OVERLAP_MODALITIES
    }
    main_intersection = set.intersection(*[modality_to_samples[m] for m in MAIN_MODALITIES])
    all6_intersection = set.intersection(*[modality_to_samples[m] for m in OVERLAP_MODALITIES])

    ge3_main = sum(1 for sid in all_samples if sum(1 for m in MAIN_MODALITIES if sid in modality_to_samples[m]) >= 3)
    ge4_main = sum(1 for sid in all_samples if sum(1 for m in MAIN_MODALITIES if sid in modality_to_samples[m]) >= 4)
    ge5_main = len(main_intersection)

    intersection_summary_df = pd.DataFrame(
        [
            {
                **only_counts,
                "main5_intersection_sample_count": len(main_intersection),
                "all6_with_rppa_intersection_sample_count": len(all6_intersection),
                "partial_fusion_main5_ge3_modalities_sample_count": ge3_main,
                "partial_fusion_main5_ge4_modalities_sample_count": ge4_main,
                "partial_fusion_main5_ge5_modalities_sample_count": ge5_main,
            }
        ]
    )

    alignment_rules_rows = []
    for modality in OVERLAP_MODALITIES:
        sample_example = next(iter(modality_to_samples[modality])) if modality_to_samples[modality] else ""
        alignment_rules_rows.append(
            {
                "modality": modality,
                "raw_sample_id_example": sample_example,
                "normalization_rule": "uppercase + keep first 16 chars as normalized_sample_id",
                "truncate_to_patient_level": True,
                "patient_level_rule": "first 12 chars as normalized_patient_id",
                "potential_ambiguity": "same patient may have multiple aliquots/sample types",
                "final_adopted_rule": "analysis primary key=normalized_sample_id; patient summaries use normalized_patient_id",
            }
        )
    alignment_rules_df = pd.DataFrame(alignment_rules_rows)

    # Save tables
    dimension_path = tables_dir / "modality_dimension_summary.csv"
    sample_inventory_path = tables_dir / "sample_inventory.csv"
    overlap_path = tables_dir / "sample_overlap_matrix.csv"
    intersection_path = tables_dir / "multiomics_intersection_summary.csv"
    rules_path = tables_dir / "patient_alignment_rules.csv"

    dimension_df.to_csv(dimension_path, index=False, encoding="utf-8")
    sample_inventory_df.to_csv(sample_inventory_path, index=False, encoding="utf-8")
    overlap_df.to_csv(overlap_path, encoding="utf-8")
    intersection_summary_df.to_csv(intersection_path, index=False, encoding="utf-8")
    alignment_rules_df.to_csv(rules_path, index=False, encoding="utf-8")

    # Figures
    sns.set_theme(style="whitegrid")

    plt.figure(figsize=(8, 6))
    sns.heatmap(overlap_df.astype(int), annot=True, fmt="d", cmap="YlGnBu")
    plt.title("TCGA-STAD Multi-omics Sample Overlap Matrix")
    plt.xlabel("Modality")
    plt.ylabel("Modality")
    heatmap_path = figures_dir / "sample_overlap_heatmap.png"
    plt.tight_layout()
    plt.savefig(heatmap_path, dpi=220)
    plt.close()

    plt.figure(figsize=(8, 5))
    bar_df = dimension_df[["modality", "raw_sample_count"]].sort_values("raw_sample_count", ascending=False)
    sns.barplot(data=bar_df, x="modality", y="raw_sample_count", hue="modality", dodge=False, legend=False)
    plt.title("Raw Sample Counts by Modality")
    plt.xlabel("Modality")
    plt.ylabel("Sample Count")
    barplot_path = figures_dir / "modality_sample_counts_barplot.png"
    plt.tight_layout()
    plt.savefig(barplot_path, dpi=220)
    plt.close()

    plt.figure(figsize=(7, 5))
    compare_df = pd.DataFrame(
        {
            "cohort": ["Main 5 modalities", "Main 5 + RPPA"],
            "intersection_sample_count": [len(main_intersection), len(all6_intersection)],
        }
    )
    sns.barplot(data=compare_df, x="cohort", y="intersection_sample_count", hue="cohort", dodge=False, legend=False)
    plt.title("Intersection Size: Main 5 vs With RPPA")
    plt.xlabel("Cohort Definition")
    plt.ylabel("Intersection Sample Count")
    compare_path = figures_dir / "main_vs_with_rppa_intersection.png"
    plt.tight_layout()
    plt.savefig(compare_path, dpi=220)
    plt.close()

    return {
        "modality_dimension_summary": dimension_path,
        "sample_inventory": sample_inventory_path,
        "sample_overlap_matrix": overlap_path,
        "multiomics_intersection_summary": intersection_path,
        "patient_alignment_rules": rules_path,
        "sample_overlap_heatmap": heatmap_path,
        "modality_sample_counts_barplot": barplot_path,
        "main_vs_with_rppa_intersection": compare_path,
    }


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[2]
    outputs = build_phase7_outputs(root)
    print("Phase 7 outputs generated:")
    for key, value in outputs.items():
        print(f"- {key}: {value}")
