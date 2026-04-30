from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _file_type(path: Path) -> str:
    ext = path.suffix.lower().strip(".")
    if not ext:
        return "no_ext"
    # Some TCGA-style files include dots in the basename without true extension.
    if ext in {"csv", "tsv", "txt", "json", "yaml", "yml", "parquet", "feather", "pkl", "gz", "zip"}:
        return ext
    return "no_ext"


def _read_tsv_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        line = f.readline().rstrip("\n\r")
    return line.split("\t")


def _sample_type_code(sample_id: str) -> str:
    parts = str(sample_id).strip().upper().split("-")
    if len(parts) >= 4 and parts[3]:
        return parts[3][:2]
    return "NA"


def _modality_or_role(path: Path) -> str:
    s = str(path).lower().replace("\\", "/")
    if "mutation" in s or "mc3" in s:
        return "mutation"
    if "gistic" in s or "/cnv/" in s or "cnv_round" in s or "cnv_" in path.name.lower():
        return "cnv"
    if "methyl" in s or "probe" in s:
        return "methylation"
    if "mirna" in s:
        return "mirna"
    if "hiseqv2" in s or "/rna/" in s or "rna_round" in s or "rna_" in path.name.lower():
        return "rna"
    if "rppa" in s:
        return "rppa"
    if "survival" in s:
        return "survival"
    if "clinical" in s:
        return "clinical"
    if "similarity" in s:
        return "similarity"
    if "partial_fusion" in s:
        return "partial_fusion"
    if "classifier" in s or "subtype_assignment" in s:
        return "classifier"
    if "final_" in s:
        return "final_reporting"
    return "other"


def _asset_level(path: Path) -> str:
    s = str(path).lower().replace("\\", "/")
    if "/data/interim/" in s:
        return "interim"
    if "/results/tables/" in s:
        return "results"
    if "probemap" in s or "patient_alignment_rules" in s:
        return "annotation"
    if "/data/raw/" in s or "/experiment_data/" in s:
        return "raw"
    return "results"


def _usage_flags(path: Path) -> dict[str, Any]:
    rel = str(path).replace("\\", "/").lower()
    modality = _modality_or_role(path)
    is_raw = ("/data/raw/" in rel) or ("/experiment_data/" in rel)
    core_final_inputs = {
        "data/interim/mutation_round1.csv",
        "data/interim/cnv_round1.csv",
        "data/interim/methylation_round2_modeling.csv",
        "data/interim/rna_round1.csv",
        "data/interim/mirna_round1.csv",
        "results/tables/baseline_main5_labels_early_fusion_kmeans_k2.csv",
    }
    core_final_reports = {
        "results/tables/final_model_comparison_summary.csv",
        "results/tables/final_dashboard_master.csv",
        "results/tables/final_key_numbers_for_ppt.csv",
        "results/tables/final_innovation_model_report.csv",
        "results/tables/final_innovation_key_takeaways.csv",
    }
    classifier_core = {
        "results/tables/subtype_classifier_cv_results.csv",
        "results/tables/subtype_assignment_confusion_matrix.csv",
        "results/tables/subtype_classifier_feature_importance.csv",
    }

    used_main = any(x in rel for x in [
        "baseline_main5_metrics.csv",
        "proposed_main5_metrics.csv",
        "ablation_results.csv",
        "improvement_model_comparison.csv",
        "multiomics_intersection_summary.csv",
        "modality_modeling_shapes.csv",
        "methylation_processing_summary.csv",
    ]) or any(rel.endswith(x) for x in core_final_inputs | core_final_reports | classifier_core)

    if is_raw and modality in {"mutation", "cnv", "methylation", "rna", "mirna", "clinical", "survival"}:
        used_main = True

    used_final = any(rel.endswith(x) for x in core_final_inputs | core_final_reports)
    if is_raw and modality in {"mutation", "cnv", "methylation", "rna", "mirna", "survival"}:
        used_final = True

    used_classifier = any(rel.endswith(x) for x in classifier_core) or any(
        rel.endswith(x)
        for x in {
            "data/interim/mutation_round1.csv",
            "data/interim/cnv_round1.csv",
            "data/interim/methylation_round2_modeling.csv",
            "data/interim/rna_round1.csv",
            "data/interim/mirna_round1.csv",
            "results/tables/baseline_main5_labels_early_fusion_kmeans_k2.csv",
        }
    )
    if is_raw and modality in {"mutation", "cnv", "methylation", "rna", "mirna"}:
        used_classifier = True

    partially_used = any(x in rel for x in [
        "preview",
        "summary",
        "partial_fusion",
        "aligned_matrices",
        "proposed_weight_inputs_preview",
        "probe",
    ]) or ("rppa_round1" in rel)

    underused = any(x in rel for x in [
        "rppa",
        "probe",
        "survival_stad_survival",
        "clinicalmatrix",
        "aligned_matrices",
        "partial_fusion",
    ]) and (not used_final)

    return {
        "currently_used_in_mainline": bool(used_main),
        "currently_used_in_final_model": bool(used_final),
        "currently_used_in_classifier": bool(used_classifier),
        "partially_used": bool(partially_used),
        "unused_but_potentially_useful": bool(underused),
    }


def _why_and_next_use(path: Path) -> tuple[str, str]:
    rel = str(path).replace("\\", "/").lower()
    modality = _modality_or_role(path)

    if "probe" in rel:
        return (
            "probe 注释已用于覆盖率统计，但未系统转为 promoter/gene/region 聚合特征",
            "构建 methylation 的 promoter/gene-level 聚合特征与区域注释特征",
        )
    if "survival_stad_survival" in rel:
        return (
            "主线主要使用 OS/PFI，DSS/DFI 尚未进入正式终点评估",
            "将 DSS/DFI 纳入多终点一致性验证与敏感性分析",
        )
    if "clinicalmatrix" in rel:
        return (
            "临床变量当前仅使用少数字段（年龄/性别/分期等）",
            "扩展 grade/治疗相关/复发相关字段做临床关联与分层解释",
        )
    if "rppa" in rel:
        return (
            "RPPA 会明显缩小 complete-case 样本，主线未纳入 final",
            "将 RPPA 作为辅助视图或弱配对 late-fusion 证据层",
        )
    if "aligned_matrices" in rel:
        return (
            "已生成对齐矩阵，但在 final 管线中未复用",
            "将其作为后续特征工程和快速重现实验的直接输入",
        )
    if "partial_fusion" in rel:
        return (
            "partial-fusion 已有结果，但尚未作为下一轮正式扩容主线",
            "优先推进 at-least-4 扩容并建立与 strict complete-case 的统一比较框架",
        )
    if modality == "methylation" and ("round2" in rel or "filter_log" in rel):
        return (
            "目前以高方差 probe-level 选择为主，聚合特征尚未展开",
            "加入 probe-to-gene/promoter 聚合与注释驱动特征降噪",
        )
    if modality == "mutation":
        return (
            "当前主要为基因层事件矩阵，尚未形成 burden/pathway/driver 汇总特征",
            "追加 burden、driver-centric 与 pathway-level alteration 特征",
        )
    if modality == "cnv":
        return (
            "当前主要直接使用 gene-level GISTIC，路径级摘要较少",
            "增加 focal event count 与 pathway CNV summary",
        )
    if modality in {"rna", "mirna"}:
        return (
            "当前主要使用筛选后的表达特征，尚未加入功能聚合表示",
            "构建 pathway/module score 与 gene-set 聚合特征",
        )
    return ("已用于当前流程或仅作中间记录", "可按需复用，不作为首要新增项")


def _collect_assets(root: Path) -> pd.DataFrame:
    paths: list[Path] = []
    for base in [root / "data" / "raw", root / "experiment_data", root / "data" / "interim", root / "results" / "tables"]:
        if not base.exists():
            continue
        for p in base.rglob("*"):
            if p.is_file() and p.name != ".gitkeep":
                paths.append(p)

    rows: list[dict[str, Any]] = []
    for p in sorted(set(paths)):
        flags = _usage_flags(p)
        why, nxt = _why_and_next_use(p)
        rows.append(
            {
                "asset_level": _asset_level(p),
                "file_name": p.name,
                "file_path": str(p),
                "file_type": _file_type(p),
                "modality_or_role": _modality_or_role(p),
                **flags,
                "why_not_fully_used": why,
                "potential_next_use": nxt,
            }
        )

    df = pd.DataFrame(rows)
    return df.sort_values(["asset_level", "modality_or_role", "file_name"]).reset_index(drop=True)


def _build_sample_salvage(root: Path) -> pd.DataFrame:
    inter = pd.read_csv(root / "results" / "tables" / "multiomics_intersection_summary.csv")
    m = inter.iloc[0]
    complete5 = int(m["main5_intersection_sample_count"])
    ge4 = int(m["partial_fusion_main5_ge4_modalities_sample_count"])
    ge3 = int(m["partial_fusion_main5_ge3_modalities_sample_count"])
    all6 = int(m["all6_with_rppa_intersection_sample_count"])

    normals = {}
    omics_raw = {
        "mutation": root / "data" / "raw" / "mutation" / "STAD_mc3_gene_level.txt",
        "cnv": root / "data" / "raw" / "cnv" / "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",
        "methylation": root / "data" / "raw" / "methylation" / "HumanMethylation450",
        "rna": root / "data" / "raw" / "rna" / "HiSeqV2",
        "mirna": root / "data" / "raw" / "mirna" / "miRNA_HiSeq_gene",
        "rppa": root / "data" / "raw" / "rppa" / "RPPA",
    }
    for mod, p in omics_raw.items():
        if not p.exists():
            normals[mod] = 0
            continue
        hdr = _read_tsv_header(p)
        samples = hdr[1:]
        normals[mod] = sum(1 for s in samples if _sample_type_code(s) == "11")

    rows = [
        {
            "opportunity": "strict_main5_complete_case",
            "current_n": complete5,
            "excluded_vs_atleast4": 0,
            "recoverable_n": 0,
            "strategy": "current_baseline",
            "can_include_via": "N/A",
            "notes": "主线当前严格 complete-case",
        },
        {
            "opportunity": "salvage_to_partial_atleast4_main5",
            "current_n": complete5,
            "excluded_vs_atleast4": ge4 - complete5,
            "recoverable_n": ge4 - complete5,
            "strategy": "partial_fusion",
            "can_include_via": "at_least_4_modalities",
            "notes": "最稳妥扩容路线，样本缺失较少",
        },
        {
            "opportunity": "salvage_to_partial_atleast3_main5",
            "current_n": complete5,
            "excluded_vs_atleast4": ge3 - complete5,
            "recoverable_n": ge3 - complete5,
            "strategy": "partial_fusion",
            "can_include_via": "at_least_3_modalities",
            "notes": "扩容更大，但异质性更高",
        },
        {
            "opportunity": "main5_with_optional_rppa_auxiliary",
            "current_n": complete5,
            "excluded_vs_atleast4": complete5 - all6,
            "recoverable_n": complete5 - all6,
            "strategy": "reduced_modality_or_auxiliary_layer",
            "can_include_via": "keep_main5_and_use_rppa_as_optional",
            "notes": "避免强制六组学带来的样本损失",
        },
    ]

    for mod, n in normals.items():
        rows.append(
            {
                "opportunity": f"normal_samples_in_{mod}",
                "current_n": n,
                "excluded_vs_atleast4": np.nan,
                "recoverable_n": n,
                "strategy": "auxiliary_analysis",
                "can_include_via": "sanity_check_or_contrastive_visualization",
                "notes": "可用于肿瘤-正常对照可视化或特征筛选 sanity check",
            }
        )

    return pd.DataFrame(rows)


def _build_feature_opportunities() -> pd.DataFrame:
    rows = [
        {
            "modality": "methylation",
            "current_state": "probe-level high-variance top features",
            "opportunity": "promoter/gene-level aggregation + CpG region annotation",
            "expected_benefit": "提高稳定性与可解释性，降低噪声",
            "implementation_difficulty": "medium",
            "priority": "high",
            "priority_score": 9,
        },
        {
            "modality": "rna",
            "current_state": "gene expression selected features",
            "opportunity": "pathway/module score and gene set aggregation",
            "expected_benefit": "增强临床关联与泛化解释",
            "implementation_difficulty": "medium",
            "priority": "high",
            "priority_score": 8,
        },
        {
            "modality": "mirna",
            "current_state": "selected expression features",
            "opportunity": "miRNA family/module aggregation and pathway linkage",
            "expected_benefit": "提升信号一致性与生物解释",
            "implementation_difficulty": "medium",
            "priority": "medium_high",
            "priority_score": 7,
        },
        {
            "modality": "mutation",
            "current_state": "gene-level event matrix",
            "opportunity": "TMB/burden + driver-centric + pathway alteration score",
            "expected_benefit": "提升临床分层与分类稳健性",
            "implementation_difficulty": "low_medium",
            "priority": "high",
            "priority_score": 8,
        },
        {
            "modality": "cnv",
            "current_state": "gene-level gistic",
            "opportunity": "focal event count + arm/pathway-level CNV summary",
            "expected_benefit": "增强结构性特征表达",
            "implementation_difficulty": "medium",
            "priority": "medium_high",
            "priority_score": 7,
        },
        {
            "modality": "rppa",
            "current_state": "available but excluded from final mainline",
            "opportunity": "late-fusion auxiliary protein evidence",
            "expected_benefit": "增强证据层而不强制缩样本",
            "implementation_difficulty": "low_medium",
            "priority": "medium",
            "priority_score": 6,
        },
    ]
    return pd.DataFrame(rows)


def _build_clinical_field_opportunities(root: Path) -> pd.DataFrame:
    clinical_path = root / "data" / "raw" / "clinical" / "TCGA.STAD.sampleMap_STAD_clinicalMatrix"
    survival_path = root / "data" / "raw" / "clinical" / "survival_STAD_survival.txt"

    clinical_cols = _read_tsv_header(clinical_path)
    survival_cols = _read_tsv_header(survival_path)

    used_clinical = {
        "age_at_initial_pathologic_diagnosis",
        "gender",
        "pathologic_stage",
        "histological_type",
    }
    used_survival = {"OS", "OS.time", "PFI", "PFI.time"}

    candidate_fields = [
        "DSS", "DSS.time", "DFI", "DFI.time",
        "neoplasm_histologic_grade",
        "pathologic_T", "pathologic_N", "pathologic_M",
        "vital_status", "days_to_death", "days_to_last_followup",
        "new_tumor_event_after_initial_treatment",
        "primary_therapy_outcome_success",
        "radiation_therapy", "additional_pharmaceutical_therapy",
        "additional_radiation_therapy",
        "h_pylori_infection",
        "family_history_of_stomach_cancer",
        "sample_type", "sample_type_id",
    ]

    rows: list[dict[str, Any]] = []

    for c in sorted(set(clinical_cols) & set(candidate_fields + list(used_clinical))):
        rows.append(
            {
                "field_name": c,
                "source": "clinical_matrix",
                "currently_used": c in used_clinical,
                "usage_status": "used" if c in used_clinical else "underused",
                "potential_analysis": (
                    "cluster-clinical association in current pipeline"
                    if c in used_clinical
                    else "stratification / covariate adjustment / subgroup interpretation"
                ),
                "why_not_fully_used": (
                    "already in use" if c in used_clinical else "尚未纳入统一临床变量建模与多变量验证"
                ),
            }
        )

    for c in [x for x in ["OS", "OS.time", "PFI", "PFI.time", "DSS", "DSS.time", "DFI", "DFI.time"] if x in survival_cols]:
        rows.append(
            {
                "field_name": c,
                "source": "survival_table",
                "currently_used": c in used_survival,
                "usage_status": "used" if c in used_survival else "underused",
                "potential_analysis": (
                    "final survival endpoint"
                    if c in used_survival
                    else "extra endpoint consistency validation"
                ),
                "why_not_fully_used": (
                    "already in use" if c in used_survival else "终点扩展尚未纳入正式报告主线"
                ),
            }
        )

    return pd.DataFrame(rows)


def _build_top3() -> pd.DataFrame:
    rows = [
        {
            "rank": 1,
            "recommendation": "推进 at-least-4 partial-fusion 正式扩容并统一与 complete-case 的评估框架",
            "target_layer": "sample",
            "expected_benefit": "以较低异质性代价增加样本，提升稳定性与统计功效",
            "expected_risk": "模态缺失引入偏差，需要缺失机制控制",
            "implementation_cost": "medium",
            "why_high_priority": "可直接挽回约 27 个样本，且已有 partial 产物基础",
        },
        {
            "rank": 2,
            "recommendation": "methylation 从 probe-level 升级到 promoter/gene-level 注释聚合特征",
            "target_layer": "feature",
            "expected_benefit": "降低噪声并提升解释性，改善聚类与分类稳健性",
            "expected_risk": "聚合策略不当可能损失局部信号",
            "implementation_cost": "medium",
            "why_high_priority": "现有 probe map 与过滤日志齐全，具备直接落地条件",
        },
        {
            "rank": 3,
            "recommendation": "将 DSS/DFI + 关键临床变量纳入多终点一致性与分层证据链",
            "target_layer": "evidence",
            "expected_benefit": "显著增强结论完整性与说服力",
            "expected_risk": "终点缺失值与事件数可能限制显著性",
            "implementation_cost": "low_medium",
            "why_high_priority": "不改模型也能增强证据强度，最适合下一轮先做",
        },
    ]
    return pd.DataFrame(rows)


def _make_figures(root: Path, assets: pd.DataFrame, salvage: pd.DataFrame, feat: pd.DataFrame, clinical: pd.DataFrame) -> None:
    fig_dir = root / "results" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # 1) data asset map
    level_counts = assets.groupby("asset_level").size().reindex(["raw", "annotation", "interim", "results"]).fillna(0)
    plt.figure(figsize=(10, 4.8))
    x = np.arange(len(level_counts))
    plt.bar(x, level_counts.values, color=["#457b9d", "#2a9d8f", "#e9c46a", "#e76f51"])
    plt.xticks(x, ["raw", "annotation", "interim", "results"])
    plt.ylabel("file count")
    plt.title("Data Asset Map: raw -> annotation -> interim -> results")
    for i, v in enumerate(level_counts.values):
        plt.text(i, v + 0.5, str(int(v)), ha="center", fontsize=9)
    plt.tight_layout()
    plt.savefig(fig_dir / "data_asset_map.png", dpi=240)
    plt.close()

    # 2) sample salvage opportunity chart
    ss = salvage[salvage["opportunity"].isin([
        "strict_main5_complete_case",
        "salvage_to_partial_atleast4_main5",
        "salvage_to_partial_atleast3_main5",
    ])].copy()
    labels = [
        "strict_main5",
        "recover_to_ge4",
        "recover_to_ge3",
    ]
    vals = [
        int(ss.loc[ss["opportunity"] == "strict_main5_complete_case", "current_n"].iloc[0]),
        int(ss.loc[ss["opportunity"] == "salvage_to_partial_atleast4_main5", "recoverable_n"].iloc[0]),
        int(ss.loc[ss["opportunity"] == "salvage_to_partial_atleast3_main5", "recoverable_n"].iloc[0]),
    ]
    plt.figure(figsize=(8, 4.8))
    plt.bar(labels, vals, color=["#264653", "#2a9d8f", "#f4a261"])
    plt.title("Sample Salvage Opportunity: complete-case vs recoverable space")
    plt.ylabel("sample count")
    for i, v in enumerate(vals):
        plt.text(i, v + 1, str(v), ha="center", fontsize=9)
    plt.tight_layout()
    plt.savefig(fig_dir / "sample_salvage_opportunity_chart.png", dpi=240)
    plt.close()

    # 3) feature opportunity priority chart
    f = feat.sort_values("priority_score", ascending=True)
    plt.figure(figsize=(9, 5.2))
    plt.barh(f["modality"], f["priority_score"], color="#3a86ff")
    plt.title("Feature Opportunity Priority by Modality")
    plt.xlabel("priority score (higher = more attractive)")
    for y, xval in zip(f["modality"], f["priority_score"]):
        plt.text(xval + 0.05, y, f"{xval}", va="center", fontsize=9)
    plt.xlim(0, 10)
    plt.tight_layout()
    plt.savefig(fig_dir / "feature_opportunity_priority_chart.png", dpi=240)
    plt.close()

    # 4) clinical field usage map
    cm = clinical.groupby(["source", "usage_status"]).size().reset_index(name="count")
    pivot = cm.pivot(index="source", columns="usage_status", values="count").fillna(0)
    for c in ["used", "underused"]:
        if c not in pivot.columns:
            pivot[c] = 0
    pivot = pivot[["used", "underused"]]
    pivot.plot(kind="bar", stacked=True, figsize=(8, 4.8), color=["#2a9d8f", "#e76f51"])
    plt.title("Clinical Field Usage Map: used vs underused")
    plt.ylabel("field count")
    plt.xticks(rotation=0)
    plt.legend(title="status")
    plt.tight_layout()
    plt.savefig(fig_dir / "clinical_field_usage_map.png", dpi=240)
    plt.close()


def _write_docs(root: Path, assets: pd.DataFrame, underused: pd.DataFrame, salvage: pd.DataFrame, top3: pd.DataFrame) -> None:
    docs = root / "docs"
    docs.mkdir(parents=True, exist_ok=True)

    raw_assets = assets[assets["asset_level"] == "raw"]
    raw_count = int(raw_assets.shape[0])
    teacher_pkg_raw_count = int(raw_assets[raw_assets["file_path"].str.contains("experiment_data", case=False, regex=False)].shape[0])
    copied_raw_count = int(raw_assets[raw_assets["file_path"].str.contains("data\\raw", case=False, regex=False)].shape[0])
    raw_used_final = int(raw_assets["currently_used_in_final_model"].sum())
    underused_count = int(underused.shape[0])
    complete5 = int(salvage.loc[salvage["opportunity"] == "strict_main5_complete_case", "current_n"].iloc[0])
    ge4_rec = int(salvage.loc[salvage["opportunity"] == "salvage_to_partial_atleast4_main5", "recoverable_n"].iloc[0])
    ge3_rec = int(salvage.loc[salvage["opportunity"] == "salvage_to_partial_atleast3_main5", "recoverable_n"].iloc[0])

    atlas = f"""# 35 Dataset Atlas（数据资产深度盘点）

## 1. 盘点结论概览
- 原始层（raw）共识别数据文件 {raw_count} 个（含 data/raw 与 experiment_data 的可读数据文件）。
- 其中老师原始包目录（experiment_data）可读数据文件：{teacher_pkg_raw_count} 个。
- 其中规范化复制目录（data/raw）可读数据文件：{copied_raw_count} 个。
- 原始层中已进入 final 主模型直接链路的文件数：{raw_used_final}。
- 识别到未充分利用或可进一步利用资产：{underused_count} 项。
- 当前 strict main5 complete-case 样本数：{complete5}。
- 可通过 partial-fusion at-least-4 潜在挽回样本：{ge4_rec}。
- 可通过 partial-fusion at-least-3 潜在挽回样本：{ge3_rec}。

## 2. 原始数据全景
核心组学与临床资源：
- mutation / cnv / methylation / rna / mirna / rppa
- clinical matrix / survival（含 OS/PFI/DSS/DFI）
- methylation probe map 注释资源

对应全清单见：
- [results/tables/data_asset_inventory.csv](results/tables/data_asset_inventory.csv)

## 3. 当前项目实际使用链路
raw -> interim -> final 的链路图：
![data asset map](../results/figures/data_asset_map.png)

关键说明：
- final 主线主要依赖 main5（mutation/cnv/methylation/rna/mirna）interim 矩阵。
- rppa 原始数据存在且可用，但因样本交集代价在 final 主线中未纳入。
- 临床与生存数据已用于 OS/PFI 相关验证，但未完全展开到 DSS/DFI 与更多临床字段。

## 4. 未充分利用资产画像
详表见：
- [results/tables/unused_or_underused_assets.csv](results/tables/unused_or_underused_assets.csv)

高价值未充分利用方向：
1. methylation probe 注释到 gene/promoter/region 的聚合特征
2. survival 的 DSS/DFI 终点
3. clinical matrix 的治疗/复发/分级字段
4. partial-fusion 扩容样本空间
5. rppa 作为辅助证据层

## 5. 样本挽回空间
![sample salvage](../results/figures/sample_salvage_opportunity_chart.png)

对应表：
- [results/tables/sample_salvage_opportunities.csv](results/tables/sample_salvage_opportunities.csv)

解释：
- 在不改动主线结论前提下，at-least-4 是最稳妥扩容入口。
- at-least-3 提供更大样本空间，但异质性与缺失机制风险更高。

## 6. 特征增强潜力
![feature opportunity](../results/figures/feature_opportunity_priority_chart.png)

对应表：
- [results/tables/feature_enrichment_opportunities.csv](results/tables/feature_enrichment_opportunities.csv)

解释：
- methylation 聚合特征、RNA/miRNA pathway 表达、mutation/CNV 汇总特征是下一轮重点。

## 7. 临床字段使用缺口
![clinical usage](../results/figures/clinical_field_usage_map.png)

对应表：
- [results/tables/clinical_field_opportunities.csv](results/tables/clinical_field_opportunities.csv)

解释：
- 当前已用字段相对有限，仍有大量可用临床变量可增强证据链完整性。

## 8. 最值得挖掘的资源
1. methylation 注释驱动聚合（噪声控制 + 解释增强）
2. partial-fusion at-least-4 样本扩容（提升统计功效）
3. DSS/DFI + 治疗相关临床字段（增强证据充分性）
"""
    (docs / "35_dataset_atlas.md").write_text(atlas, encoding="utf-8")

    gap = """# 36 Data Utilization Gap Analysis

## 1. 当前主线已用到了什么
- main5 五组学 complete-case 矩阵（interim）
- baseline/proposed/improvement 结果链路
- OS/PFI 生存分离与部分临床变量关联
- subtype classifier 内部验证

## 2. 还缺什么
- 多终点（DSS/DFI）一致性证据
- 更系统的临床字段整合与分层分析
- 注释驱动的聚合特征（尤其 methylation）
- partial-fusion 正式扩容主线

## 3. 缺口类型判断
- 数据本身没有：外部独立队列（当前目录未见外部验证数据）
- 尚未充分利用：
  1. probe map 到聚合特征
  2. DSS/DFI 与治疗相关临床字段
  3. RPPA 辅助证据层
  4. weakly paired 样本扩容

## 4. 哪些缺口最值得补
优先级高：
1. at-least-4 partial-fusion 扩容
2. methylation 注释聚合
3. 临床与生存终点扩展（DSS/DFI + 关键临床字段）

对应资产依据：
- [results/tables/unused_or_underused_assets.csv](results/tables/unused_or_underused_assets.csv)
- [results/tables/sample_salvage_opportunities.csv](results/tables/sample_salvage_opportunities.csv)
- [results/tables/clinical_field_opportunities.csv](results/tables/clinical_field_opportunities.csv)
"""
    (docs / "36_data_utilization_gap_analysis.md").write_text(gap, encoding="utf-8")

    rec_lines = [
        "# 37 Recommendations for Next Optimization Round",
        "",
        "## 1. 建议清单（3-5件事）",
        "详表：",
        "- [results/tables/top3_next_steps_recommended.csv](results/tables/top3_next_steps_recommended.csv)",
        "",
    ]
    for _, r in top3.sort_values("rank").iterrows():
        rec_lines.extend(
            [
                f"### {int(r['rank'])}. {r['recommendation']}",
                f"- 目标层：{r['target_layer']}",
                f"- 预期收益：{r['expected_benefit']}",
                f"- 风险：{r['expected_risk']}",
                f"- 成本：{r['implementation_cost']}",
                "",
                "更可能改善指标：",
                (
                    "- stability / clinical relevance / evidence strength"
                    if int(r["rank"]) == 1
                    else "- cluster quality / classification performance / interpretability"
                    if int(r["rank"]) == 2
                    else "- clinical relevance / evidence strength"
                ),
                "",
            ]
        )

    rec_lines.extend(
        [
            "## 2. 结论",
            "下一轮最值得优先做的是：先做 at-least-4 partial-fusion 扩容与统一评估框架，同时并行推进 methylation 注释聚合与 DSS/DFI 终点补强。",
        ]
    )
    (docs / "37_recommendations_for_next_optimization_round.md").write_text("\n".join(rec_lines), encoding="utf-8")


def main() -> None:
    root = _project_root()
    out_tables = root / "results" / "tables"
    out_tables.mkdir(parents=True, exist_ok=True)

    assets = _collect_assets(root)
    assets.to_csv(out_tables / "data_asset_inventory.csv", index=False, encoding="utf-8")

    underused = assets[(assets["partially_used"] == True) | (assets["unused_but_potentially_useful"] == True)].copy()
    underused.to_csv(out_tables / "unused_or_underused_assets.csv", index=False, encoding="utf-8")

    salvage = _build_sample_salvage(root)
    salvage.to_csv(out_tables / "sample_salvage_opportunities.csv", index=False, encoding="utf-8")

    feat = _build_feature_opportunities()
    feat.to_csv(out_tables / "feature_enrichment_opportunities.csv", index=False, encoding="utf-8")

    clinical = _build_clinical_field_opportunities(root)
    clinical.to_csv(out_tables / "clinical_field_opportunities.csv", index=False, encoding="utf-8")

    top3 = _build_top3()
    top3.to_csv(out_tables / "top3_next_steps_recommended.csv", index=False, encoding="utf-8")

    _make_figures(root, assets, salvage, feat, clinical)
    _write_docs(root, assets, underused, salvage, top3)

    print("Dataset atlas generated:")
    print("- results/tables/data_asset_inventory.csv")
    print("- results/tables/unused_or_underused_assets.csv")
    print("- results/tables/sample_salvage_opportunities.csv")
    print("- results/tables/feature_enrichment_opportunities.csv")
    print("- results/tables/clinical_field_opportunities.csv")
    print("- results/tables/top3_next_steps_recommended.csv")
    print("- docs/35_dataset_atlas.md")
    print("- docs/36_data_utilization_gap_analysis.md")
    print("- docs/37_recommendations_for_next_optimization_round.md")
    print("- results/figures/data_asset_map.png")
    print("- results/figures/sample_salvage_opportunity_chart.png")
    print("- results/figures/feature_opportunity_priority_chart.png")
    print("- results/figures/clinical_field_usage_map.png")


if __name__ == "__main__":
    main()
