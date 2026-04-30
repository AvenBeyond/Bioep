# 01 数据总览

统计依据说明：
- 文件完整性与总行列依据：[results/logs/raw_file_inventory.csv](results/logs/raw_file_inventory.csv)、[results/logs/raw_file_inventory.json](results/logs/raw_file_inventory.json)、[docs/11_data_reading_integrity_check.md](docs/11_data_reading_integrity_check.md)。
- 样本方向、样本数、交集、主分析范围依据：
  [results/tables/modality_dimension_summary.csv](results/tables/modality_dimension_summary.csv)、
  [results/tables/sample_inventory.csv](results/tables/sample_inventory.csv)、
  [results/tables/sample_overlap_matrix.csv](results/tables/sample_overlap_matrix.csv)、
  [results/tables/multiomics_intersection_summary.csv](results/tables/multiomics_intersection_summary.csv)。

## 1. 组学维度与方向（Phase 7 实测）

| modality | raw_file_name | sample_axis_orientation | raw_sample_count | raw_feature_count | value_type | missingness_summary |
|---|---|---|---:|---:|---|---|
| mutation | STAD_mc3_gene_level.txt | column | 439 | 40543 | binary_or_event_indicator | estimated_missing=0.000%_from_first_20000_rows |
| cnv | Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes | column | 441 | 24776 | discrete_gistic_-2_to_2 | estimated_missing=0.000%_from_first_20000_rows |
| methylation | HumanMethylation450 | column | 398 | 485577 | continuous_beta_0_to_1 | estimated_missing=18.789%_from_first_20000_rows |
| rna | HiSeqV2 | column | 450 | 20530 | continuous_expression | estimated_missing=0.000%_from_first_20000_rows |
| mirna | miRNA_HiSeq_gene | column | 428 | 2178 | continuous_expression_with_na | estimated_missing=63.723%_from_first_2178_rows |
| rppa | RPPA | column | 357 | 227 | continuous_protein_expression | estimated_missing=3.850%_from_first_227_rows |

说明：样本均位于列（第一行为样本 ID，第一列为特征 ID）。

## 2. 样本 ID 规范化结果

- normalized_sample_id：TCGA barcode 前 16 位（例：TCGA-3M-AB46-01）。
- normalized_patient_id：TCGA barcode 前 12 位（例：TCGA-3M-AB46）。
- sample_type_candidate：由 barcode 第 14-15 位解析，当前主要为 01（Primary Tumor）与 11（Solid Tissue Normal）。

关键统计：
- 全部样本记录行数：2513。
- sample type 分布：Primary Tumor 2435，Solid Tissue Normal 78。
- 各组学 Primary Tumor 样本：mutation 439、cnv 441、methylation 396、rna 415、mirna 387、rppa 357。
- 同一 patient 出现多样本的 patient 数：rna 32、mirna 41、methylation 2（其余为 0）。

## 3. 样本交集概况（真实统计）

两两交集矩阵见：[results/tables/sample_overlap_matrix.csv](results/tables/sample_overlap_matrix.csv)。

关键交集：
- 主分析 5 组学（mutation+cnv+methylation+rna+mirna）交集：364。
- 6 组学（再加 rppa）交集：288。
- 纳入 rppa 后主队列交集减少：76（约 20.9%）。
- patient-level 交集：主 5 组学 367，含 rppa 的 6 组学 291。

partial fusion 覆盖摘要（主 5 组学）：
- 覆盖至少 3 组学：436。
- 覆盖至少 4 组学：391。
- 覆盖 5 组学：364。

## 4. 主分析范围建议（数据层）

- 主分析默认队列：5 组学交集的 Primary Tumor 样本（当前为 364）。
- RPPA：当前不纳入主分析队列，保留为补充/消融分析。
- 理由：纳入 RPPA 会显著压缩交集样本数，不利于主分析稳定性。

## 5. placeholder 替换状态

- 已替换：
  [results/tables/sample_inventory_placeholder.csv](results/tables/sample_inventory_placeholder.csv) 不再作为正式引用。
- 正式文件：
  [results/tables/sample_inventory.csv](results/tables/sample_inventory.csv)。

## 6. 风险与下一步

- methylation 维度极高，完整探针级聚合不宜一步到位，当前已用分块扫描与预览筛选策略。
- mirna 缺失较高，后续需要评估不同缺失处理方案对聚类稳定性的影响。
- 后续所有预处理变更以 [results/logs/preprocessing_dimension_changes.csv](results/logs/preprocessing_dimension_changes.csv) 追踪。
