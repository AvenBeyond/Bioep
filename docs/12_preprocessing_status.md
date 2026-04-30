# 12 预处理落地状态

## 1. 状态说明
- 已真实运行：代码已执行并产出文件与日志。
- 已实现未执行：代码已写好但尚未在当前项目数据上运行。
- 仅接口预留：函数签名或逻辑框架存在，但未形成有效处理链。

## 2. 各组学当前进展

| modality | 当前状态 | 输入 | 输出 | 备注 |
|---|---|---|---|---|
| mutation | 已真实运行 | data/raw/mutation/STAD_mc3_gene_level.txt | data/interim/mutation_round1.csv | 非零事件二值化 + 低频过滤 |
| cnv | 已真实运行 | data/raw/cnv/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes | data/interim/cnv_round1.csv | -2..2 裁剪离散化 + 方差过滤 |
| rna | 已真实运行 | data/raw/rna/HiSeqV2 | data/interim/rna_round1.csv | 低表达过滤 + top variance + z-score |
| mirna | 已真实运行 | data/raw/mirna/miRNA_HiSeq_gene | data/interim/mirna_round1.csv | NA处理 + 低表达过滤 + top variance + z-score |
| methylation | 已真实运行（round2 可建模版） | data/raw/methylation/HumanMethylation450 | data/interim/methylation_round2_modeling.csv | 已完成缺失过滤+高方差筛选+probe-level建模矩阵；gene/promoter 聚合接口保留 |
| rppa | 已真实运行 | data/raw/rppa/RPPA | data/interim/rppa_round1.csv | 中位数填补 + 非零率过滤 + z-score；补充组学轨道 |

## 3. 维度变更日志
- 日志文件：[results/logs/preprocessing_dimension_changes.csv](results/logs/preprocessing_dimension_changes.csv)
- 已记录字段：
  - modality
  - input_shape
  - output_shape
  - sample_count_before
  - sample_count_after
  - feature_count_before
  - feature_count_after
  - filtering_steps
  - read_mode
  - notes

## 4. 当前真实结果摘要
- mutation：439x40543 -> 364x4598
- cnv：441x24776 -> 364x24776
- rna：450x20530 -> 364x5000
- mirna：428x2178 -> 364x1000
- methylation round1：398x485577 -> 364x2000（preview 规模）
- methylation round2：485577 probes -> 5000 features（sample=364，`final_matrix_ready=True`）
- rppa：357x227 -> 357x227

## 5. 风险与未完成项
- methylation 已可进入 main5 正式实验，但当前为 probe-level 高方差特征，不等价于正式 gene/promoter 生物学聚合。
- mirna 缺失较高，后续需比较多种缺失处理策略对聚类稳定性的影响。
- patient-level 多样本策略尚未进入最终建模聚合逻辑，当前仍按 sample-level 主键运行。

## 6. Phase 16 关键产物
- [data/interim/methylation_round2_modeling.csv](data/interim/methylation_round2_modeling.csv)
- [results/tables/methylation_processing_summary.csv](results/tables/methylation_processing_summary.csv)
- [results/tables/methylation_feature_filter_log.csv](results/tables/methylation_feature_filter_log.csv)
- [results/figures/methylation_missingness_summary.png](results/figures/methylation_missingness_summary.png)
- [results/figures/methylation_variance_distribution.png](results/figures/methylation_variance_distribution.png)

## 7. 下一步计划
1. 在现有 probe-level matrix 基础上探索 gene/promoter 聚合对稳定性与解释性的影响。
2. 结合消融结果判断 methylation 在最终主模型中的净贡献。
3. 保留 round1 与 final 结果的清晰边界，不混用结论。
