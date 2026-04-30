# 13 建模就绪性与 Baseline（Round1 + Main5 Final）

## 1. 建模就绪性审计结论（Phase 11）
正式审计表：
- [results/tables/modeling_readiness_summary.csv](results/tables/modeling_readiness_summary.csv)

核心结论：
- 可直接用于 round1 baseline 的组学：mutation、cnv、rna、mirna。
- methylation 当前状态：仅有 chunked summary + preview（无完整 sample-by-feature 建模矩阵），暂不纳入 baseline round1。
- rppa 当前状态：有可建模矩阵，但根据既有决策定位为补充轨道，不进入主 baseline round1。

## 2. round1 实际采用组学组合
- 实际运行组合：main4 = mutation + cnv + rna + mirna。
- 原因：
  - 主5中 methylation 尚未完成可建模矩阵落地；
  - 为保证本轮 baseline 真实可运行与可复核，先采用稳定可用 4 组学。
- 下一轮目标：补齐 methylation 可建模矩阵后，升级为 main5 对照运行。

## 3. 相似网络构建与输出
已完成：
- 代码：[src/feature_engineering/build_similarity_matrices.py](src/feature_engineering/build_similarity_matrices.py)
- 对齐矩阵输出目录：data/interim/aligned_matrices
- 相似矩阵输出目录：[results/tables/similarity_matrices](results/tables/similarity_matrices)
- 审计输出：
  - [results/tables/modality_modeling_shapes.csv](results/tables/modality_modeling_shapes.csv)
  - [results/tables/similarity_build_log.csv](results/tables/similarity_build_log.csv)

已生成图：
- [results/figures/sample_similarity_heatmap_mutation.png](results/figures/sample_similarity_heatmap_mutation.png)
- [results/figures/sample_similarity_heatmap_cnv.png](results/figures/sample_similarity_heatmap_cnv.png)
- [results/figures/sample_similarity_heatmap_rna.png](results/figures/sample_similarity_heatmap_rna.png)
- [results/figures/sample_similarity_heatmap_mirna.png](results/figures/sample_similarity_heatmap_mirna.png)

## 4. Baseline round1 实际运行结果（Phase 12）
结果表：
- [results/tables/baseline_round1_metrics.csv](results/tables/baseline_round1_metrics.csv)
- [results/tables/baseline_round1_cluster_sizes.csv](results/tables/baseline_round1_cluster_sizes.csv)
- [results/tables/fused_similarity_summary.csv](results/tables/fused_similarity_summary.csv)

标签文件（按方法与 K 分开落盘）：
- [results/tables](results/tables) 下 `baseline_round1_labels_*` 系列文件。

主要观察：
- early_fusion_kmeans 在 k=2 时稳定性较高（consensus_stability≈0.894），cluster size 也相对平衡（191 vs 173）。
- equal_weight_fusion 虽在 silhouette 上数值较高（k=2 时≈0.126），但簇规模极端失衡（363 vs 1），提示当前实现下存在退化分簇风险。
- 因此本轮更建议将 `early_fusion_kmeans, k=2` 作为首轮候选基线。

## 5. 稳定性评估（Phase 13 部分）
稳定性结果表：
- [results/tables/cluster_stability_round1.csv](results/tables/cluster_stability_round1.csv)

已真实计算指标：
- repeated subsampling（6 次）
- mean NMI
- mean ARI
- consensus_stability（NMI 与 ARI 均值）

说明：
- 当前为首轮轻量稳定性版本，已真实运行；
- 下一轮计划增加更重的 consensus clustering 与 bootstrap 版本。

## 6. 临床关联首轮（Phase 13 部分）
临床关联表：
- [results/tables/clinical_association_round1.csv](results/tables/clinical_association_round1.csv)

图：
- KM（OS/PFI）：`km_<method>_<modalityset>_bestk_<endpoint>.png`
- [results/figures/clinical_association_heatmap_round1.png](results/figures/clinical_association_heatmap_round1.png)

结果解释边界：
- 当前结论属于首轮发现性分析（discovery in internal-validation context），不是外部独立验证。

## 7. 仍待补强项
1. methylation 完整 sample-by-feature 可建模矩阵（当前仅 summary/preview）。
2. equal-weight fusion 的退化簇问题（需改进 affinity 构造或引入图正则/密度约束）。
3. 稳定性评估从 6 次扩展到更高重复与更完整共识矩阵流程。
4. 临床变量清洗（尤其 stage/Lauren 编码）继续标准化。

---

## 8. Phase 16-17 更新：main5 正式就绪与 baseline

### 8.1 methylation 就绪状态更新
- 新增可建模矩阵：`data/interim/methylation_round2_modeling.csv`
- 审计表：
  - [results/tables/methylation_processing_summary.csv](results/tables/methylation_processing_summary.csv)
  - [results/tables/modeling_readiness_main5.csv](results/tables/modeling_readiness_main5.csv)
- 结论：`final_matrix_ready=True`，main5 可正式运行。

### 8.2 main5 baseline 结果（独立于 round1 main4）
- 结果文件：
  - [results/tables/baseline_main5_metrics.csv](results/tables/baseline_main5_metrics.csv)
  - [results/tables/baseline_main5_cluster_sizes.csv](results/tables/baseline_main5_cluster_sizes.csv)
  - [results/tables/fused_similarity_summary_main5.csv](results/tables/fused_similarity_summary_main5.csv)
- 关键观察：
  - `early_fusion_kmeans@k=2` 在 main5 下 silhouette≈0.053，簇平衡优于等权融合。
  - `equal_weight_fusion` 在 main5 下多 K 出现近似退化（silhouette<0，cluster_balance 接近 0）。
  - 说明引入 methylation 后并未自动提升所有 baseline，需谨慎解释“更多组学=更好”这一假设。

### 8.3 warning 记录
- 在 main5 运行中仍出现 spectral graph connectivity warning（graph not fully connected）。
- 已按风险项记录，相关结论以“内部发现性”表达，不作为临床部署级证据。

## 9. 结论边界（必须区分）
1. round1 main4：用于建立可运行基线与方法对比起点。
2. main5 final：在补强 methylation 后得到的正式阶段结果。
3. 两者不可直接混写为同一轮结果，报告中必须显式标注 `main4/round1` 与 `main5/final`。
