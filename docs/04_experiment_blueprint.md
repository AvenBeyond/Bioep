# 04 实验蓝图

## 1. 数据读入
- 统一由 [main.py](main.py) + [src/data_loading/load_all_modalities.py](src/data_loading/load_all_modalities.py) 读取。
- 已完成 Phase 7 真实样本解析与交集统计，正式产物见 [results/tables](results/tables)。

## 2. 样本 ID 统一策略
- 原始样本 ID 多为 `TCGA-XX-XXXX-01` 风格。
- 正式规则：
   - normalized_sample_id = 前 16 位（例：TCGA-3M-AB46-01）
   - normalized_patient_id = 前 12 位（例：TCGA-3M-AB46）
- patient-level 统计用于报告与重复样本识别；样本级建模主键仍为 normalized_sample_id。

## 3. 各组学预处理设计
- mutation：非零事件二值化 + 低频特征过滤（min_event_rate）。
- cnv：裁剪到 -2..2 并离散化 + 方差过滤。
- rna：低表达过滤 + 高变异筛选（top_var）+ z-score。
- mirna：NA 转换与填补 + 低表达过滤 + 高变异筛选 + z-score。
- methylation：分块扫描 + 缺失统计 + top-variance 预览 + probe map 匹配覆盖评估。
- rppa：中位数填补 + 低非零率过滤 + z-score，定位为补充组学。

已执行入口：[src/preprocessing/run_round1_preprocessing.py](src/preprocessing/run_round1_preprocessing.py)。

## 4. 特征筛选/降维
- 首轮优先采用过滤法（方差/非零率/缺失率）。
- 后续可加 PCA/NMF/稀疏表示用于降维对比。

当前已落地：
- mutation 输出 4598 特征。
- cnv 输出 24776 特征。
- rna 输出 5000 特征。
- mirna 输出 1000 特征。
- methylation 输出 summary + preview（2000 高变异探针预览）。

## 5. 相似网络构建
- 在 [src/feature_engineering/build_similarity_matrices.py](src/feature_engineering/build_similarity_matrices.py) 预留三类接口：
  - Pearson similarity
  - Cosine affinity
  - Euclidean-based Gaussian affinity

## 6. baseline 设计
1. Early fusion clustering
   - 各组学特征拼接后聚类
2. Equal-weight network fusion
   - 各组学相似网络等权融合

## 7. proposed 设计
在 [src/clustering/proposed_weighted_fusion.py](src/clustering/proposed_weighted_fusion.py) 预留：
- 权重计算接口（$w_m$）
- partial fusion 逻辑
- K 搜索与综合评分接口

## 8. 评估指标
- 聚类质量：Silhouette/CH/DBI
- 稳定性：NMI/ARI/Jaccard（重采样）
- 临床相关性：KM/log-rank/Cox
- 分类器内部验证：Accuracy/Macro-F1/Weighted-F1/AUROC/混淆矩阵

## 9. 消融实验
- equal-weight vs weighted-fusion
- weighted-fusion 去除临床项
- 去除 partial fusion
- 去除部分组学
- 是否纳入 RPPA

## 10. 关键图表输出计划
- 样本交集图（UpSet/heatmap）
- K 搜索综合评分图
- 生存曲线图（OS/PFI）
- 亚型临床关联图
- 消融结果图

## 11. 主分析队列定义草案（Phase 8）
- 主分析组学：mutation + cnv + methylation + rna + mirna。
- 主分析样本集：5 组学交集 + Primary Tumor（当前 364，且与 all sample 交集一致）。
- RPPA 定位：补充组学，不纳入主分析交集定义。

## 12. partial fusion 触发条件（草案）
- 当样本未达到 5 组学完整覆盖时，若覆盖组学数 >= 3，则允许进入 partial fusion 分支。
- 当前覆盖统计（主 5 组学）：>=3 为 436，>=4 为 391，=5 为 364。
- 解释策略：主结论以完整 5 组学为主，partial fusion 用于鲁棒性补充与样本利用率分析。

## 13. 预处理实际落地顺序
1. 先跑样本交集与对齐统计。
2. 再跑 mutation/cnv/rna/mirna 的可执行预处理。
3. methylation 采用 chunked_read 先产出 summary 与 preview。
4. RPPA 单独生成补充轨道预处理结果。
5. 统一记录到 [results/logs/preprocessing_dimension_changes.csv](results/logs/preprocessing_dimension_changes.csv)。

## 14. Baseline Round1 实际执行路径（Phase 11-12）
1. 建模就绪性审计：
   - [results/tables/modeling_readiness_summary.csv](results/tables/modeling_readiness_summary.csv)
2. 本轮实际组学组合：main4（mutation+cnv+rna+mirna）
3. K 搜索：2..6
4. baseline 方法：
   - early_fusion_kmeans
   - early_fusion_spectral
   - equal_weight_fusion
5. 输出：
   - [results/tables/baseline_round1_metrics.csv](results/tables/baseline_round1_metrics.csv)
   - [results/tables/baseline_round1_cluster_sizes.csv](results/tables/baseline_round1_cluster_sizes.csv)
   - [results/tables/fused_similarity_summary.csv](results/tables/fused_similarity_summary.csv)

## 15. 首轮临床验证策略（Phase 13）
- 采用“候选最佳 k（按 silhouette）”进行 KM（OS/PFI）与临床变量关联。
- 结果定位为 discovery in internal-validation context。
- 输出：
  - [results/tables/clinical_association_round1.csv](results/tables/clinical_association_round1.csv)
  - [results/figures/clinical_association_heatmap_round1.png](results/figures/clinical_association_heatmap_round1.png)
