# 05 验证与测试协议

## 1. 现实前提：无独立外部测试集
当前项目数据均来自 TCGA-STAD，同源数据用于亚型发现与内部验证。文档与汇报中必须避免将内部交叉验证表述为“外部测试”。

## 1.1 当前样本集定义（Phase 8）
- 主分析样本集：5 组学（mutation/cnv/methylation/rna/mirna）交集样本，当前 n=364。
- 含 RPPA 的 6 组学交集样本：n=288。
- 说明：RPPA 纳入会导致交集减少 76（约 20.9%），故当前将 RPPA 作为补充样本集与补充分析轨道。

## 2. 结论类型边界
- 发现性分析（discovery）：
  - 基于完整可用队列进行亚型发现、生物学解释、临床关联探索。
- 泛化评估（internal validation）：
  - 通过内部重采样与交叉验证评估方法稳定性与泛化趋势。

## 3. 无监督内部验证
1. K 搜索范围：2..6
2. repeated subsampling + consensus clustering
3. bootstrap / repeated resampling
4. 稳定性指标：NMI、ARI、Jaccard
5. 聚类质量指标：Silhouette、CH、DBI
6. leave-one-modality-out robustness
7. perturbation robustness（轻噪声/随机删样本）

执行边界：
- 主结果优先在主分析样本集（n=364）上报告。
- 若报告 partial fusion 结果，必须明确其样本覆盖阈值与样本规模（如 >=3 组学覆盖 n=436）。

## 4. 临床相关性验证
1. Kaplan-Meier survival analysis
2. log-rank test
3. Cox regression
4. subtype 与 stage/grade/Lauren 等关联
5. 主终点 OS/PFI，DSS/DFI 作为补充（视完整性）

## 5. subtype predictor 内部验证
1. 以最终 subtype label 构建分类任务
2. 候选模型：RF/SVM/XGBoost
3. nested CV 或 repeated stratified k-fold CV
4. 所有预处理严格在训练折内完成，防止泄漏
5. 输出 Accuracy/Macro-F1/Weighted-F1/AUROC/Confusion Matrix

执行边界：
- internal validation 的训练/验证拆分以主分析样本集为主。
- 若在补充样本集（含 RPPA 或 partial fusion）上运行，需与主分析结果分开报告，不可混写。

## 6. 如何降低过拟合和乐观偏差
- 严格区分 discovery 与 validation
- 重采样和重复交叉验证
- 消融与鲁棒性测试
- 报告指标时给出均值与方差/置信区间

## 7. 局限性声明
- 缺乏外部独立验证队列
- 弱配对和缺失可能影响稳定性
- 临床变量编码异质性需要谨慎解释

## 8. partial fusion 结果解释规范
- partial fusion 的目标是缓解弱配对导致的样本浪费，不等价于“更强模型”。
- 若 partial fusion 与主分析结果一致，可作为稳健性支持；若不一致，以主分析完整覆盖结果为主，并讨论差异来源。
- RPPA 相关结果在当前阶段定位为补充证据，不用于主结论主线。

## 9. baseline round1 已完成验证项（Phase 12-13）
- 已完成稳定性（首轮轻量版）：
  - repeated subsampling（6 次）
  - mean NMI / mean ARI / consensus_stability
  - 结果表：[results/tables/cluster_stability_round1.csv](results/tables/cluster_stability_round1.csv)
- 已完成内部发现性临床分析（首轮）：
  - OS 与 PFI 的 KM/log-rank
  - 首轮 Cox 摘要（可用时）
  - age/sex/stage 的关联检验
  - 结果表：[results/tables/clinical_association_round1.csv](results/tables/clinical_association_round1.csv)

## 10. 当前未完成与下一轮增强
1. methylation 尚未形成完整可建模矩阵，main5 baseline 仍待补齐。
2. 稳定性评估当前为轻量版本，后续需扩展到更完整的 consensus clustering 与 bootstrap。
3. 临床变量清洗与编码标准化仍需继续（尤其 Lauren 相关变量解析）。
4. 当前临床结果属于首轮发现性分析，不作为“临床有效性已证实”的结论。
