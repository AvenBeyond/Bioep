# 34 单一图文成果汇报（最终提交版）

## 报告定位
本文件为本项目唯一成果汇报文档，整合最终结论、关键图表与逐图解释。

报告主题：基于多组学数据的胃腺癌分型与分型预测系统（研究原型）。

---

## 1. 一句话结论
在当前数据与验证边界下，最终主模型为 early_fusion_kmeans@main5,k=2；该分型在内部评估中表现为“稳定性较高、临床分离显著、可用于内部预测”，创新改进分支提供了增量证据但不替代主线结论。

---

## 2. 核心指标总览

来源表：
- [results/tables/final_dashboard_master.csv](results/tables/final_dashboard_master.csv)
- [results/tables/final_model_comparison_summary.csv](results/tables/final_model_comparison_summary.csv)
- [results/tables/subtype_classifier_cv_results.csv](results/tables/subtype_classifier_cv_results.csv)

关键数值（最终口径）：
- 最终主模型：early_fusion_kmeans（main5, k=2）
- silhouette：0.0527
- stability：0.8566
- OS log-rank p：0.00140
- PFI log-rank p：0.0240
- 最优分类器（内部）：SVM，macro-F1=0.9377，accuracy=0.9396

---

## 3. 图1：四维综合能力概览

![图1-四维综合能力概览](../results/figures/final_dashboard_overview.png)

图文件：
- [results/figures/final_dashboard_overview.png](results/figures/final_dashboard_overview.png)

图意说明：
- 该图将 clustering quality、robustness、clinical relevance、classification performance 四个维度做归一化展示。
- clinical relevance 柱高明显，原因是 OS 的 p 值很小，转换后信号强，说明分型与生存结局关联明确。
- robustness 保持较高，体现为主模型在重采样场景下稳定性较好。
- clustering quality 相对保守，提示该系统不是依赖“几何分离极强”取胜，而是依赖多证据平衡。

结论解释：
- 该图支持“综合最优而非单指标最优”的模型选择逻辑。

---

## 4. 图2：方法对比（紧凑版）

![图2-方法对比紧凑图](../results/figures/final_method_comparison_compact.png)

图文件：
- [results/figures/final_method_comparison_compact.png](results/figures/final_method_comparison_compact.png)

图意说明：
- 对 final baseline、proposed weighted fusion、improvement best 三类候选进行多维归一化对比。
- final baseline 在稳定性和临床维度综合更均衡。
- improvement best 在部分结构性指标上有亮点，但整体尚未形成稳健全面超越。
- proposed weighted fusion 在当前轮次未表现出替代主线的综合优势。

结论解释：
- 该图给出“为什么最终仍选 baseline”的直接证据链。

---

## 5. 图3：最终 OS 生存曲线

![图3-最终OS生存曲线](../results/figures/final_km_curve_os.png)

图文件：
- [results/figures/final_km_curve_os.png](results/figures/final_km_curve_os.png)

图意说明：
- 两个亚型的生存曲线存在持续分离趋势，提示分型具有临床风险分层意义。
- 结合 OS p=0.00140，可认为在内部数据上具有统计学显著差异。

结论解释：
- 分型结果不仅“可聚类”，还具有临床关联解释力，这是最终保留主模型的重要依据。

风险说明：
- 当前证据来自内部数据，外部队列泛化能力尚未验证。

---

## 6. 图4：分类器性能对比

![图4-分类器性能对比](../results/figures/subtype_classifier_performance.png)

图文件：
- [results/figures/subtype_classifier_performance.png](results/figures/subtype_classifier_performance.png)

图意说明：
- SVM 在 accuracy 与 macro-F1 上均优于 RandomForest。
- XGBoost 在当前流程中未形成稳定可用结果。

结论解释：
- discovery 产生的 subtype 标签具备可学习结构，assignment 在内部验证下可行。

---

## 7. 图5：分类混淆矩阵

![图5-分类混淆矩阵](../results/figures/subtype_classifier_confusion_matrix.png)

图文件：
- [results/figures/subtype_classifier_confusion_matrix.png](results/figures/subtype_classifier_confusion_matrix.png)

图意说明：
- SVM 的错分在两类中分布较均衡，没有明显偏向单一亚型。
- 这说明模型没有出现严重类偏置，分类边界整体稳定。

结论解释：
- 内部预测性能可用于“研究原型级”新样本赋值流程。

---

## 8. 图6：特征重要性（Top20）

![图6-特征重要性Top20](../results/figures/subtype_feature_importance_top20.png)

图文件：
- [results/figures/subtype_feature_importance_top20.png](results/figures/subtype_feature_importance_top20.png)

图意说明：
- 头部特征贡献高于长尾特征，呈现“少数强特征 + 一组辅助特征”结构。
- 该结果可用于后续生物学解释或 marker 候选筛查。

结论解释：
- 模型并非纯黑箱，具备一定可解释性基础。

风险说明：
- 重要性不等于因果，后续仍需独立验证与机制分析。

---

## 9. 图7：创新改进分支综合对比

![图7-创新改进分支对比](../results/figures/final_innovation_comparison.png)

图文件：
- [results/figures/final_innovation_comparison.png](results/figures/final_innovation_comparison.png)
- [results/tables/final_innovation_model_report.csv](results/tables/final_innovation_model_report.csv)
- [results/tables/final_innovation_key_takeaways.csv](results/tables/final_innovation_key_takeaways.csv)

图意说明：
- improvement_best（shared_embedding_prototype）在部分结构性指标中表现良好。
- proposed_weighted_fusion 当前轮次未达到替代主线的综合效果。

结论解释：
- 创新改进分支的价值是“增量证据与方向探索”，而非强行替换主结论。

---

## 10. 交付口径（可直接朗读）
我们最终交付了一个“多组学分型发现 + 新样本分型预测”的完整研究系统原型。
在当前内部证据下，主模型确定为 early_fusion_kmeans@main5,k=2，并在稳定性、临床分离和可预测性上形成了可复现的综合优势。
创新改进分支已被系统评估并保留为增强证据，但当前不替代主线结论。

---

## 11. 当前边界与下一步
- 当前边界：缺少外部独立队列验证，尚不构成临床部署级证据。
- 下一步建议：
  1. 增加外部队列验证（跨平台数据）
  2. 固化持久化分类模型导出与加载接口
  3. 继续推进关键特征的生物学解释与验证

---

## 12. 本文档与其他文档关系
本文件为最终“单一成果汇报”文档。
若仅提交一份报告，请提交本文件。