# 16 Subtype Classifier Internal Validation

## 1. 目标与边界
- 目标：验证最终候选 subtype labels 在内部数据上的可预测性稳定性。
- 边界：仅为 internal validation，不是 external test，也不构成临床部署结论。

## 2. 方法设置
- 模型：RandomForest / SVM / XGBoost（可用时）。
- 验证：RepeatedStratifiedKFold（训练折内完成预处理与特征选择，避免信息泄漏）。
- 输出：
	- [results/tables/subtype_classifier_cv_results.csv](results/tables/subtype_classifier_cv_results.csv)
	- [results/tables/subtype_classifier_feature_importance.csv](results/tables/subtype_classifier_feature_importance.csv)
	- [results/tables/subtype_assignment_confusion_matrix.csv](results/tables/subtype_assignment_confusion_matrix.csv)

## 3. 本轮实际情况
1. RF 与 SVM 完成评估并给出稳定指标。
2. XGBoost 在当前标签分布与折分条件下失败，已在结果 notes 明确标注。
3. AUROC 在类别极少折上可能不可定义，已按 NaN/警告处理。

## 4. 结论表达建议
- 该部分用于支撑“分型 assignment 具有内部可预测性”。
- 不应写成“模型具有外部泛化能力”。
