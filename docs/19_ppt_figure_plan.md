# 19 PPT 图表计划（最终推荐主图顺序）

目标：仅保留一套答辩主图顺序，支持 8-10 页核心结果讲述。

## Page 1 背景与任务定义
- 文件名：results/figures/final_pipeline_flowchart.png
- 想证明什么：我们实现的是完整系统，不是零散实验。
- 一句话讲稿：本项目从多组学原始数据出发，形成了发现分型到新样本赋值的完整闭环。

## Page 2 样本队列与数据规模
- 文件名：results/figures/main_vs_with_rppa_intersection.png
- 想证明什么：主分析队列选择有数据依据。
- 一句话讲稿：主5组学交集样本更稳，加入 RPPA 会显著压缩队列规模。

## Page 3 主模型比较（总览）
- 文件名：results/figures/final_method_comparison_compact.png
- 想证明什么：主模型选择是多维指标平衡结果。
- 一句话讲稿：final baseline 在稳健性与可解释性上更均衡，创新分支作为增强证据保留。

## Page 4 聚类质量总览
- 文件名：results/figures/final_dashboard_overview.png
- 想证明什么：结果不是单一指标判断。
- 一句话讲稿：我们从聚类质量、鲁棒性、临床相关性、分类性能四个维度联合评估。

## Page 5 鲁棒性与可重复性
- 文件名：results/figures/baseline_main5_metric_comparison.png
- 想证明什么：不同方法和 K 的稳定差异真实存在。
- 一句话讲稿：main5 条件下并非所有方法都受益，主模型是在可重复性约束下确定。

## Page 6 临床意义
- 文件名：results/figures/final_km_curve_os.png
- 想证明什么：最终分型与生存结局存在可解释分离。
- 一句话讲稿：在当前内部评估中，最终分型对 OS 呈现显著分离信号。

## Page 7 创新与消融
- 文件名：results/figures/ablation_comparison.png
- 想证明什么：创新方法被真实验证而非主观宣称。
- 一句话讲稿：创新设计对结果有影响，但并未在本轮稳健超越主模型。

## Page 8 分类系统（assignment）
- 文件名：results/figures/subtype_classifier_performance.png
- 想证明什么：分型标签具备内部可预测性。
- 一句话讲稿：subtype assignment 模块可用，但当前定位仍是研究原型。

## Page 9 最终系统总结
- 文件名：results/figures/final_result_storyboard.png
- 想证明什么：项目交付已达到系统化与答辩级可读性。
- 一句话讲稿：我们最终交付的是可复现的胃腺癌分型与分型预测系统，而不是一组孤立实验。
