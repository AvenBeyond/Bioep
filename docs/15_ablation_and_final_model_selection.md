# 15 消融与最终模型选择

## 1. 消融范围
已完成：
- 方法消融：early fusion / equal-weight fusion / proposed weighted fusion
- 权重项消融：no-C / no-S / no-M / only-Q / full
- 组学消融：drop methylation / mirna / mutation / cnv / rna

结果文件：
- [results/tables/ablation_results.csv](results/tables/ablation_results.csv)
- [results/tables/ablation_rank_summary.csv](results/tables/ablation_rank_summary.csv)

## 2. 当前可支持的结论
1. 消融已显示不同权重项与组学组合对结果有显著影响，创新模块并非“无条件提升”。
2. 在本轮数据与实现下，proposed full 并未稳定优于最佳 baseline，应如实保留为“真实结果但不夸大”。
3. 当前最终主候选更偏向稳健性优先（main5 early_fusion_kmeans@k=2）。

## 3. 对答辩可直接使用的话术
- “我们没有把 proposed 当作先验最优，而是通过系统消融让数据告诉我们哪些组件在当前 cohort 下有效。”
- “若 proposed 未稳定领先，我们也不强行下结论，这保证了研究可信度。”
