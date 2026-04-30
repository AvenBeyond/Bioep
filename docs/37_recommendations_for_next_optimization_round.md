# 37 Recommendations for Next Optimization Round

## 1. 建议清单（3-5件事）
详表：
- [results/tables/top3_next_steps_recommended.csv](results/tables/top3_next_steps_recommended.csv)

### 1. 推进 at-least-4 partial-fusion 正式扩容并统一与 complete-case 的评估框架
- 目标层：sample
- 预期收益：以较低异质性代价增加样本，提升稳定性与统计功效
- 风险：模态缺失引入偏差，需要缺失机制控制
- 成本：medium

更可能改善指标：
- stability / clinical relevance / evidence strength

### 2. methylation 从 probe-level 升级到 promoter/gene-level 注释聚合特征
- 目标层：feature
- 预期收益：降低噪声并提升解释性，改善聚类与分类稳健性
- 风险：聚合策略不当可能损失局部信号
- 成本：medium

更可能改善指标：
- cluster quality / classification performance / interpretability

### 3. 将 DSS/DFI + 关键临床变量纳入多终点一致性与分层证据链
- 目标层：evidence
- 预期收益：显著增强结论完整性与说服力
- 风险：终点缺失值与事件数可能限制显著性
- 成本：low_medium

更可能改善指标：
- clinical relevance / evidence strength

## 2. 结论
下一轮最值得优先做的是：先做 at-least-4 partial-fusion 扩容与统一评估框架，同时并行推进 methylation 注释聚合与 DSS/DFI 终点补强。