# 24 Partial Fusion 正式化结果

- 已实现 strict complete-case 与 partial cohort（至少4组学）的正式实验分支。
- partial 分支采用样本级可用组学重归一化权重，不因单组学缺失直接丢样本。
- 本轮在当前处理后的 main5 队列中，partial 主要提升了方法可扩展性，指标增益需以 `partial_fusion_metrics.csv` 为准。
- 若未超过 complete-case baseline，保持如实记录，不更改主结论。
