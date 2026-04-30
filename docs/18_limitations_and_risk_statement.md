# 18 局限与风险声明

说明：以下风险声明针对“基于多组学数据的胃腺癌分型与分型预测系统”当前研究版本。

1. 无外部独立验证：
	- 当前全部评估基于 TCGA-STAD 内部数据，结论属于 internal evidence。
2. 样本量与簇分布限制：
	- 某些方法会产生极小簇，导致生存模型与分类器出现数值不稳定或不可定义指标。
3. methylation 处理策略局限：
	- 本轮采用 probe-level 缺失+方差筛选，尚未完成正式 gene/promoter 聚合，生物学解释仍有限。
4. graph connectivity 风险：
	- main5/proposed 运行中仍可见 spectral graph not fully connected warning，提示图结构可能影响谱聚类稳定性。
5. Cox 数值告警：
	- complete separation 与 ill-conditioned 相关警告提示部分生存比较存在估计不稳，结论应以发现性表达。
6. RPPA 未纳入主模型：
	- 这是样本覆盖与稳定性权衡结果，不代表 RPPA 无价值；其定位为补充轨道。
7. partial fusion 对照未完整实现：
	- 本轮已保留可追溯占位记录，完整版本仍需后续开发。
