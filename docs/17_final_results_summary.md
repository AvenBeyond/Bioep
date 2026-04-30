# 17 最终结果摘要

统一定位：本项目最终交付为“基于多组学数据的胃腺癌分型与分型预测系统”，由 subtype discovery 与 subtype assignment 两层组成。

## 1. 最终主结论（当前版本）
- 最终主候选：`early_fusion_kmeans@main5,k=2`。
- 该选择是“当前数据与实现条件下更稳妥”的结论，而非宣称其绝对最优。

## 2. 方法比较结论
1. main5 baseline 相比 main4 提供了更完整组学输入，但并未保证所有方法都提升。
2. equal-weight fusion 在 main5 下出现明显退化簇现象（簇极不平衡）。
3. proposed weighted fusion 已完成正式运行，但在本轮设置下未稳定优于最佳 baseline。

## 3. 为什么最终选这个模型
- 选择依据来自：
	- baseline/proposed 直接对比
	- 消融排名
	- 临床分离指标
	- 分类器内部验证可行性
- 当前 proposed 的表现不够稳健，因此不作为最终主结论模型。

## 4. 创新点是否有效
- 有效性证据：消融显示权重项与组学组合会改变结果。
- 但“创新方法必然优于 baseline”在本轮不成立，应如实保留为待继续优化的问题。

## 5. 局限性与后续
- 详见 [docs/18_limitations_and_risk_statement.md](docs/18_limitations_and_risk_statement.md)。
- 后续优先项：partial fusion 正式接口、methylation 聚合增强、外部队列验证。
