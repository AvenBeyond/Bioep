# 14 Proposed Method Readiness（Round0）

## 1. 当前目标定位
本阶段不是最终版 weighted fusion 实验，而是为创新方法做输入准备与权重预估 dry run。

## 2. 已完成输入准备
已具备输入：
- 各组学相似矩阵/affinity 矩阵（main4）
  - 目录：[results/tables/similarity_matrices](results/tables/similarity_matrices)
- 建模样本对齐与组学维度审计
  - [results/tables/modeling_readiness_summary.csv](results/tables/modeling_readiness_summary.csv)
  - [results/tables/modality_modeling_shapes.csv](results/tables/modality_modeling_shapes.csv)
- baseline round1 指标可作为 `Q_m` 先验来源
  - [results/tables/baseline_round1_metrics.csv](results/tables/baseline_round1_metrics.csv)

## 3. 权重框架与当前实现
已推进代码：
- [src/clustering/proposed_weighted_fusion.py](src/clustering/proposed_weighted_fusion.py)

当前框架：
$$
 w_m = \mathrm{softmax}(\alpha Q_m + \beta C_m + \delta S_m - \gamma M_m)
$$

round0 输入预览：
- [results/tables/proposed_weight_inputs_preview.csv](results/tables/proposed_weight_inputs_preview.csv)

round0 权重估计：
- [results/tables/proposed_weight_estimates_round0.csv](results/tables/proposed_weight_estimates_round0.csv)

说明：
- 当前为 prototype/dry run，不应当作为最终实验结果对外宣称。

## 4. 当前缺口
1. `C_m`（临床相关项）目前仍是简化代理值，需与正式临床分离度评估绑定。
2. `Q_m` 目前是轻量 proxy，需要替换为更稳健的模态质量度量（如模态内稳定性 + 与生存分离一致性）。
3. partial fusion 仍需接入真实缺失图结构与样本覆盖策略，不仅是“可用即加权”。
4. 主5（含 methylation）版本尚未 ready，当前仅 main4 权重预估。

## 5. 距离正式运行还差什么
1. 完成 methylation 可建模矩阵并并入 main5。
2. 完整落地 `Q_m/C_m/S_m/M_m` 计算链并固定参数来源。
3. 将 weighted fusion 与 baseline 在同一 cohort 上做公平比较（含稳定性与临床端点）。
4. 增加消融：去临床项、去 partial fusion、去某组学。

## 6. 下一轮建议
1. 先补强 methylation 可建模矩阵，再跑 main5 baseline 对照。
2. 用 main5 与 main4 双轨验证 proposed 权重策略对稳定性和临床分离度的影响。
3. 将 round0 权重预估升级为 round1 可解释实验结果。
