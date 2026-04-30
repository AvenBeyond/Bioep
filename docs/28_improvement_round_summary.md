# 28 改进轮总结（Phase 23-29）

统一定位：改进轮属于“基于多组学数据的胃腺癌分型与分型预测系统”的增强分支，不单独定义为新主线。

## 借鉴到的前沿思路
- 多核/多视图加权、组学交互项、弱配对融合、图连通性感知调优、轻量共享嵌入。

## 已实现且可复现的改进
- interaction-aware weighted fusion v2
- partial fusion formal experiments
- consensus ensemble subtyping
- graph connectivity-aware tuning
- shared embedding prototype（exploratory）

## 有效与无效改进
- 以 `improvement_model_comparison.csv` 与 `improvement_key_takeaways.csv` 为准。
- 若改进未超过既有 final baseline，保持原主结论不变。
- 本轮相对稳健的改进候选为 shared embedding prototype（稳定性与簇平衡更好）。
- interaction-aware v2 与 partial equal/weighted fusion 未体现稳健净增益。

## 是否替换原 final 主模型
- 不替换。当前仍保留既有 final baseline 作为主结果。
- 改进轮结果作为“前沿启发的增量证据”在答辩中呈现。

## 答辩建议表述
- 本轮改进强调“可解释 + 可复现 + 可比较”，并对未增益分支如实披露，避免过度结论。
