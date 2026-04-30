# 23 前沿改进映射（frontier-inspired）

## 借鉴并落地的前沿思想
- weighted multi-kernel learning：转译为现有 affinity 的质量加权 + 交互项加权，不引入重型优化器。
- omics-omics interaction-aware integration：在 v2 增加组学对交互矩阵与交互权重。
- weakly paired multi-omics clustering：将 partial fusion 从占位升级为正式重加权实验。
- graph-aware clustering：对谱聚类构图参数做系统网格调优，而非仅记录 warning。
- lightweight shared embedding：采用浅层降维+共享投影原型，作为探索分支。

## 不适合当前项目的方向
- 端到端深图神经网络与大规模对比学习：样本规模与课程周期下复现风险高。
- 高复杂度贝叶斯非参数多视图模型：调参与解释成本高。

## 最终实现的增量升级
- Phase 24: interaction-aware weighted fusion v2
- Phase 25: partial-fusion formal experiments
- Phase 26: consensus ensemble subtyping
- Phase 27: graph connectivity-aware tuning
- Phase 28: shared-embedding prototype（探索）

## 模块映射
- `src/clustering/proposed_weighted_fusion_v2.py`: v2 交互加权融合
- `src/clustering/partial_fusion_experiments.py`: partial-fusion 正式化
- `src/clustering/consensus_ensemble.py`: 共识分型
- `src/feature_engineering/graph_connectivity_tuning.py`: 图连通性调优
- `src/clustering/shared_embedding_prototype.py`: 轻量共享嵌入

## 预期收益与风险
- 预期收益：提高稳定性与簇平衡，在不牺牲可解释性的前提下争取临床分离改进。
- 风险：部分改进可能只提升稳定性不提升生存分离；partial cohort 受上游预处理样本覆盖约束。
