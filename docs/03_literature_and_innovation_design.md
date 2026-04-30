# 03 文献综述与创新设计（重点文档）

## A. 胃腺癌多组学分型背景
胃腺癌（STAD）具有显著分子异质性，单组学难以全面描述肿瘤生物学。TCGA 提供 mutation/CNV/RNA/methylation/RPPA 等多模态数据，为发现稳定且具有临床意义的亚型提供基础。

## B. 为什么单组学或简单拼接不够
- 单组学：仅覆盖特定生物层，可能遗漏关键机制。
- 简单拼接（early fusion）：
  - 忽略各组学噪声水平和可靠性差异
  - 易受高维组学主导
  - 对缺失组学样本不友好
因此需要“可解释的加权融合 + 弱配对友好策略 + 稳健 K 选择”。

## C. 2024-2026 相关文献综述（首轮）
说明：本轮通过 PubMed E-utilities 在线检索获取条目。后续需补充 DOI 与全文细读笔记。

### 1) SeOMLR: one-step multi-view latent representation with self-weighted ensemble learning for multi-omics cancer subtyping
- 年份：2026
- 来源：Bioinformatics (Oxford, England)
- 任务：多组学癌症分型
- 核心方法：自加权多视图集成学习
- 优点：显式建模视图权重，贴近本项目加权思想
- 局限：方法复杂度高，对参数和实现细节依赖大
- 借鉴：可作为“质量加权融合”理论支撑，实际实现可采用更轻量权重公式

### 2) Survival-Informed Multi-Omics Kernel Fusion for Cancer Subtyping
- 年份：2026
- 来源：IEEE TCBBS
- 任务：结合生存信息的多组学亚型识别
- 核心方法：生存引导核融合
- 优点：强调临床终点导向，和本项目“临床相关权重”一致
- 局限：核方法参数较多，计算和解释门槛偏高
- 借鉴：将“临床分离度 C_m”纳入权重项而非端到端复杂优化

### 3) Adaptive multi-view information bottleneck for multi-omics data clustering
- 年份：2026
- 来源：Briefings in Bioinformatics
- 任务：多视图聚类
- 核心方法：信息瓶颈 + 自适应多视图融合
- 优点：多视图权重自适应，鲁棒性较强
- 局限：模型抽象层较高，对课程项目解释成本偏大
- 借鉴：保留“自适应权重”思想，放弃重型深度框架

### 4) Decoupled contrastive multi-view clustering with adaptive false negative elimination for cancer subtyping
- 年份：2025
- 来源：PLoS Computational Biology
- 任务：多视图对比聚类用于癌症分型
- 核心方法：对比学习 + 多视图聚类
- 优点：在视图一致性方面表现好
- 局限：训练稳定性和调参成本较高
- 借鉴：作为“前沿对照方法”讨论，不作为主线实现

### 5) GL-Fusion: graph-level structure fusion and locus-level feature fusion for cancer subtype classification
- 年份：2025
- 来源：IEEE JBHI
- 任务：多组学亚型分类
- 核心方法：图结构融合 + 位点级融合
- 优点：对多尺度信息建模细致
- 局限：工程复杂，课程周期内复现风险较高
- 借鉴：图结构思路可用于未来扩展，不进入第一轮主实现

### 6) Cross-cohort multi-omics analysis identifies novel clusters driven by EMT signatures in gastric cancer
- 年份：2026
- 来源：Cancer Cell International
- 任务：胃癌多组学聚类与生物学解释
- 核心方法：跨队列整合 + EMT 相关分群
- 优点：胃癌场景相关性强，具备临床叙事价值
- 局限：跨队列数据准备复杂，当前课程项目难完整复刻
- 借鉴：在结果讨论中增强胃癌特异性生物学解释

### 7) Integrative proteogenomics maps multifactorial aetiology, progression and therapeutic vulnerabilities in gastric cancer
- 年份：2026
- 来源：Gut
- 任务：胃癌多组学机制图谱
- 核心方法：蛋白组学 + 基因组学整合
- 优点：高影响力、临床转化导向强
- 局限：数据与方法体系重，不适合课程级端到端复刻
- 借鉴：用于支撑“多组学整合必要性”与“临床意义”论证

### 8) Machine Learning-Based Prognostic Model for Gastric Cancer Using Integrated Multi-Omics Data
- 年份：2025
- 来源：Cancer Investigation
- 任务：胃癌预后建模
- 核心方法：多组学融合 + 机器学习
- 优点：与课程项目规模接近，方法可落地
- 局限：偏监督预后，不直接等价无监督分型
- 借鉴：用于 subtype assignment 内部验证模块设计

## D. 本项目最终创新点（固定版本）
1. 创新点 1：临床引导的质量加权多组学相似网络融合
2. 创新点 2：弱配对样本友好的部分融合策略
3. 创新点 3：兼顾聚类质量、稳定性和临床区分度的多目标 K 选择机制

## E. 数学定义与流程

### 1) 组学权重定义
对于组学 $m$，定义权重：

$$
w_m = softmax(\alpha Q_m + \beta C_m + \delta S_m - \gamma M_m)
$$

其中：
- $Q_m$：组学内聚类质量
- $C_m$：临床相关性得分
- $S_m$：样本覆盖度
- $M_m$：缺失惩罚

### 2) 融合相似网络

$$
S_{fused} = \sum_m w_m S_m
$$

### 3) 多目标 K 评分

$$
Score(K)=0.35\cdot Silhouette + 0.25\cdot ConsensusStability + 0.20\cdot ClinicalSeparation + 0.20\cdot SubtypeBalance
$$

### 4) 方法流程（答辩版）
1. 各组学独立预处理与相似网络构建
2. 计算每组学质量/临床/覆盖/缺失指标
3. 基于 softmax 计算权重并融合网络
4. 若样本弱配对，启用 partial fusion fallback
5. 在 K=2..6 搜索并用多目标函数选最优 K
6. 输出亚型并做临床/稳定性评估

## F. 为什么不直接采用重型端到端深度学习
- 课程项目关注可解释、可复现、可答辩，不是纯 SOTA 竞赛。
- TCGA-STAD 样本规模有限，重模型高方差风险高。
- 深度模型调参成本高，短周期内稳定复现难度大。
- 本项目方案保留创新性，同时降低工程风险。

## G. 答辩问答草稿（至少 10 组）
1. 问：为什么不用简单拼接？
   答：简单拼接默认各组学同质量，且对缺失样本不友好。我们通过权重与 partial fusion 提升鲁棒性。
2. 问：权重怎么保证不是拍脑袋？
   答：权重由 $Q_m,C_m,S_m,M_m$ 四项组成，均可被量化并可做消融验证。
3. 问：没有外部测试集如何证明可靠？
   答：明确区分 discovery 与 internal validation，采用 repeated resampling、consensus 稳定性、nested CV 等减少乐观偏差。
4. 问：为什么 K 范围选 2..6？
   答：结合 STAD 样本规模和课程可解释性，避免过细分导致小类不稳定。
5. 问：为什么 RPPA 默认不进主分析？
   答：RPPA 样本覆盖通常更低，先作为可选组学，待交集统计后再决定。
6. 问：创新点和现有文献区别？
   答：我们将“临床引导 + 缺失惩罚 + 覆盖度”统一到轻量可解释权重公式中，并增加课程可落地的多目标 K 机制。
7. 问：如何证明 partial fusion 有价值？
   答：设置启用/禁用 partial fusion 的对照与稳定性比较。
8. 问：如何避免信息泄漏？
   答：分类器验证阶段所有 scaler/feature selection/model fitting 均在训练折内完成。
9. 问：结果最关键看什么？
   答：同时看聚类质量、稳定性、生存区分度和亚型平衡度，而非单一指标。
10. 问：项目局限性是什么？
    答：缺乏独立外部验证；部分组学缺失较多；结果偏发现性，需后续外部队列验证。

## H. 可直接放入 PPT 的创新点精炼总结
- 我们不是“简单拼接多组学”，而是“先评估每组学可信度，再做临床引导加权融合”。
- 我们显式处理 TCGA 常见弱配对问题，避免仅保留全交集样本导致样本浪费。
- 我们用多目标 K 选择替代单指标选 K，让亚型结果更稳定、更临床可解释。

## I. 参考文献（首轮，后续补 DOI）
1. SeOMLR. Bioinformatics, 2026.
2. Survival-Informed Multi-Omics Kernel Fusion for Cancer Subtyping. IEEE TCBBS, 2026.
3. Adaptive multi-view information bottleneck for multi-omics data clustering. Briefings in Bioinformatics, 2026.
4. Decoupled contrastive multi-view clustering for cancer subtyping. PLoS Comput Biol, 2025.
5. GL-Fusion for cancer subtype classification. IEEE JBHI, 2025.
6. Cross-cohort multi-omics EMT clusters in gastric cancer. Cancer Cell International, 2026.
7. Integrative proteogenomics in gastric cancer. Gut, 2026.
8. ML prognostic model using integrated multi-omics in gastric cancer. Cancer Investigation, 2025.
