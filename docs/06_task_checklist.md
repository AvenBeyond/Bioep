# 06 任务清单（持续更新）

状态说明：未开始 / 进行中 / 已完成 / 阻塞

## Phase 1 盘点与目录规范化
- [已完成] 检查原始目录与真实数据文件位置
- [已完成] 创建规范目录结构
- [已完成] 复制原始数据到 data/raw
- [已完成] 生成 docs/00
- [已完成] 生成 configs/paths.yaml
- [已完成] README 初稿

## Phase 2 数据概览与项目规划
- [已完成] 读取文件头与维度
- [已完成] 生成 docs/01
- [已完成] 生成 docs/02

## Phase 3 文献与创新设计
- [已完成] 在线检索 2024-2026 相关论文（首轮）
- [已完成] 生成 docs/03
- [已完成] 固化创新点与数学定义

## Phase 4 实验蓝图与验证协议
- [已完成] 生成 docs/04
- [已完成] 生成 docs/05
- [已完成] 生成 docs/09
- [已完成] 生成 docs/10

## Phase 5 工程骨架初始化
- [已完成] 创建 src 目录与模块骨架
- [已完成] 创建 tests 骨架
- [已完成] 创建 main.py
- [已完成] 创建 requirements.txt 与 .gitignore
- [已完成] main 状态检查逻辑实现

## Phase 6 第一轮自检
- [已完成] 原始数据读取完整性审计（streamed_scan + 哈希校验）
- [已完成] 生成 raw_file_inventory.csv/json 与 raw_copy_verification.csv
- [已完成] 生成 docs/11_data_reading_integrity_check.md
- [已完成] 安装最小运行依赖并运行 pytest + main smoke
- [已完成] 记录通过/失败项到 docs/08
- [已完成] 通过命令行完成本地虚拟环境初始化（不依赖 UI）
- [已完成] 记录 python/pip 版本与 pip list 到 environment_setup.log
- [已完成] 在 .venv 中执行 `python main.py` 与 `pytest -q` 并通过

## Phase 7 真实样本盘点与交集统计
- [已完成] 生成 modality_dimension_summary.csv
- [已完成] 生成 sample_inventory.csv（替代 placeholder）
- [已完成] 生成 sample_overlap_matrix.csv
- [已完成] 生成 multiomics_intersection_summary.csv
- [已完成] 生成 patient_alignment_rules.csv
- [已完成] 生成 3 张样本交集与样本量图

## Phase 8 主分析范围决策
- [已完成] 基于真实交集统计给出主分析队列定义
- [已完成] 基于数字证据明确 RPPA 定位（补充分析）
- [已完成] 更新 docs/04、docs/05、docs/07、docs/09

## Phase 9 preprocess 第一轮真实落地
- [已完成] mutation 真实预处理与 interim 输出
- [已完成] cnv 真实预处理与 interim 输出
- [已完成] rna 真实预处理与 interim 输出
- [已完成] mirna 真实预处理与 interim 输出
- [已完成] methylation 分块扫描 + 缺失统计 + preview 输出
- [已完成] rppa 补充组学预处理输出
- [已完成] 生成 preprocessing_dimension_changes.csv
- [已完成] 生成 docs/12_preprocessing_status.md

## Phase 10 第一轮阶段总结
- [已完成] 更新 README 当前状态与下一步
- [已完成] 更新 docs/06、docs/07、docs/08
- [已完成] 明确下一步进入相似网络构建与 baseline 初跑

## Phase 11 建模就绪性审计
- [已完成] 生成 modeling_readiness_summary.csv
- [已完成] 形成 docs/13 建模就绪性计划
- [已完成] 决定 round1 采用 main4（methylation 暂缓）

## Phase 12 相似网络构建 + baseline round1
- [已完成] 构建并保存各组学 similarity/affinity matrices
- [已完成] 运行 early_fusion（kmeans/spectral, K=2..6）
- [已完成] 运行 equal_weight_fusion（K=2..6）
- [已完成] 输出 labels / metrics / cluster sizes / fused summary / embedding 图

## Phase 13 稳定性与临床关联 round1
- [已完成] 计算 subsampling 稳定性（NMI/ARI/consensus）
- [已完成] 完成最佳 k 的 OS/PFI KM 首轮分析
- [已完成] 完成 age/sex/stage 的首轮关联检验

## Phase 14 proposed method readiness
- [已完成] 完成 weighted fusion 输入预览与权重估计 dry run
- [已完成] 生成 docs/14 proposed readiness 文档

## Phase 15 阶段总结
- [已完成] 更新 README
- [已完成] 更新 docs/06、docs/07、docs/08、docs/09、docs/13、docs/14
- [已完成] 明确下一步：补强 methylation + 正式 weighted fusion + 消融实验

## Phase 16 methylation 补强
- [已完成] 生成 methylation round2 可建模矩阵
- [已完成] 生成 methylation processing summary/filter log
- [已完成] 输出缺失与方差分布图

## Phase 17 main5 baseline 正式实验
- [已完成] 完成 main5 readiness 复核
- [已完成] 完成 main5 similarity 构建
- [已完成] 完成 main5 early/equal-weight baseline（K=2..6）

## Phase 18 weighted fusion 正式运行
- [已完成] 计算 Q/C/S/M 组件并生成 final weights
- [已完成] 运行 proposed weighted fusion（K=2..6）
- [已完成] 产出正式 proposed 指标、标签、簇规模、融合矩阵摘要

## Phase 19 系统消融实验
- [已完成] 方法消融（early/equal/proposed）
- [已完成] 权重项消融（no-C/no-S/no-M/only-Q/full）
- [已完成] 组学消融（drop methylation/mirna/mutation/cnv/rna）
- [已完成] partial fusion 对照以可追溯占位方式记录（接口待后续增强）

## Phase 20 subtype classifier 内部验证
- [已完成] 运行 RF/SVM/XGBoost（XGBoost 失败已如实记录）
- [已完成] 输出内部 CV 结果、混淆矩阵、特征重要性

## Phase 21 最终结果整合
- [已完成] 生成 final_model_comparison_summary.csv
- [已完成] 生成 final_key_numbers_for_ppt.csv
- [已完成] 形成 final results / limitations 文档

## Phase 22 答辩材料与提交版整理
- [已完成] 生成 PPT 图表规划与报告表格规划
- [已完成] 更新答辩讲述文档
- [已完成] 更新 README（round1 vs final 边界）
- [已完成] 新增 weighted/ablation/classifier smoke tests

## 当前优先级
1. 复核最终表述与答辩稿，确保“发现性结论”边界清晰
2. 若时间允许，补充 partial fusion 正式接口并替换占位对照
3. 作为扩展项：外部队列验证（GEO/ICGC）
4. 作为扩展项：methylation probe->gene/promoter 生物学聚合

## 需要人工确认（少数）
- 是否允许后续引入外部队列（GEO/ICGC）作为 stretch goal
- RPPA 在当前阶段不纳入主分析是否需要在答辩中固定为默认策略（建议是）

## Final 状态总览
- [已完成] final/report-ready 主线（Phase 16-22）
- [可选扩展] partial fusion 完整实现
- [stretch goal] 外部队列验证（GEO/ICGC）

## Phase 23-29 frontier-inspired 改进轮
- [已完成] Phase 23 前沿映射到实现计划
- [已完成] Phase 24 interaction-aware weighted fusion v2
- [已完成] Phase 25 partial-fusion 正式实验
- [已完成] Phase 26 consensus ensemble
- [已完成] Phase 27 graph connectivity-aware tuning
- [已完成] Phase 28 shared embedding prototype（exploratory）
- [已完成] Phase 29 改进轮统一汇总

## Final 系统化收口阶段
- [已完成] 统一项目定位为“胃腺癌分型与分型预测系统”
- [已完成] 生成系统总览文档（docs/29）
- [已完成] 生成最终 dashboard 总表与总览图板
- [已完成] 固定答辩主图顺序（8-10 张）
- [已完成] 新增最终入口脚本 run_final_system.py
- [已完成] 生成 assignment 说明、结果解释、课程核对清单（docs/30-32）
- [已完成] README 收口为项目总入口说明

## Targeted Enhancement Round
- [Done] Phase A methylation representation benchmark
- [Done] Phase B partial-fusion ge4 formal comparison
- [Done] Phase C multi-endpoint + extended clinical evidence
- [Done] Integrated replacement decision

## Targeted Enhancement Round
- [已完成] Phase A methylation 注释聚合增强
- [已完成] Phase B partial-fusion ge4 正式增强
- [已完成] Phase C DSS/DFI + 扩展临床字段补强
- [已完成] 统一对照与是否替换 final 的判定
