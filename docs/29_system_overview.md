# 29 最终系统总览

## 系统目标
本项目最终交付为：基于多组学数据的胃腺癌分型与分型预测系统。

系统目标不是只输出单次聚类结果，而是形成可复现的两层流程：
1. subtype discovery：在既有队列中发现稳定、可解释、具临床相关性的分型。
2. subtype assignment：基于已发现分型标签，训练内部预测器用于新样本 subtype 赋值（研究原型）。

## 输入数据类型
系统当前主分析输入为 TCGA-STAD 主5组学：
- mutation
- cnv
- methylation
- rna
- mirna

临床矩阵与生存数据用于分型临床相关性评估（如 log-rank）。

## 数据预处理模块
预处理模块负责：
- 样本 ID 规范化与跨组学对齐
- 组学内缺失值与尺度处理
- methylation 的缺失率过滤与方差筛选（round2 建模矩阵）
- interim 矩阵与过程日志落盘

这一步保证后续 discovery 与 assignment 共用同一数据边界，避免比较不公平。

## 分型发现模块（subtype discovery）
分型发现模块包含：
- baseline 方法：early fusion（kmeans/spectral）、equal-weight fusion
- 创新方法：proposed weighted fusion
- 改进方法：interaction-aware v2、consensus、graph tuning、shared embedding 等
- 统一评价：聚类质量、稳定性、临床相关性、簇平衡、消融证据

系统当前最终主结果保持为 main5 下 early_fusion_kmeans (k=2)。

## 最终主模型
当前主模型：early_fusion_kmeans@main5,k=2。

选择原因是“综合稳健性最优”，而不是单一指标最高：
- 聚类质量处于可接受范围
- 稳定性与簇平衡更稳
- OS/PFI 生存分离在当前内部评估下具有统计意义
- 支持后续 subtype assignment 内部验证

## 创新方法与改进方法的位置
创新与改进方法在系统中的角色是“增强证据”和“边界探索”，而非强制替换主模型：
- proposed weighted fusion：验证了加权机制可解释性，但本轮未稳定超越主模型
- frontier-inspired 改进（Phase 23-29）：用于验证交互融合、partial/consensus/graph/shared-embedding 的可行性与潜在增益

因此，系统结论是“主线稳定 + 创新可追溯”，不是“为了创新而替换结果”。

## subtype classifier 模块（subtype assignment）
assignment 模块基于 discovery 得到的最终 subtype 标签进行内部训练与验证：
- RandomForest / SVM / XGBoost（失败分支如实记录）
- 重复分层交叉验证
- 输出准确率、macro-F1、混淆矩阵与特征重要性

该模块回答的是“分型标签是否可预测”，用于支持系统完整性。

## 输出结果形式
系统输出包括：
- 结构化结果表：模型比较、关键数字、dashboard 总表
- 答辩级图表：总览图、方法比较图、流程图、storyboard
- 运行入口：run_final_system.py（summary/report/predict-subtype prototype）
- 解释文档：系统说明、结果解读、课程要求核对

## 适用范围与局限性
适用范围：
- 课程项目与研究原型验证
- 多组学分型与分型预测流程演示

当前局限：
- 仍缺外部独立验证队列
- 一些谱方法分支存在图连通性风险
- subtype assignment 目前为研究原型，不是临床部署工具
- 新样本预测对输入特征格式有前置约束
