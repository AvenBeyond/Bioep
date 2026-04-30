# 30 分型预测模块说明（subtype assignment）

## 1. 为什么系统分成“发现分型”和“分型预测”两层
在多组学任务中，先做 discovery 再做 assignment 是标准工程路径：
- discovery 负责回答“是否存在可解释、可重复的分型结构”。
- assignment 负责回答“对新样本是否能复用已有分型规则”。

如果跳过 discovery，assignment 将失去标签来源与生物学解释基础。

## 2. 当前最终用于 subtype assignment 的模型
当前内部验证表现最佳的分类器为 SVM（main5_concat 特征）。

它使用 discovery 阶段最终主模型生成的 subtype 标签作为监督信号。

## 3. 输入需要哪些组学
默认期望主5组学特征：mutation/cnv/methylation/rna/mirna。

在 prototype 入口中，新样本输入文件应尽量与训练特征模板一致（包含组学前缀列名）。

## 4. 输入缺失组学怎么办
当前 assignment 入口可作为 prototype 运行，但存在前置约束：
- 需要与训练模板有足够特征重叠（脚本当前要求至少 50 个重叠特征）。
- 若重叠不足，系统会拒绝预测并提示输入不满足条件。

这代表当前版本是“研究级可控输入”，不是“任意临床输入即插即用”。

## 5. 输出是什么
prototype 输出至少包含：
- sample_id
- predicted_subtype
- model_used
- disclaimer（研究原型声明）

输出位置示例：results/tables/prototype_predicted_subtypes.csv。

## 6. 当前 classifier 的内部验证表现
以现有结果为准：
- SVM：accuracy 约 0.94，macro-F1 约 0.94（内部交叉验证）

应强调这是内部证据，不等价于外部泛化性能。

## 7. 为什么这仍然是研究原型，而不是临床部署工具
主要原因：
- 缺外部独立验证与前瞻性验证
- 输入标准化与缺失组学处理尚未临床级产品化
- 结果解释与风险控制流程尚未满足临床合规要求

因此当前可称为：胃腺癌多组学分型与分型预测研究系统原型。
