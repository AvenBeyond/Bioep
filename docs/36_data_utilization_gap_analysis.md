# 36 Data Utilization Gap Analysis

## 1. 当前主线已用到了什么
- main5 五组学 complete-case 矩阵（interim）
- baseline/proposed/improvement 结果链路
- OS/PFI 生存分离与部分临床变量关联
- subtype classifier 内部验证

## 2. 还缺什么
- 多终点（DSS/DFI）一致性证据
- 更系统的临床字段整合与分层分析
- 注释驱动的聚合特征（尤其 methylation）
- partial-fusion 正式扩容主线

## 3. 缺口类型判断
- 数据本身没有：外部独立队列（当前目录未见外部验证数据）
- 尚未充分利用：
  1. probe map 到聚合特征
  2. DSS/DFI 与治疗相关临床字段
  3. RPPA 辅助证据层
  4. weakly paired 样本扩容

## 4. 哪些缺口最值得补
优先级高：
1. at-least-4 partial-fusion 扩容
2. methylation 注释聚合
3. 临床与生存终点扩展（DSS/DFI + 关键临床字段）

对应资产依据：
- [results/tables/unused_or_underused_assets.csv](results/tables/unused_or_underused_assets.csv)
- [results/tables/sample_salvage_opportunities.csv](results/tables/sample_salvage_opportunities.csv)
- [results/tables/clinical_field_opportunities.csv](results/tables/clinical_field_opportunities.csv)
