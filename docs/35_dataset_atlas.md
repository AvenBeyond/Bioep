# 35 Dataset Atlas（数据资产深度盘点）

## 1. 盘点结论概览
- 原始层（raw）共识别数据文件 16 个（含 data/raw 与 experiment_data 的可读数据文件）。
- 其中老师原始包目录（experiment_data）可读数据文件：8 个。
- 其中规范化复制目录（data/raw）可读数据文件：8 个。
- 原始层中已进入 final 主模型直接链路的文件数：12。
- 识别到未充分利用或可进一步利用资产：55 项。
- 当前 strict main5 complete-case 样本数：364。
- 可通过 partial-fusion at-least-4 潜在挽回样本：27。
- 可通过 partial-fusion at-least-3 潜在挽回样本：72。

## 2. 原始数据全景
核心组学与临床资源：
- mutation / cnv / methylation / rna / mirna / rppa
- clinical matrix / survival（含 OS/PFI/DSS/DFI）
- methylation probe map 注释资源

对应全清单见：
- [results/tables/data_asset_inventory.csv](results/tables/data_asset_inventory.csv)

## 3. 当前项目实际使用链路
raw -> interim -> final 的链路图：
![data asset map](../results/figures/data_asset_map.png)

关键说明：
- final 主线主要依赖 main5（mutation/cnv/methylation/rna/mirna）interim 矩阵。
- rppa 原始数据存在且可用，但因样本交集代价在 final 主线中未纳入。
- 临床与生存数据已用于 OS/PFI 相关验证，但未完全展开到 DSS/DFI 与更多临床字段。

## 4. 未充分利用资产画像
详表见：
- [results/tables/unused_or_underused_assets.csv](results/tables/unused_or_underused_assets.csv)

高价值未充分利用方向：
1. methylation probe 注释到 gene/promoter/region 的聚合特征
2. survival 的 DSS/DFI 终点
3. clinical matrix 的治疗/复发/分级字段
4. partial-fusion 扩容样本空间
5. rppa 作为辅助证据层

## 5. 样本挽回空间
![sample salvage](../results/figures/sample_salvage_opportunity_chart.png)

对应表：
- [results/tables/sample_salvage_opportunities.csv](results/tables/sample_salvage_opportunities.csv)

解释：
- 在不改动主线结论前提下，at-least-4 是最稳妥扩容入口。
- at-least-3 提供更大样本空间，但异质性与缺失机制风险更高。

## 6. 特征增强潜力
![feature opportunity](../results/figures/feature_opportunity_priority_chart.png)

对应表：
- [results/tables/feature_enrichment_opportunities.csv](results/tables/feature_enrichment_opportunities.csv)

解释：
- methylation 聚合特征、RNA/miRNA pathway 表达、mutation/CNV 汇总特征是下一轮重点。

## 7. 临床字段使用缺口
![clinical usage](../results/figures/clinical_field_usage_map.png)

对应表：
- [results/tables/clinical_field_opportunities.csv](results/tables/clinical_field_opportunities.csv)

解释：
- 当前已用字段相对有限，仍有大量可用临床变量可增强证据链完整性。

## 8. 最值得挖掘的资源
1. methylation 注释驱动聚合（噪声控制 + 解释增强）
2. partial-fusion at-least-4 样本扩容（提升统计功效）
3. DSS/DFI + 治疗相关临床字段（增强证据充分性）
