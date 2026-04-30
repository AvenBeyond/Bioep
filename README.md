# TCGA-STAD 多组学胃腺癌分型与分型预测系统

## 项目简介
本项目最终交付为：基于多组学数据的胃腺癌分型与分型预测系统。

系统包含两层：
1. subtype discovery（无监督分型发现）
2. subtype assignment（新样本分型预测，当前为研究原型）

方法体系包含：
- 传统基线：early fusion、equal-weight fusion
- 创新方法：proposed weighted fusion
- 改进分支：frontier-inspired methods（interaction-aware、consensus、graph tuning、shared embedding）

## 当前目录说明
- [configs](configs): 配置文件
- [data](data): 规范化数据目录（raw/interim/processed）
- [docs](docs): 项目文档体系
- [src](src): 主代码骨架
- [tests](tests): 基础工程测试
- [results](results): 输出目录（图表、日志、模型）
- [notebooks](notebooks): 探索与可视化补充
- [experiment_data](experiment_data): 原始来源目录（只读，不覆盖）

## 数据路径说明
已将数据从 [experiment_data](experiment_data) 复制归档到 [data/raw](data/raw)：
- mutation: [data/raw/mutation/STAD_mc3_gene_level.txt](data/raw/mutation/STAD_mc3_gene_level.txt)
- cnv: [data/raw/cnv/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes](data/raw/cnv/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes)
- rna: [data/raw/rna/HiSeqV2](data/raw/rna/HiSeqV2)
- methylation: [data/raw/methylation/HumanMethylation450](data/raw/methylation/HumanMethylation450)
- methylation probe map: [data/raw/methylation/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy](data/raw/methylation/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy)
- mirna: [data/raw/mirna/miRNA_HiSeq_gene](data/raw/mirna/miRNA_HiSeq_gene)
- rppa: [data/raw/rppa/RPPA](data/raw/rppa/RPPA)
- clinical: [data/raw/clinical/TCGA.STAD.sampleMap_STAD_clinicalMatrix](data/raw/clinical/TCGA.STAD.sampleMap_STAD_clinicalMatrix)
- survival: [data/raw/clinical/survival_STAD_survival.txt](data/raw/clinical/survival_STAD_survival.txt)

## 当前进度
第一轮已完成：
- 目录规范化与原始数据归档
- 文档体系 00-10 创建
- 创新点与验证协议固化
- 代码骨架初始化
- main 状态检查与基础测试框架

第二轮（数据分析准备）已完成：
- 真实样本盘点与交集统计（Phase 7）
- 主分析范围与 RPPA 定位决策（Phase 8）
- preprocess 第一轮真实落地（Phase 9）
- 阶段总结与文档同步（Phase 10）

第三轮（baseline round1）已完成：
- 建模就绪性审计（Phase 11）
- 相似网络构建与 baseline 首轮运行（Phase 12）
- 稳定性与临床关联首轮（Phase 13）
- proposed method readiness（round0 dry run, Phase 14）
- 阶段总结（Phase 15）

第四轮（final/report-ready）已完成：
- methylation round2 可建模矩阵补强（Phase 16）
- main5 baseline 正式实验（Phase 17）
- proposed weighted fusion 正式运行（Phase 18）
- 系统消融实验（Phase 19）
- subtype classifier 内部验证（Phase 20）
- 最终结果整合与答辩材料固化（Phase 21-22）

关键数字：
- 主5组学交集样本数：364
- 加入 RPPA 后6组学交集样本数：288
- 纳入 RPPA 导致交集下降：76（约20.9%）

进度详情见 [docs/08_progress_log.md](docs/08_progress_log.md)。

## 如何运行
1. 安装依赖
```bash
pip install -r requirements.txt
```
2. 运行项目状态检查
```bash
python main.py
```
3. 运行原始数据完整性审计（推荐先执行）
```bash
python main.py --integrity
```
4. 运行基础测试
```bash
pytest -q
```

5. 运行 Phase 7 样本对齐与交集统计
```bash
python -m src.data_loading.sample_alignment_analysis
```

6. 运行 Phase 9 round1 预处理
```bash
python -m src.preprocessing.run_round1_preprocessing
```

7. 运行 Phase 11-15 baseline round1 管线
```bash
python -m src.pipelines.run_phase11_15
```

8. 运行 Phase 16-22 final/report-ready 管线
```bash
python -m src.pipelines.run_phase16_22
```

完整性审计产物：
- [results/logs/raw_file_inventory.csv](results/logs/raw_file_inventory.csv)
- [results/logs/raw_file_inventory.json](results/logs/raw_file_inventory.json)
- [results/logs/raw_copy_verification.csv](results/logs/raw_copy_verification.csv)
- [docs/11_data_reading_integrity_check.md](docs/11_data_reading_integrity_check.md)

说明：头部预览仅用于结构识别，不代表已完成全量读取。后续 docs 中涉及维度/样本数的描述，应以完整性日志为依据。

样本对齐与预处理关键产物：
- [results/tables/modality_dimension_summary.csv](results/tables/modality_dimension_summary.csv)
- [results/tables/sample_inventory.csv](results/tables/sample_inventory.csv)
- [results/tables/sample_overlap_matrix.csv](results/tables/sample_overlap_matrix.csv)
- [results/tables/multiomics_intersection_summary.csv](results/tables/multiomics_intersection_summary.csv)
- [results/logs/preprocessing_dimension_changes.csv](results/logs/preprocessing_dimension_changes.csv)
- [docs/12_preprocessing_status.md](docs/12_preprocessing_status.md)

baseline round1 关键产物：
- [results/tables/modeling_readiness_summary.csv](results/tables/modeling_readiness_summary.csv)
- [results/tables/baseline_round1_metrics.csv](results/tables/baseline_round1_metrics.csv)
- [results/tables/cluster_stability_round1.csv](results/tables/cluster_stability_round1.csv)
- [results/tables/clinical_association_round1.csv](results/tables/clinical_association_round1.csv)
- [docs/13_modeling_readiness_and_baseline_plan.md](docs/13_modeling_readiness_and_baseline_plan.md)
- [docs/14_proposed_method_readiness.md](docs/14_proposed_method_readiness.md)

final/report-ready 关键产物：
- [results/tables/methylation_processing_summary.csv](results/tables/methylation_processing_summary.csv)
- [results/tables/modeling_readiness_main5.csv](results/tables/modeling_readiness_main5.csv)
- [results/tables/baseline_main5_metrics.csv](results/tables/baseline_main5_metrics.csv)
- [results/tables/proposed_main5_metrics.csv](results/tables/proposed_main5_metrics.csv)
- [results/tables/ablation_results.csv](results/tables/ablation_results.csv)
- [results/tables/subtype_classifier_cv_results.csv](results/tables/subtype_classifier_cv_results.csv)
- [results/tables/final_model_comparison_summary.csv](results/tables/final_model_comparison_summary.csv)
- [docs/17_final_results_summary.md](docs/17_final_results_summary.md)
- [docs/18_limitations_and_risk_statement.md](docs/18_limitations_and_risk_statement.md)

结果边界说明：
- round1/exploratory：`baseline_round1_*`, `proposed_weight_estimates_round0.csv`
- final/report-ready：`*_main5_*`, `proposed_weight_estimates_final.csv`, `ablation_*`, `subtype_classifier_*`, `final_*`

placeholder 状态：
- [results/tables/sample_inventory_placeholder.csv](results/tables/sample_inventory_placeholder.csv) 保留为历史文件，不再作为正式引用。
- 正式清单文件为 [results/tables/sample_inventory.csv](results/tables/sample_inventory.csv)。

## 下一步计划
1. 扩展项：补齐 partial fusion 正式接口与对照
2. 扩展项：methylation probe->gene/promoter 聚合验证
3. stretch goal：外部队列验证（GEO/ICGC）
4. 持续优化：图连通性与生存模型数值稳定性


第五轮（frontier-inspired improvement）已完成：
- Phase 23-29 在不覆盖 final 主线前提下新增 v2/improvement 实验。
- 关键汇总见 `results/tables/improvement_model_comparison.csv` 与 `docs/28_improvement_round_summary.md`。
- 运行命令：`python -m src.pipelines.run_phase23_29`。
- 当前结论：改进轮结果作为增强证据保留，原 final baseline 主结论不变。

## 最终主结论
- 当前最终主模型：early_fusion_kmeans@main5,k=2。
- 主分析组学：mutation/cnv/methylation/rna/mirna。
- 创新方法与改进方法的角色：提供增量证据与边界探索，不强行替换主模型。

## 最终系统入口
统一入口脚本：`python -m src.pipelines.run_final_system`

支持模式：
1. `--summary`
	- 输出最终主模型、关键指标、classifier 概要、主要局限。
2. `--report`
	- 自动生成并检查最终 dashboard 表与总览图。
3. `--predict-subtype --input <csv> [--output <csv>]`
	- 新样本 subtype assignment prototype（研究原型，不用于临床诊断）。

## 最终图表与文档导航
系统总览：
- [docs/29_system_overview.md](docs/29_system_overview.md)

分型预测说明：
- [docs/30_how_subtype_assignment_works.md](docs/30_how_subtype_assignment_works.md)

最终结果解释：
- [docs/31_final_result_interpretation.md](docs/31_final_result_interpretation.md)

课程要求核对：
- [docs/32_requirement_checklist.md](docs/32_requirement_checklist.md)

最终 dashboard 产物：
- [results/tables/final_dashboard_master.csv](results/tables/final_dashboard_master.csv)
- [results/figures/final_dashboard_overview.png](results/figures/final_dashboard_overview.png)
- [results/figures/final_method_comparison_compact.png](results/figures/final_method_comparison_compact.png)
- [results/figures/final_pipeline_flowchart.png](results/figures/final_pipeline_flowchart.png)
- [results/figures/final_result_storyboard.png](results/figures/final_result_storyboard.png)

## 边界声明
- 最终主结果：基于当前内部证据确定，可复现、可答辩。
- 研究原型部分：新样本预测接口、外部泛化能力与临床部署仍需后续验证。

## Round 6: Targeted Enhancement Completed
- Added results/tables/targeted_enhancement_comparison.csv
- Added docs/38_targeted_enhancement_summary.md

## Round 6 Stable Survival Update
- Targeted enhancement survival validation now uses stable mode (KM + log-rank first).
- Cox is conditional and may be skipped for separation/instability.
- New outputs: `results/tables/multi_endpoint_clinical_validation_stable.csv`, `results/tables/targeted_enhancement_survival_summary.csv`.
- Real figures: `results/figures/targeted_km_os_real.png`, `results/figures/targeted_km_pfi_real.png`, `results/figures/multi_endpoint_validation_comparison_stable.png`.

## 第六轮（targeted enhancement）已完成
- 已按 atlas Top3 路线完成 A/B/C 三阶段增强，不新增复杂模型家族。
- 新增核心对照：`results/tables/targeted_enhancement_comparison.csv`。
- 新增总结文档：`docs/38_targeted_enhancement_summary.md`。

## Round 6 Stable Survival Update
- Targeted enhancement survival validation now uses stable mode (KM + log-rank first).
- Cox is conditional and may be skipped for separation/instability.
- New outputs: `results/tables/multi_endpoint_clinical_validation_stable.csv`, `results/tables/targeted_enhancement_survival_summary.csv`.
- Real figures: `results/figures/targeted_km_os_real.png`, `results/figures/targeted_km_pfi_real.png`, `results/figures/multi_endpoint_validation_comparison_stable.png`.
