# 08 进度日志

## 2026-04-17 | Phase 1 完成
- 已完成内容：
  - 盘点原始目录结构
  - 创建规范目录树
  - 将 9 个原始数据文件复制至 data/raw 对应子目录
  - 创建路径配置与 README 初稿
- 产出文件：
  - [docs/00_current_directory_inventory.md](docs/00_current_directory_inventory.md)
  - [configs/paths.yaml](configs/paths.yaml)
  - [README.md](README.md)
- 问题与备注：
  - 多个数据以“目录+内部文件”形式存在，已按真实数据文件复制

## 2026-04-17 | Phase 2 完成
- 已完成内容：
  - 统计各数据文件行列规模
  - 抽取头部预览并识别主要编码模式
  - 形成数据问题与预处理建议
- 产出文件：
  - [docs/01_data_overview.md](docs/01_data_overview.md)
  - [docs/02_project_master_plan.md](docs/02_project_master_plan.md)

## 2026-04-17 | Phase 3 完成
- 已完成内容：
  - 在线检索 2024-2026 相关文献（PubMed 首轮）
  - 形成 8 篇相关方法/场景论文总结
  - 固化创新点、公式与答辩 QA
- 产出文件：
  - [docs/03_literature_and_innovation_design.md](docs/03_literature_and_innovation_design.md)
- 问题与备注：
  - DOI 与深度全文解读待下一轮补充

## 2026-04-17 | Phase 4 完成
- 已完成内容：
  - 实验蓝图、验证协议、答辩笔记、分工模板建立
- 产出文件：
  - [docs/04_experiment_blueprint.md](docs/04_experiment_blueprint.md)
  - [docs/05_validation_and_testing_protocol.md](docs/05_validation_and_testing_protocol.md)
  - [docs/09_defense_notes.md](docs/09_defense_notes.md)
  - [docs/10_group_contribution_template.md](docs/10_group_contribution_template.md)

## 2026-04-17 | Phase 5 完成
- 已完成内容：
  - 创建 src 模块骨架
  - 创建 tests 基础测试
  - 创建 main.py 统一入口
  - 创建 default/path 配置、依赖与忽略规则
- 产出文件：
  - [main.py](main.py)
  - [src](src)
  - [tests](tests)
  - [configs/default.yaml](configs/default.yaml)
  - [requirements.txt](requirements.txt)

## 2026-04-17 | Phase 6 进行中
- 已完成内容：
  - 准备运行 smoke test
- 待完成内容：
  - 执行 pytest
  - 记录通过/失败项
  - 更新 docs/06 与 docs/08

## 2026-04-17 | 完整性审计机制补充完成
- 已完成内容：
  - 新增原始数据完整性扫描模块：[src/data_loading/raw_integrity_checks.py](src/data_loading/raw_integrity_checks.py)
  - 新增主入口审计命令：`python main.py --integrity`
  - 对 9 个 raw 文件完成流式全量扫描（`read_mode=streamed_scan`）
  - 对 9 个源/目标文件完成大小一致性与 SHA256 校验（9/9 通过）
  - 新增预处理维度变更日志接口并接入 6 个 preprocess 骨架
- 产出文件：
  - [results/logs/raw_file_inventory.csv](results/logs/raw_file_inventory.csv)
  - [results/logs/raw_file_inventory.json](results/logs/raw_file_inventory.json)
  - [results/logs/raw_copy_verification.csv](results/logs/raw_copy_verification.csv)
  - [docs/11_data_reading_integrity_check.md](docs/11_data_reading_integrity_check.md)
- 问题与处理：
  - 依赖安装阶段遇到 `xgboost.dll` 文件占用导致整包安装中断，已先安装最小运行依赖（pandas/PyYAML/pytest）并修复 pandas 版本冲突。
- 后续接续点：
  - 完成 pytest 与 main smoke 自检闭环并记录结果。

## 2026-04-17 | Phase 6 自检闭环完成
- 已完成内容：
  - 执行 `python main.py` 状态检查：raw 数据发现 9/9
  - 执行 `python -m pytest -q`：5 passed
- 产出文件：
  - [results/tables/raw_data_inventory.csv](results/tables/raw_data_inventory.csv)
  - [results/tables/sample_inventory_placeholder.csv](results/tables/sample_inventory_placeholder.csv)
- 当前结论：
  - 第一轮“目录规范化 + 文档先行 + 代码骨架 + 完整性校验 + 基础自检”已闭环。

## 2026-04-17 | 命令行虚拟环境初始化（补充要求）
- 已完成内容：
  - 在 `D:\Bioep` 通过命令行执行环境初始化流程，不依赖 UI。
  - Python 可用性检查：`python command: C:\Users\huawei\AppData\Local\Programs\Python\Python312\python.exe`，版本 `Python 3.12.4`。
  - `.venv` 检查结果：已存在，按要求跳过创建步骤。
  - 虚拟环境激活：成功（`Activate.ps1` 执行成功）。
  - 基础工具升级：`pip/setuptools/wheel` 已升级（pip=26.0.1, setuptools=82.0.1, wheel=0.46.3）。
  - 依赖安装：`requirements.txt` 安装完成（含 lifelines/matplotlib/seaborn）。
  - 执行验证：`python main.py` 退出码 0；`pytest -q` 结果 `5 passed`。
- 产出文件：
  - [results/logs/environment_setup.log](results/logs/environment_setup.log)
- 备注：
  - 本次未触发创建/激活失败分支，因此失败命令与替代方案日志条目为“无失败”。

## 2026-04-17 | Phase 7 真实样本盘点与交集统计完成
- 已完成内容：
  - 新增并执行样本对齐脚本：[src/data_loading/sample_alignment_analysis.py](src/data_loading/sample_alignment_analysis.py)
  - 生成正式表格：
    - [results/tables/modality_dimension_summary.csv](results/tables/modality_dimension_summary.csv)
    - [results/tables/sample_inventory.csv](results/tables/sample_inventory.csv)
    - [results/tables/sample_overlap_matrix.csv](results/tables/sample_overlap_matrix.csv)
    - [results/tables/multiomics_intersection_summary.csv](results/tables/multiomics_intersection_summary.csv)
    - [results/tables/patient_alignment_rules.csv](results/tables/patient_alignment_rules.csv)
  - 生成正式图形：
    - [results/figures/sample_overlap_heatmap.png](results/figures/sample_overlap_heatmap.png)
    - [results/figures/modality_sample_counts_barplot.png](results/figures/modality_sample_counts_barplot.png)
    - [results/figures/main_vs_with_rppa_intersection.png](results/figures/main_vs_with_rppa_intersection.png)
- 核心数字：
  - 主5组学交集=364
  - 含RPPA 6组学交集=288
  - partial fusion 覆盖（主5）>=3组学=436，>=4组学=391
  - 主5加入RPPA后交集下降 76（约20.9%）

## 2026-04-17 | Phase 8 主分析范围决策完成
- 已完成内容：
  - 明确主分析队列：主5组学交集 Primary Tumor 样本。
  - 明确 RPPA 定位：补充分析、消融分析、生物学解释补充，不纳入主分析。
  - 决策与依据已写入 [docs/07_decision_log.md](docs/07_decision_log.md)。

## 2026-04-17 | Phase 9 preprocess 第一轮真实落地完成
- 已完成内容：
  - 新增运行入口：[src/preprocessing/run_round1_preprocessing.py](src/preprocessing/run_round1_preprocessing.py)
  - mutation/cnv/rna/mirna/rppa 完成真实可执行预处理并输出 interim 文件。
  - methylation 完成 chunked_read 全量扫描统计与 top-variance preview，不夸大为完整聚合。
  - 生成维度变更日志：[results/logs/preprocessing_dimension_changes.csv](results/logs/preprocessing_dimension_changes.csv)
- 产出文件：
  - [data/interim/mutation_round1.csv](data/interim/mutation_round1.csv)
  - [data/interim/cnv_round1.csv](data/interim/cnv_round1.csv)
  - [data/interim/rna_round1.csv](data/interim/rna_round1.csv)
  - [data/interim/mirna_round1.csv](data/interim/mirna_round1.csv)
  - [data/interim/rppa_round1.csv](data/interim/rppa_round1.csv)
  - [data/interim/methylation_round1_summary.csv](data/interim/methylation_round1_summary.csv)
  - [data/interim/methylation_round1_preview_topvar.csv](data/interim/methylation_round1_preview_topvar.csv)

## 2026-04-17 | Phase 10 第一轮阶段总结
- 本阶段结论：
  - 工程检查阶段已完成，当前已进入真实数据分析准备阶段。
  - 样本对齐与主分析范围已定稿，预处理第一轮已落地并可复现运行。
- 下一步：
  - 使用 round1 interim 数据进入相似网络构建与 baseline 初跑。

## 2026-04-18 | Phase 11 建模就绪性审计完成
- 已完成内容：
  - 生成 [results/tables/modeling_readiness_summary.csv](results/tables/modeling_readiness_summary.csv)
  - 审计结论：methylation 当前仅 summary/preview，不是完整建模矩阵。
  - round1 baseline 组学组合确定为 main4（mutation+cnv+rna+mirna）。

## 2026-04-18 | Phase 12 相似网络构建 + baseline round1 完成
- 已完成内容：
  - 构建 aligned matrices 与 similarity/affinity matrices。
  - 运行 baseline：early_fusion(kmeans/spectral) + equal_weight_fusion，K=2..6。
  - 输出 labels、metrics、cluster_sizes、fused summary、embedding 图。
- 产出文件（核心）：
  - [results/tables/modality_modeling_shapes.csv](results/tables/modality_modeling_shapes.csv)
  - [results/tables/similarity_build_log.csv](results/tables/similarity_build_log.csv)
  - [results/tables/baseline_round1_metrics.csv](results/tables/baseline_round1_metrics.csv)
  - [results/tables/baseline_round1_cluster_sizes.csv](results/tables/baseline_round1_cluster_sizes.csv)
  - [results/tables/fused_similarity_summary.csv](results/tables/fused_similarity_summary.csv)
- 备注：
  - spectral 过程中出现 `Graph is not fully connected` 警告，流程可完成但提示图连通性风险，已在后续计划中记录。

## 2026-04-18 | Phase 13 稳定性与临床关联 round1 完成
- 已完成内容：
  - 稳定性（subsampling 6次）：NMI/ARI/consensus。
  - 候选最佳 k 的 OS/PFI KM 与首轮临床变量关联分析。
- 产出文件：
  - [results/tables/cluster_stability_round1.csv](results/tables/cluster_stability_round1.csv)
  - [results/tables/clinical_association_round1.csv](results/tables/clinical_association_round1.csv)
  - [results/figures/clinical_association_heatmap_round1.png](results/figures/clinical_association_heatmap_round1.png)

## 2026-04-18 | Phase 14 proposed readiness 完成
- 已完成内容：
  - 生成创新方法输入预览与权重 dry run（非最终结果）。
- 产出文件：
  - [results/tables/proposed_weight_inputs_preview.csv](results/tables/proposed_weight_inputs_preview.csv)
  - [results/tables/proposed_weight_estimates_round0.csv](results/tables/proposed_weight_estimates_round0.csv)
  - [docs/14_proposed_method_readiness.md](docs/14_proposed_method_readiness.md)

## 2026-04-18 | Phase 15 阶段总结
- 本轮结论：
  - baseline 已真实跑通，首轮指标、稳定性与临床关联已落盘。
  - 当前最稳妥的 round1 主候选为 `early_fusion_kmeans@k=2`（稳定性与簇平衡更好）。
  - equal_weight_fusion 出现极端不平衡簇，暂不作为主结论。
- 下一步接续点：
  - 补齐 methylation 可建模矩阵并升级到 main5 baseline。
  - 在 main5 上推进正式 weighted fusion 与消融实验。

## 2026-04-18 | 回归测试补记
- 运行命令：
  - `.venv\\Scripts\\python.exe -m pytest -q`
- 首次结果：
  - 1 个失败（`tests/test_similarity_build.py`），原因为 matplotlib 默认 Tk 后端在当前环境缺失 `tk.tcl`。
- 修复措施：
  - 在绘图模块统一切换到 `Agg` 后端（headless 兼容）。
- 二次结果：
  - `7 passed, 19 warnings`（全部通过）。

## 2026-04-18 | Phase 16 methylation 补强完成
- 已完成内容：
  - 实现并运行 methylation round2 可建模流程（两阶段 chunk 处理）。
  - 形成缺失率过滤 + 高方差筛选后的 sample-by-feature 矩阵。
- 关键输出：
  - [data/interim/methylation_round2_modeling.csv](data/interim/methylation_round2_modeling.csv)
  - [results/tables/methylation_processing_summary.csv](results/tables/methylation_processing_summary.csv)
  - [results/tables/methylation_feature_filter_log.csv](results/tables/methylation_feature_filter_log.csv)
  - [results/figures/methylation_missingness_summary.png](results/figures/methylation_missingness_summary.png)
  - [results/figures/methylation_variance_distribution.png](results/figures/methylation_variance_distribution.png)
- 结论：
  - `final_matrix_ready=True`，methylation 已可进入 main5 正式实验。

## 2026-04-18 | Phase 17 main5 baseline 正式实验完成
- 已完成内容：
  - 完成 main5 readiness 复核。
  - 重建 main5 similarity/affinity。
  - 运行 early_fusion 与 equal_weight_fusion（K=2..6）。
- 核心输出：
  - [results/tables/modeling_readiness_main5.csv](results/tables/modeling_readiness_main5.csv)
  - [results/tables/similarity_build_log_main5.csv](results/tables/similarity_build_log_main5.csv)
  - [results/tables/baseline_main5_metrics.csv](results/tables/baseline_main5_metrics.csv)
  - [results/tables/baseline_main5_cluster_sizes.csv](results/tables/baseline_main5_cluster_sizes.csv)
  - [results/tables/fused_similarity_summary_main5.csv](results/tables/fused_similarity_summary_main5.csv)
- 观察：
  - main5 不等价于“必然提升”；equal-weight fusion 出现明显退化分簇。

## 2026-04-18 | Phase 18 weighted fusion 正式运行完成
- 已完成内容：
  - 计算 Q_m/C_m/S_m/M_m 并得到正式权重。
  - 执行 proposed weighted fusion（K=2..6）。
- 输出：
  - [results/tables/proposed_weight_components.csv](results/tables/proposed_weight_components.csv)
  - [results/tables/proposed_weight_estimates_final.csv](results/tables/proposed_weight_estimates_final.csv)
  - [results/tables/proposed_main5_metrics.csv](results/tables/proposed_main5_metrics.csv)
  - [results/tables/proposed_main5_cluster_sizes.csv](results/tables/proposed_main5_cluster_sizes.csv)
- 如实记录：
  - 本轮 proposed 在当前设置下并未优于最佳 baseline，保持为“结果真实但不强结论”。

## 2026-04-18 | Phase 19 消融实验完成
- 已完成内容：
  - 方法消融、权重项消融、组学消融。
  - partial fusion 对照因接口未完全落地，保留可追溯占位并写入 notes。
- 输出：
  - [results/tables/ablation_results.csv](results/tables/ablation_results.csv)
  - [results/tables/ablation_rank_summary.csv](results/tables/ablation_rank_summary.csv)
  - [results/figures/ablation_comparison.png](results/figures/ablation_comparison.png)

## 2026-04-18 | Phase 20 分类器内部验证完成
- 已完成内容：
  - 训练 RF/SVM/XGBoost（XGBoost 在当前标签分布下失败并已记录）。
  - 采用 repeated stratified k-fold 内部验证（含泄漏安全 pipeline）。
- 输出：
  - [results/tables/subtype_classifier_cv_results.csv](results/tables/subtype_classifier_cv_results.csv)
  - [results/tables/subtype_classifier_feature_importance.csv](results/tables/subtype_classifier_feature_importance.csv)
  - [results/tables/subtype_assignment_confusion_matrix.csv](results/tables/subtype_assignment_confusion_matrix.csv)
- 风险记录：
  - 部分类别样本过少导致 AUROC 不稳定，已在结果表 notes 体现。

## 2026-04-18 | Phase 21-22 最终整合与答辩包装完成
- 已完成内容：
  - 生成最终比较与答辩关键数字。
  - 生成 15-20 号文档并更新答辩材料规划。
- 输出：
  - [results/tables/final_model_comparison_summary.csv](results/tables/final_model_comparison_summary.csv)
  - [results/tables/final_key_numbers_for_ppt.csv](results/tables/final_key_numbers_for_ppt.csv)
  - [docs/15_ablation_and_final_model_selection.md](docs/15_ablation_and_final_model_selection.md)
  - [docs/16_subtype_classifier_validation.md](docs/16_subtype_classifier_validation.md)
  - [docs/17_final_results_summary.md](docs/17_final_results_summary.md)
  - [docs/18_limitations_and_risk_statement.md](docs/18_limitations_and_risk_statement.md)
  - [docs/19_ppt_figure_plan.md](docs/19_ppt_figure_plan.md)
  - [docs/20_report_table_plan.md](docs/20_report_table_plan.md)

## 2026-04-18 | warning 与数值稳定性记录（必须披露）
- spectral graph connectivity warning 在 main5/proposed 过程中仍存在。
- Cox 在部分簇划分上出现 complete separation / ill-conditioned 相关告警。
- 处理方式：
  - 不忽略；在结论与 limitations 中明确写入，解释为“发现性证据，需谨慎外推”。

## 2026-04-18 | 最终回归测试（Phase 16-22 后）
- 执行命令：
  - `.venv\\Scripts\\python.exe -m pytest -q`
- 结果：
  - `10 passed, 40 warnings`（全部通过）
- 说明：
  - warning 主要来自 spectral graph connectivity、lifelines 绘图弃用提示、以及 classifier smoke 中 `k>n_features` 的提示。
  - 当前作为课程研究版本可复现通过，warning 已在 limitations 与 defense 文档中披露。

## 2026-04-18 19:52 | Phase 23-29 改进轮完成
- 已完成：frontier-inspired mapping、interaction-aware v2、partial formalization、consensus ensemble、graph tuning、shared embedding prototype。
- 关键产物：`results/tables/improvement_model_comparison.csv`、`results/tables/improvement_key_takeaways.csv`、`docs/28_improvement_round_summary.md`。
- 原 final 主线未覆盖，全部新增结果采用 improvement/v2 独立命名。

## 2026-04-18 | 改进轮结论对齐
- 结论：shared embedding prototype 在稳定性与簇平衡方面更稳健，但整体仍不足以替换既有 final baseline 主结论。
- 如实记录：interaction-aware v2 与 partial equal/weighted fusion 未体现稳健净提升。

## 2026-04-18 | Final 系统化收口完成
- 已完成：final system 统一口径、dashboard 总表与四张总览图、最终入口脚本、assignment 说明与答辩解释文档。
- 新增核心产物：
  - `results/tables/final_dashboard_master.csv`
  - `results/figures/final_dashboard_overview.png`
  - `results/figures/final_method_comparison_compact.png`
  - `results/figures/final_pipeline_flowchart.png`
  - `results/figures/final_result_storyboard.png`
  - `docs/29_system_overview.md` 到 `docs/32_requirement_checklist.md`

## 2026-04-20 Progress
- Targeted enhancement finished. Selected candidate: old_final_baseline (k=2).

## 2026-04-20 | Targeted Enhancement Stable Survival Pass
- 已停止全量不稳定 Cox 路径，改为分层精算策略。
- 已完成 3 candidates × OS/PFI 的真实 KM + log-rank 主证据产出。
- Cox 采用门控 + penalized 单次回退，失败即降级为 NA，不再卡住整轮。
- DSS/DFI：在稳定条件下已补算。

## 2026-04-20 | Targeted Enhancement Round 完成
- 严格执行 A->B->C 顺序：methylation 聚合增强、partial ge4 正式增强、多终点临床补强。
- 未新增复杂模型家族，未覆盖旧 final 结果。
- 本轮选中候选：phaseA_best_methylation (k=2)。
- 详见 docs/38 与 targeted_enhancement_comparison.csv。

## 2026-04-20 | Targeted Enhancement Stable Survival Pass
- 已停止全量不稳定 Cox 路径，改为分层精算策略。
- 已完成 3 candidates × OS/PFI 的真实 KM + log-rank 主证据产出。
- Cox 采用门控 + penalized 单次回退，失败即降级为 NA，不再卡住整轮。
- DSS/DFI：在稳定条件下已补算。
