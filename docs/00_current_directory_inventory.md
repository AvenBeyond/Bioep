# 00 当前目录盘点与归档映射

## 1. 原始目录盘点（实际读取）
项目根目录下原始数据目录：
- [experiment_data](experiment_data)

原始文件/目录：
- [experiment_data/mc3_gene_level_STAD_mc3_gene_level.txt](experiment_data/mc3_gene_level_STAD_mc3_gene_level.txt)
  - 数据文件位于子路径：[experiment_data/mc3_gene_level_STAD_mc3_gene_level.txt/STAD_mc3_gene_level.txt](experiment_data/mc3_gene_level_STAD_mc3_gene_level.txt/STAD_mc3_gene_level.txt)
- [experiment_data/TCGA.STAD.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes](experiment_data/TCGA.STAD.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes)
  - 数据文件：[experiment_data/TCGA.STAD.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes](experiment_data/TCGA.STAD.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes)
- [experiment_data/TCGA.STAD.sampleMap_HiSeqV2](experiment_data/TCGA.STAD.sampleMap_HiSeqV2)
  - 数据文件：[experiment_data/TCGA.STAD.sampleMap_HiSeqV2/HiSeqV2](experiment_data/TCGA.STAD.sampleMap_HiSeqV2/HiSeqV2)
- [experiment_data/TCGA.STAD.sampleMap_HumanMethylation450](experiment_data/TCGA.STAD.sampleMap_HumanMethylation450)
  - 数据文件：[experiment_data/TCGA.STAD.sampleMap_HumanMethylation450/HumanMethylation450](experiment_data/TCGA.STAD.sampleMap_HumanMethylation450/HumanMethylation450)
- [experiment_data/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy](experiment_data/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy)
- [experiment_data/TCGA.STAD.sampleMap_miRNA_HiSeq_gene](experiment_data/TCGA.STAD.sampleMap_miRNA_HiSeq_gene)
  - 数据文件：[experiment_data/TCGA.STAD.sampleMap_miRNA_HiSeq_gene/miRNA_HiSeq_gene](experiment_data/TCGA.STAD.sampleMap_miRNA_HiSeq_gene/miRNA_HiSeq_gene)
- [experiment_data/TCGA.STAD.sampleMap_RPPA](experiment_data/TCGA.STAD.sampleMap_RPPA)
  - 数据文件：[experiment_data/TCGA.STAD.sampleMap_RPPA/RPPA](experiment_data/TCGA.STAD.sampleMap_RPPA/RPPA)
- [experiment_data/TCGA.STAD.sampleMap_STAD_clinicalMatrix](experiment_data/TCGA.STAD.sampleMap_STAD_clinicalMatrix)
- [experiment_data/survival_STAD_survival.txt](experiment_data/survival_STAD_survival.txt)

## 2. 新项目规范目录（已创建）
- [data/raw](data/raw)
- [data/interim](data/interim)
- [data/processed](data/processed)
- [src](src)
- [results](results)
- [tests](tests)
- [notebooks](notebooks)
- [docs](docs)
- [configs](configs)

## 3. 原始文件到 raw 归档映射（已复制）
- mutation:
  - 源：[experiment_data/mc3_gene_level_STAD_mc3_gene_level.txt/STAD_mc3_gene_level.txt](experiment_data/mc3_gene_level_STAD_mc3_gene_level.txt/STAD_mc3_gene_level.txt)
  - 目标：[data/raw/mutation/STAD_mc3_gene_level.txt](data/raw/mutation/STAD_mc3_gene_level.txt)
- cnv:
  - 源：[experiment_data/TCGA.STAD.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes](experiment_data/TCGA.STAD.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes)
  - 目标：[data/raw/cnv/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes](data/raw/cnv/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes)
- rna:
  - 源：[experiment_data/TCGA.STAD.sampleMap_HiSeqV2/HiSeqV2](experiment_data/TCGA.STAD.sampleMap_HiSeqV2/HiSeqV2)
  - 目标：[data/raw/rna/HiSeqV2](data/raw/rna/HiSeqV2)
- methylation:
  - 源：[experiment_data/TCGA.STAD.sampleMap_HumanMethylation450/HumanMethylation450](experiment_data/TCGA.STAD.sampleMap_HumanMethylation450/HumanMethylation450)
  - 目标：[data/raw/methylation/HumanMethylation450](data/raw/methylation/HumanMethylation450)
- methylation probe map:
  - 源：[experiment_data/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy](experiment_data/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy)
  - 目标：[data/raw/methylation/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy](data/raw/methylation/probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy)
- mirna:
  - 源：[experiment_data/TCGA.STAD.sampleMap_miRNA_HiSeq_gene/miRNA_HiSeq_gene](experiment_data/TCGA.STAD.sampleMap_miRNA_HiSeq_gene/miRNA_HiSeq_gene)
  - 目标：[data/raw/mirna/miRNA_HiSeq_gene](data/raw/mirna/miRNA_HiSeq_gene)
- rppa:
  - 源：[experiment_data/TCGA.STAD.sampleMap_RPPA/RPPA](experiment_data/TCGA.STAD.sampleMap_RPPA/RPPA)
  - 目标：[data/raw/rppa/RPPA](data/raw/rppa/RPPA)
- clinical:
  - 源：[experiment_data/TCGA.STAD.sampleMap_STAD_clinicalMatrix](experiment_data/TCGA.STAD.sampleMap_STAD_clinicalMatrix)
  - 目标：[data/raw/clinical/TCGA.STAD.sampleMap_STAD_clinicalMatrix](data/raw/clinical/TCGA.STAD.sampleMap_STAD_clinicalMatrix)
- survival:
  - 源：[experiment_data/survival_STAD_survival.txt](experiment_data/survival_STAD_survival.txt)
  - 目标：[data/raw/clinical/survival_STAD_survival.txt](data/raw/clinical/survival_STAD_survival.txt)

## 4. 命名不一致与潜在问题
- 部分数据放在目录中，真实数据文件位于子层级，易误复制空目录。
- 命名风格混杂：`sampleMap_` 前缀、全称文件名、无扩展名文件并存。
- clinical/survival 不是标准同构矩阵，后续读取逻辑需单独处理。
- 超大文件（例如 methylation）不适合直接在编辑器完整加载，需采用分块读取策略。

## 5. 本阶段操作记录
- 已完成目录创建与数据复制，不修改 [experiment_data](experiment_data) 原始数据。
- 已写入路径别名到 [configs/paths.yaml](configs/paths.yaml)。

## 6. 复制一致性与读取完整性（新增）
- 已输出复制一致性日志：[results/logs/raw_copy_verification.csv](results/logs/raw_copy_verification.csv)
  - 校验项：源/目标文件大小一致性、SHA256 哈希一致性。
- 已输出读取完整性日志：[results/logs/raw_file_inventory.csv](results/logs/raw_file_inventory.csv) 与 [results/logs/raw_file_inventory.json](results/logs/raw_file_inventory.json)
  - 校验项：全量行扫描状态、分隔符/编码猜测、总行数、列统计、预览信息与异常备注。
- 完整说明文档：[docs/11_data_reading_integrity_check.md](docs/11_data_reading_integrity_check.md)
