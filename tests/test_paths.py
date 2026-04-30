from pathlib import Path


def test_required_raw_files_exist():
    root = Path(__file__).resolve().parents[1]
    required = [
        root / "data" / "raw" / "mutation" / "STAD_mc3_gene_level.txt",
        root / "data" / "raw" / "cnv" / "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",
        root / "data" / "raw" / "rna" / "HiSeqV2",
        root / "data" / "raw" / "methylation" / "HumanMethylation450",
        root / "data" / "raw" / "methylation" / "probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy",
        root / "data" / "raw" / "mirna" / "miRNA_HiSeq_gene",
        root / "data" / "raw" / "rppa" / "RPPA",
        root / "data" / "raw" / "clinical" / "TCGA.STAD.sampleMap_STAD_clinicalMatrix",
        root / "data" / "raw" / "clinical" / "survival_STAD_survival.txt",
    ]
    for p in required:
        assert p.exists(), f"Missing required raw data file: {p}"
