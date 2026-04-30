"""Phase 23-29 frontier-inspired improvement runner.

This pipeline adds improvement/v2 experiments without overwriting final results.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

from src.clustering.consensus_ensemble import run_consensus_ensemble
from src.clustering.partial_fusion_experiments import run_partial_fusion_experiments
from src.clustering.proposed_weighted_fusion import estimate_final_weights
from src.clustering.proposed_weighted_fusion_v2 import run_weighted_fusion_v2
from src.clustering.shared_embedding_prototype import run_shared_embedding_prototype
from src.evaluation.clinical_association import evaluate_clinical_association
from src.evaluation.cluster_metrics import estimate_stability_by_subsampling
from src.feature_engineering.build_similarity_matrices import euclidean_affinity
from src.feature_engineering.graph_connectivity_tuning import run_graph_connectivity_tuning
from src.utils.config_utils import load_project_config


MAIN5_MODALITIES = ["mutation", "cnv", "methylation", "rna", "mirna"]


def _normalize_sample_id(sample_id: str) -> str:
    sid = str(sample_id).strip().upper()
    return sid[:16] if len(sid) >= 16 else sid


def _read_interim(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).map(_normalize_sample_id)
    return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def _load_main5_matrices(project_root: Path) -> dict[str, pd.DataFrame]:
    interim = project_root / "data" / "interim"
    mats = {
        "mutation": _read_interim(interim / "mutation_round1.csv"),
        "cnv": _read_interim(interim / "cnv_round1.csv"),
        "methylation": _read_interim(interim / "methylation_round2_modeling.csv"),
        "rna": _read_interim(interim / "rna_round1.csv"),
        "mirna": _read_interim(interim / "mirna_round1.csv"),
    }
    cohort = sorted(set.intersection(*[set(m.index) for m in mats.values()]))
    return {m: df.loc[cohort].copy() for m, df in mats.items()}


def _load_or_build_affinities(project_root: Path, mats: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    sim_dir = project_root / "results" / "tables" / "similarity_matrices_main5"
    sim_dir.mkdir(parents=True, exist_ok=True)
    out = {}
    for m, df in mats.items():
        p = sim_dir / f"{m}_euclidean_affinity_main5.csv"
        if p.exists():
            a = pd.read_csv(p, index_col=0)
            a.index = a.index.astype(str).map(_normalize_sample_id)
            a.columns = a.columns.astype(str).map(_normalize_sample_id)
            out[m] = a.loc[df.index, df.index]
        else:
            a = euclidean_affinity(df, sigma=1.0)
            a.to_csv(p, encoding="utf-8")
            out[m] = a
    return out


def _compute_weight_components(mats: dict[str, pd.DataFrame], affinity_mats: dict[str, pd.DataFrame], project_root: Path) -> pd.DataFrame:
    clinical = pd.read_csv(project_root / "data" / "raw" / "clinical" / "TCGA.STAD.sampleMap_STAD_clinicalMatrix", sep="\t", dtype=str)
    survival = pd.read_csv(project_root / "data" / "raw" / "clinical" / "survival_STAD_survival.txt", sep="\t", dtype=str)

    n_ref = max(df.shape[0] for df in mats.values())
    rows = []
    for m, df in mats.items():
        x = df.values
        labels = KMeans(n_clusters=2, n_init=20, random_state=42).fit_predict(PCA(n_components=min(30, x.shape[1], x.shape[0] - 1), random_state=42).fit_transform(x))
        q = float(max(0.0, np.nanmean(np.var(x, axis=0))))

        s = pd.Series(labels, index=df.index, name="cluster")
        os_stat = evaluate_clinical_association(s, clinical, survival, endpoint="OS")
        pfi_stat = evaluate_clinical_association(s, clinical, survival, endpoint="PFI")
        pvals = [p for p in [os_stat.get("logrank_p", np.nan), pfi_stat.get("logrank_p", np.nan)] if pd.notna(p)]
        c_score = float(np.mean([-np.log10(max(float(p), 1e-12)) for p in pvals])) if pvals else 0.0

        s_score = float(df.shape[0] / n_ref)
        m_score = float((1.0 - np.mean(np.diag(affinity_mats[m].values))) ** 2)
        rows.append({"modality": m, "Q_m": q, "C_m": c_score, "S_m": s_score, "M_m": m_score})

    comp = pd.DataFrame(rows)
    for col in ["Q_m", "C_m", "S_m", "M_m"]:
        v = comp[col]
        rng = float(v.max() - v.min())
        if rng < 1e-12:
            comp[col] = 0.0
        else:
            comp[col] = (v - v.min()) / rng
    return comp


def _choose_best_label(metrics_path: Path, labels_pattern: str, base_dir: Path) -> tuple[int, pd.Series]:
    mdf = pd.read_csv(metrics_path)
    best_k = int(mdf.sort_values("silhouette", ascending=False).iloc[0]["k"])
    lp = base_dir / labels_pattern.format(k=best_k)
    s = pd.read_csv(lp, index_col=0).iloc[:, 0]
    s.index = s.index.astype(str).map(_normalize_sample_id)
    s.name = "cluster"
    return best_k, s


def _estimate_stability(x: np.ndarray, labels: np.ndarray, k: int) -> tuple[float, float, float]:
    x2 = PCA(n_components=min(30, x.shape[1], x.shape[0] - 1), random_state=42).fit_transform(x)

    def _recluster(z: np.ndarray) -> np.ndarray:
        return KMeans(n_clusters=k, n_init=20, random_state=42).fit_predict(z)

    stat = estimate_stability_by_subsampling(
        x=x2,
        full_labels=labels,
        recluster_fn=_recluster,
        runs=6,
        subsample_ratio=0.8,
        random_state=42,
    )
    return float(stat["mean_nmi"]), float(stat["mean_ari"]), float(stat["consensus_stability"])


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def run_phase23_29() -> None:
    project_root = Path(__file__).resolve().parents[2]
    _, default_cfg = load_project_config(project_root)

    tables_dir = project_root / "results" / "tables"
    figures_dir = project_root / "results" / "figures"
    docs_dir = project_root / "docs"
    logs_dir = project_root / "results" / "logs"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    k_values = list(default_cfg.get("k_range", [2, 3, 4, 5, 6]))

    # Phase 23
    _write_text(
        docs_dir / "23_frontier_improvement_plan.md",
        """# 23 前沿改进映射（frontier-inspired）

## 借鉴并落地的前沿思想
- weighted multi-kernel learning：转译为现有 affinity 的质量加权 + 交互项加权，不引入重型优化器。
- omics-omics interaction-aware integration：在 v2 增加组学对交互矩阵与交互权重。
- weakly paired multi-omics clustering：将 partial fusion 从占位升级为正式重加权实验。
- graph-aware clustering：对谱聚类构图参数做系统网格调优，而非仅记录 warning。
- lightweight shared embedding：采用浅层降维+共享投影原型，作为探索分支。

## 不适合当前项目的方向
- 端到端深图神经网络与大规模对比学习：样本规模与课程周期下复现风险高。
- 高复杂度贝叶斯非参数多视图模型：调参与解释成本高。

## 最终实现的增量升级
- Phase 24: interaction-aware weighted fusion v2
- Phase 25: partial-fusion formal experiments
- Phase 26: consensus ensemble subtyping
- Phase 27: graph connectivity-aware tuning
- Phase 28: shared-embedding prototype（探索）

## 模块映射
- `src/clustering/proposed_weighted_fusion_v2.py`: v2 交互加权融合
- `src/clustering/partial_fusion_experiments.py`: partial-fusion 正式化
- `src/clustering/consensus_ensemble.py`: 共识分型
- `src/feature_engineering/graph_connectivity_tuning.py`: 图连通性调优
- `src/clustering/shared_embedding_prototype.py`: 轻量共享嵌入

## 预期收益与风险
- 预期收益：提高稳定性与簇平衡，在不牺牲可解释性的前提下争取临床分离改进。
- 风险：部分改进可能只提升稳定性不提升生存分离；partial cohort 受上游预处理样本覆盖约束。
""",
    )

    mats = _load_main5_matrices(project_root)
    affinity_mats = _load_or_build_affinities(project_root, mats)
    combined = pd.concat([mats[m] for m in MAIN5_MODALITIES], axis=1)

    # Phase 24
    comp = _compute_weight_components(mats, affinity_mats, project_root)
    weights_v2 = estimate_final_weights(comp, alpha=1.0, beta=1.0, delta=1.0, gamma=1.0)
    v2 = run_weighted_fusion_v2(
        affinity_mats=affinity_mats,
        modality_weights=weights_v2,
        combined_feature_matrix=combined,
        interaction_pairs=[("rna", "methylation"), ("rna", "mirna"), ("cnv", "rna"), ("mutation", "rna")],
        lambda_interaction=0.25,
        k_values=k_values,
        out_tables_dir=tables_dir,
        out_figures_dir=figures_dir,
        prefix="proposed_v2",
        random_state=42,
    )

    # Phase 25
    partial = run_partial_fusion_experiments(
        combined_feature_matrix=combined,
        affinity_mats=affinity_mats,
        modality_weights=weights_v2.set_index("modality")["final_weight"],
        out_tables_dir=tables_dir,
        out_figures_dir=figures_dir,
        k_values=k_values,
        random_state=42,
    )
    _write_text(
        docs_dir / "24_partial_fusion_formalization.md",
        """# 24 Partial Fusion 正式化结果

- 已实现 strict complete-case 与 partial cohort（至少4组学）的正式实验分支。
- partial 分支采用样本级可用组学重归一化权重，不因单组学缺失直接丢样本。
- 本轮在当前处理后的 main5 队列中，partial 主要提升了方法可扩展性，指标增益需以 `partial_fusion_metrics.csv` 为准。
- 若未超过 complete-case baseline，保持如实记录，不更改主结论。
""",
    )

    # Phase 26
    labels_candidates: dict[str, pd.Series] = {}
    base_t = tables_dir
    try:
        bk, s = _choose_best_label(base_t / "baseline_main5_metrics.csv", "baseline_main5_labels_early_fusion_kmeans_k{k}.csv", base_t)
        labels_candidates[f"main5_early_kmeans_k{bk}"] = s
    except Exception:
        pass
    try:
        bk, s = _choose_best_label(base_t / "baseline_main5_metrics.csv", "baseline_main5_labels_early_fusion_spectral_k{k}.csv", base_t)
        labels_candidates[f"main5_early_spectral_k{bk}"] = s
    except Exception:
        pass
    try:
        bk, s = _choose_best_label(base_t / "baseline_main5_metrics.csv", "baseline_main5_labels_equal_weight_fusion_k{k}.csv", base_t)
        labels_candidates[f"main5_equal_weight_k{bk}"] = s
    except Exception:
        pass
    try:
        bk, s = _choose_best_label(base_t / "proposed_main5_metrics.csv", "proposed_main5_labels_k{k}.csv", base_t)
        labels_candidates[f"proposed_v1_k{bk}"] = s
    except Exception:
        pass
    try:
        bk, s = _choose_best_label(base_t / "proposed_v2_metrics.csv", "proposed_v2_labels_k{k}.csv", base_t)
        labels_candidates[f"proposed_v2_k{bk}"] = s
    except Exception:
        pass

    if len(labels_candidates) >= 2:
        consensus = run_consensus_ensemble(
            label_candidates=labels_candidates,
            combined_feature_matrix=combined,
            out_tables_dir=tables_dir,
            out_figures_dir=figures_dir,
            n_clusters=2,
            random_state=42,
        )
    else:
        consensus = {"metrics": pd.DataFrame(), "cluster_sizes": pd.DataFrame(), "labels": pd.Series(dtype=int)}

    _write_text(
        docs_dir / "25_consensus_ensemble_results.md",
        """# 25 Consensus Ensemble 结果

- 候选输入包含 main5 baseline、equal-weight、proposed v1 与 proposed v2（可用即纳入）。
- 采用 co-association matrix + 谱聚类生成共识标签。
- 若 ensemble 稳定性提升但生存分离一般，按证据链如实陈述。
""",
    )

    # Phase 27
    graph = run_graph_connectivity_tuning(
        x_df=combined,
        k_cluster=2,
        out_tables_dir=tables_dir,
        out_figures_dir=figures_dir,
        random_state=42,
    )
    _write_text(
        docs_dir / "26_graph_connectivity_tuning.md",
        """# 26 Graph Connectivity-Aware Tuning

- 对 n_neighbors / mutual_knn / local_scaling / threshold strategy 进行了系统网格搜索。
- 记录了连通分量、孤立点、图密度与聚类质量变化。
- 已输出默认构图与调优后构图的可比结果表，支持答辩时解释 warning 处理策略。
""",
    )

    # Phase 28
    shared = run_shared_embedding_prototype(
        modality_tables=mats,
        k_values=k_values,
        out_tables_dir=tables_dir,
        out_figures_dir=figures_dir,
        random_state=42,
    )
    _write_text(
        docs_dir / "27_shared_embedding_prototype.md",
        """# 27 Shared Embedding Prototype（探索）

- 采用“各组学先降维 + 共享投影”的轻量原型，不使用重型端到端深模型。
- 目标是探索稳定性增益与可解释性折中。
- 若无明显增益，保留为 exploratory，不替代主线结果。
""",
    )

    # Phase 29 unified comparison
    clinical = pd.read_csv(project_root / "data" / "raw" / "clinical" / "TCGA.STAD.sampleMap_STAD_clinicalMatrix", sep="\t", dtype=str)
    survival = pd.read_csv(project_root / "data" / "raw" / "clinical" / "survival_STAD_survival.txt", sep="\t", dtype=str)
    previous_final = pd.read_csv(tables_dir / "final_model_comparison_summary.csv")
    prev_best = previous_final.loc[previous_final["selected_as_final"] == True]
    prev_best_name = prev_best.iloc[0]["method"] if not prev_best.empty else "unknown_previous_final"

    rows = []

    def add_row(method: str, variant: str, cohort: str, k: int, sil: float, labels: pd.Series, notes: str) -> None:
        arr = labels.loc[combined.index].values
        nmi, ari, cons = _estimate_stability(combined.values, arr, k=max(2, k))
        os_stat = evaluate_clinical_association(labels, clinical, survival, endpoint="OS")
        pfi_stat = evaluate_clinical_association(labels, clinical, survival, endpoint="PFI")
        counts = labels.value_counts().values
        balance = float(counts.min() / counts.max()) if counts.size > 0 else 0.0
        rows.append(
            {
                "method": method,
                "variant": variant,
                "cohort": cohort,
                "k": int(k),
                "silhouette": float(sil) if pd.notna(sil) else np.nan,
                "stability_nmi": nmi,
                "stability_ari": ari,
                "consensus_stability": cons,
                "logrank_os_p": float(os_stat.get("logrank_p", np.nan)),
                "logrank_pfi_p": float(pfi_stat.get("logrank_p", np.nan)),
                "cluster_balance_score": balance,
                "selected_as_best_improvement": False,
                "compare_to_previous_final": f"vs {prev_best_name}",
                "notes": notes,
            }
        )

    # proposed v2 best
    v2m = v2["metrics"].sort_values("silhouette", ascending=False).iloc[0]
    v2k = int(v2m["k"])
    add_row("proposed_weighted_fusion", "v2_interaction", "complete_case_main5", v2k, float(v2m["silhouette"]), v2["labels"][v2k], "extra interaction term over v1")

    # partial best per method
    pm = partial["metrics"].sort_values("silhouette", ascending=False)
    for method in pm["method"].unique():
        r = pm[pm["method"] == method].iloc[0]
        k = int(r["k"])
        cohort = str(r["cohort"])
        lbl = partial["labels"][(method, cohort, k)]
        add_row(method, "partial_formal", cohort, k, float(r["silhouette"]), lbl, "sample-wise renormalized partial fusion")

    # consensus
    if not consensus["metrics"].empty:
        cm = consensus["metrics"].iloc[0]
        ck = int(cm["k"])
        add_row("consensus_ensemble", "coassociation", "complete_case_main5", ck, float(cm["silhouette"]), consensus["labels"], "consensus over strong single models")

    # graph tuning best
    gb = graph["best"].iloc[0]
    rows.append(
        {
            "method": "spectral_graph_tuned",
            "variant": "connectivity_aware",
            "cohort": "complete_case_main5",
            "k": 2,
            "silhouette": float(gb["silhouette"]),
            "stability_nmi": np.nan,
            "stability_ari": np.nan,
            "consensus_stability": np.nan,
            "logrank_os_p": np.nan,
            "logrank_pfi_p": np.nan,
            "cluster_balance_score": np.nan,
            "selected_as_best_improvement": False,
            "compare_to_previous_final": f"vs {prev_best_name}",
            "notes": "graph tuning table-level comparison",
        }
    )

    # shared embedding best
    sm = shared["metrics"].sort_values("silhouette", ascending=False).iloc[0]
    sk = int(sm["k"])
    add_row("shared_embedding_prototype", "exploratory", "complete_case_main5", sk, float(sm["silhouette"]), shared["labels"][sk], "lightweight shared projection")

    cmp_df = pd.DataFrame(rows)
    if not cmp_df.empty:
        sil = cmp_df["silhouette"].astype(float).fillna(-1.0)
        cons = cmp_df["consensus_stability"].astype(float).fillna(0.0)
        bal = cmp_df["cluster_balance_score"].astype(float).fillna(0.0)

        # Penalize severely imbalanced/unstable candidates to avoid degenerate picks.
        penalty = ((bal < 0.10).astype(float) * 0.25) + ((cons < 0.10).astype(float) * 0.15)
        cmp_df["improvement_score"] = 0.45 * sil + 0.30 * cons + 0.25 * bal - penalty

        idx = cmp_df["improvement_score"].astype(float).idxmax()
        cmp_df.loc[idx, "selected_as_best_improvement"] = True
    cmp_df.to_csv(tables_dir / "improvement_model_comparison.csv", index=False, encoding="utf-8")

    best = cmp_df.loc[cmp_df["selected_as_best_improvement"] == True].iloc[0]
    previous_sil = np.nan
    if not prev_best.empty and "silhouette" in prev_best.columns:
        previous_sil = float(prev_best.iloc[0]["silhouette"])

    keep_previous = bool(pd.isna(previous_sil) or float(best["silhouette"]) <= previous_sil)
    takeaways = pd.DataFrame(
        [
            {"item": "most_effective_improvement", "value": f"{best['method']}::{best['variant']}"},
            {"item": "main_gain_metrics", "value": "silhouette/stability/balance (see comparison table)"},
            {"item": "non_helpful_improvement", "value": "methods without clear gains remain exploratory"},
            {
                "item": "change_final_model",
                "value": "keep previous final baseline" if keep_previous else "candidate improvement can be considered for final",
            },
        ]
    )
    takeaways.to_csv(tables_dir / "improvement_key_takeaways.csv", index=False, encoding="utf-8")

    _write_text(
        docs_dir / "28_improvement_round_summary.md",
        """# 28 改进轮总结（Phase 23-29）

## 借鉴到的前沿思路
- 多核/多视图加权、组学交互项、弱配对融合、图连通性感知调优、轻量共享嵌入。

## 已实现且可复现的改进
- interaction-aware weighted fusion v2
- partial fusion formal experiments
- consensus ensemble subtyping
- graph connectivity-aware tuning
- shared embedding prototype（exploratory）

## 有效与无效改进
- 以 `improvement_model_comparison.csv` 与 `improvement_key_takeaways.csv` 为准。
- 若改进未超过既有 final baseline，保持原主结论不变。

## 答辩建议表述
- 本轮改进强调“可解释 + 可复现 + 可比较”，并对未增益分支如实披露，避免过度结论。
""",
    )

    # Sync core logs/docs
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    append_lines = (
        f"\n## {now} | Phase 23-29 改进轮完成\n"
        "- 已完成：frontier-inspired mapping、interaction-aware v2、partial formalization、consensus ensemble、graph tuning、shared embedding prototype。\n"
        "- 关键产物：`results/tables/improvement_model_comparison.csv`、`results/tables/improvement_key_takeaways.csv`、`docs/28_improvement_round_summary.md`。\n"
        "- 原 final 主线未覆盖，全部新增结果采用 improvement/v2 独立命名。\n"
    )
    with (docs_dir / "08_progress_log.md").open("a", encoding="utf-8") as f:
        f.write(append_lines)

    with (docs_dir / "07_decision_log.md").open("a", encoding="utf-8") as f:
        f.write(
            "\n## 2026-04-18 决策 21\n"
            "- 决策：在不推翻 final 主线前提下开展 frontier-inspired 增量改进（v2/improvement 命名独立落盘）。\n"
            "- 原因：兼顾可答辩创新与结果可比性，避免覆盖既有 final 证据链。\n"
        )

    with (docs_dir / "06_task_checklist.md").open("a", encoding="utf-8") as f:
        f.write(
            "\n## Phase 23-29 frontier-inspired 改进轮\n"
            "- [已完成] Phase 23 前沿映射到实现计划\n"
            "- [已完成] Phase 24 interaction-aware weighted fusion v2\n"
            "- [已完成] Phase 25 partial-fusion 正式实验\n"
            "- [已完成] Phase 26 consensus ensemble\n"
            "- [已完成] Phase 27 graph connectivity-aware tuning\n"
            "- [已完成] Phase 28 shared embedding prototype（exploratory）\n"
            "- [已完成] Phase 29 改进轮统一汇总\n"
        )

    with (docs_dir / "09_defense_notes.md").open("a", encoding="utf-8") as f:
        f.write(
            "\n## 11. 改进轮答辩要点（Phase 23-29）\n"
            "- 改进轮是增量实验，不覆盖原 final 结论。\n"
            "- 若改进仅提升稳定性、未提升生存分离，也作为正当结果如实报告。\n"
            "- 若最佳模型仍为 baseline，应明确说明‘创新尝试有效但未取代主模型’。\n"
        )

    with (project_root / "README.md").open("a", encoding="utf-8") as f:
        f.write(
            "\n\n第五轮（frontier-inspired improvement）已完成：\n"
            "- Phase 23-29 在不覆盖 final 主线前提下新增 v2/improvement 实验。\n"
            "- 关键汇总见 `results/tables/improvement_model_comparison.csv` 与 `docs/28_improvement_round_summary.md`。\n"
            "- 运行命令：`python -m src.pipelines.run_phase23_29`。\n"
        )

    with (logs_dir / "phase23_29_run.log").open("a", encoding="utf-8") as f:
        f.write(f"[{now}] Phase 23-29 run completed.\n")

    print("Phase 23-29 run completed.")


if __name__ == "__main__":
    run_phase23_29()
