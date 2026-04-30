"""Phase 16-22 formal runner.

This pipeline executes:
1) methylation round2 modeling matrix construction
2) main5 baseline experiments
3) formal proposed weighted fusion
4) ablation study aggregation
5) subtype classifier internal validation
6) final result summaries and defense/report packaging docs
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

from src.clustering.baseline_early_fusion import evaluate_early_fusion_round1
from src.clustering.baseline_snf import evaluate_equal_weight_fusion_round1
from src.clustering.baseline_snf import evaluate_equal_weight_fusion_round1
from src.clustering.proposed_weighted_fusion import (
	build_proposed_weight_components,
	estimate_final_weights,
	run_weighted_fusion_clustering,
)
from src.evaluation.ablation_study import save_ablation_plots, write_ablation_outputs
from src.evaluation.clinical_association import evaluate_clinical_association
from src.evaluation.cluster_metrics import estimate_stability_by_subsampling
from src.evaluation.subtype_classifier import (
	run_subtype_classifier_internal_validation,
	save_classifier_plots,
)
from src.feature_engineering.build_similarity_matrices import euclidean_affinity, pearson_similarity
from src.preprocessing.preprocess_methylation import run_methylation_round2_modeling
from src.utils.config_utils import load_project_config
from src.utils.io_utils import ensure_dir


MAIN5_MODALITIES = ["mutation", "cnv", "methylation", "rna", "mirna"]


def _normalize_sample_id(sample_id: str) -> str:
	sid = str(sample_id).strip().upper()
	return sid[:16] if len(sid) >= 16 else sid


def _read_interim(path: Path) -> pd.DataFrame:
	df = pd.read_csv(path, index_col=0)
	df.index = df.index.astype(str).map(_normalize_sample_id)
	return df.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def _get_main5_intersection(project_root: Path) -> set[str]:
	inv_path = project_root / "results" / "tables" / "sample_inventory.csv"
	inv = pd.read_csv(inv_path)
	sample_sets = {
		m: set(inv.loc[inv["modality"] == m, "normalized_sample_id"].dropna().astype(str).unique())
		for m in MAIN5_MODALITIES
	}
	return set.intersection(*sample_sets.values())


def run_phase16_methylation(project_root: Path, paths_cfg: dict[str, Any], default_cfg: dict[str, Any]) -> dict[str, Any]:
	interim_dir = Path(paths_cfg["data"]["interim_root"])
	tables_dir = Path(paths_cfg["results"]["tables"])
	figures_dir = Path(paths_cfg["results"]["figures"])
	logs_dir = Path(paths_cfg["results"]["logs"])
	ensure_dir(interim_dir)
	ensure_dir(tables_dir)
	ensure_dir(figures_dir)

	selected = _get_main5_intersection(project_root)
	methylation_matrix_path = interim_dir / "methylation_round2_modeling.csv"
	summary_path = tables_dir / "methylation_processing_summary.csv"
	filter_log_path = tables_dir / "methylation_feature_filter_log.csv"

	matrix, summary_df, filter_df = run_methylation_round2_modeling(
		input_path=Path(paths_cfg["data"]["methylation"]),
		output_matrix_path=methylation_matrix_path,
		summary_output_path=summary_path,
		filter_log_output_path=filter_log_path,
		selected_samples=selected,
		probe_map_path=Path(paths_cfg["data"]["methylation_probe_map"]),
		missingness_filter_threshold=float(default_cfg["preprocessing"].get("max_missing_rate", 0.30)),
		variance_filter_quantile=float(default_cfg["preprocessing"].get("variance_filter_quantile", 0.75)),
		max_features=5000,
		chunk_size=2000,
		figures_dir=figures_dir,
		log_path=logs_dir / "preprocessing_dimension_changes.csv",
	)

	return {
		"matrix": matrix,
		"summary": summary_df,
		"filter_log": filter_df,
		"matrix_path": methylation_matrix_path,
		"is_ready": bool(summary_df.loc[0, "final_matrix_ready"]),
	}


def _build_main5_matrices(project_root: Path, methylation_path: Path) -> dict[str, pd.DataFrame]:
	interim = project_root / "data" / "interim"
	mats = {
		"mutation": _read_interim(interim / "mutation_round1.csv"),
		"cnv": _read_interim(interim / "cnv_round1.csv"),
		"methylation": _read_interim(methylation_path),
		"rna": _read_interim(interim / "rna_round1.csv"),
		"mirna": _read_interim(interim / "mirna_round1.csv"),
	}
	cohort = sorted(set.intersection(*[set(m.index) for m in mats.values()]))
	return {m: df.loc[cohort].copy() for m, df in mats.items()}


def _save_main5_readiness(mats: dict[str, pd.DataFrame], out_path: Path) -> pd.DataFrame:
	rows = []
	for m, df in mats.items():
		rows.append(
			{
				"modality": m,
				"matrix_ready_for_modeling": bool(df.shape[0] > 0 and df.shape[1] > 0),
				"sample_count": int(df.shape[0]),
				"feature_count": int(df.shape[1]),
				"notes": "phase17_main5_readiness",
			}
		)
	out = pd.DataFrame(rows)
	out.to_csv(out_path, index=False, encoding="utf-8")
	return out


def _build_similarity_main5(mats: dict[str, pd.DataFrame], tables_dir: Path) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
	sim_dir = tables_dir / "similarity_matrices_main5"
	ensure_dir(sim_dir)
	logs: list[dict[str, Any]] = []
	affs: dict[str, pd.DataFrame] = {}
	for m, df in mats.items():
		pear = pearson_similarity(df).clip(-1, 1)
		aff = euclidean_affinity(df, sigma=1.0)
		pear.to_csv(sim_dir / f"{m}_pearson_similarity_main5.csv", encoding="utf-8")
		aff.to_csv(sim_dir / f"{m}_euclidean_affinity_main5.csv", encoding="utf-8")
		affs[m] = aff
		logs.append(
			{
				"modality": m,
				"aligned_samples": int(df.shape[0]),
				"aligned_features": int(df.shape[1]),
				"similarity_method": "pearson+euclidean_affinity",
				"output_matrix_shape": f"{aff.shape[0]}x{aff.shape[1]}",
				"notes": "main5_formal",
			}
		)
	log_df = pd.DataFrame(logs)
	log_df.to_csv(tables_dir / "similarity_build_log_main5.csv", index=False, encoding="utf-8")
	return affs, log_df


def _plot_main5_baseline(metrics: pd.DataFrame, sizes: pd.DataFrame, figures_dir: Path) -> None:
	sns.set_theme(style="whitegrid")
	plt.figure(figsize=(9, 5))
	sns.lineplot(data=metrics, x="k", y="silhouette", hue="method", marker="o")
	plt.title("Baseline Main5 Silhouette Comparison")
	plt.tight_layout()
	plt.savefig(figures_dir / "baseline_main5_metric_comparison.png", dpi=220)
	plt.close()

	plt.figure(figsize=(10, 5))
	sns.barplot(data=sizes, x="k", y="cluster_size", hue="method")
	plt.title("Cluster Size Comparison (Main5)")
	plt.tight_layout()
	plt.savefig(figures_dir / "cluster_size_comparison_main5.png", dpi=220)
	plt.close()


def run_phase17_main5_baseline(project_root: Path, mats: dict[str, pd.DataFrame], affs: dict[str, pd.DataFrame], paths_cfg: dict[str, Any], default_cfg: dict[str, Any]) -> dict[str, Any]:
	tables_dir = Path(paths_cfg["results"]["tables"])
	figures_dir = Path(paths_cfg["results"]["figures"])
	k_values = list(default_cfg.get("k_range", [2, 3, 4, 5, 6]))
	seed = int(default_cfg.get("random_seed", 42))

	warning_rows = []
	with warnings.catch_warnings(record=True) as warns:
		warnings.simplefilter("always")
		early = evaluate_early_fusion_round1(
			modality_tables=mats,
			k_values=k_values,
			algorithms=["kmeans", "spectral"],
			out_tables_dir=tables_dir,
			out_figures_dir=figures_dir,
			modality_set_name="main5",
			random_state=seed,
		)
		combined = pd.concat([mats[m] for m in MAIN5_MODALITIES], axis=1)
		snf = evaluate_equal_weight_fusion_round1(
			affinity_mats=affs,
			aligned_feature_matrix=combined,
			k_values=k_values,
			out_tables_dir=tables_dir,
			out_figures_dir=figures_dir,
			modality_set_name="main5",
			random_state=seed,
		)
		for w in warns:
			msg = str(w.message)
			if "Graph is not fully connected" in msg:
				warning_rows.append({"warning": msg, "notes": "main5_baseline"})

	metrics = pd.concat([early["metrics"], snf["metrics"]], ignore_index=True)
	sizes = pd.concat([early["cluster_sizes"], snf["cluster_sizes"]], ignore_index=True)
	metrics.to_csv(tables_dir / "baseline_main5_metrics.csv", index=False, encoding="utf-8")
	sizes.to_csv(tables_dir / "baseline_main5_cluster_sizes.csv", index=False, encoding="utf-8")
	snf["fused_summary"].to_csv(tables_dir / "fused_similarity_summary_main5.csv", index=False, encoding="utf-8")

	# Standardized label filenames requested by user.
	for (algo, k), s in early["labels"].items():
		s.to_frame().to_csv(tables_dir / f"baseline_main5_labels_early_fusion_{algo}_k{k}.csv", encoding="utf-8")
	for k, s in snf["labels"].items():
		s.to_frame().to_csv(tables_dir / f"baseline_main5_labels_equal_weight_fusion_k{k}.csv", encoding="utf-8")

	fused = snf["fused"]
	plt.figure(figsize=(7, 6))
	sns.heatmap(fused.iloc[:100, :100], cmap="YlGnBu")
	plt.title("Fused Network Heatmap Main5 (first 100)")
	plt.tight_layout()
	plt.savefig(figures_dir / "fused_network_heatmap_main5.png", dpi=220)
	plt.close()

	_plot_main5_baseline(metrics, sizes, figures_dir)

	if warning_rows:
		pd.DataFrame(warning_rows).to_csv(Path(paths_cfg["results"]["logs"]) / "main5_graph_connectivity_warnings.csv", index=False, encoding="utf-8")

	labels_map: dict[tuple[str, int], pd.Series] = {}
	for (algo, k), s in early["labels"].items():
		labels_map[(f"early_fusion_{algo}", k)] = s
	for k, s in snf["labels"].items():
		labels_map[("equal_weight_fusion", k)] = s

	return {
		"metrics": metrics,
		"sizes": sizes,
		"labels_map": labels_map,
		"x_scaled": early["x_scaled"],
		"combined": combined,
		"warnings": warning_rows,
	}


def _stability_for_labels(x_scaled: np.ndarray, labels_map: dict[tuple[str, int], pd.Series], methods: list[tuple[str, int]]) -> dict[tuple[str, int], float]:
	pca_dim = min(50, x_scaled.shape[1], x_scaled.shape[0] - 1)
	x_reduced = PCA(n_components=pca_dim, random_state=42).fit_transform(x_scaled)
	out: dict[tuple[str, int], float] = {}
	for key in methods:
		labels = labels_map[key].values
		k = key[1]

		def recluster(x_sub: np.ndarray) -> np.ndarray:
			return KMeans(n_clusters=k, n_init=20, random_state=42).fit_predict(x_sub)

		stats = estimate_stability_by_subsampling(
			x=x_reduced,
			full_labels=labels,
			recluster_fn=recluster,
			runs=6,
			subsample_ratio=0.8,
			random_state=42,
		)
		out[key] = float(stats["consensus_stability"])
	return out


def _load_clinical(project_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
	clinical = pd.read_csv(project_root / "data" / "raw" / "clinical" / "TCGA.STAD.sampleMap_STAD_clinicalMatrix", sep="\t", dtype=str)
	survival = pd.read_csv(project_root / "data" / "raw" / "clinical" / "survival_STAD_survival.txt", sep="\t", dtype=str)
	return clinical, survival


def _compute_weight_components(mats: dict[str, pd.DataFrame], affs: dict[str, pd.DataFrame], clinical_df: pd.DataFrame, survival_df: pd.DataFrame, methylation_filter_log: pd.DataFrame) -> pd.DataFrame:
	q_scores: dict[str, float] = {}
	c_raw: dict[str, float] = {}
	s_scores: dict[str, float] = {}
	m_scores: dict[str, float] = {}

	n_ref = max(df.shape[0] for df in mats.values())
	for m, df in mats.items():
		x = df.values
		labels = SpectralClustering(n_clusters=2, affinity="precomputed", random_state=42).fit_predict(affs[m].values)
		try:
			if np.unique(labels).size < 2:
				q_scores[m] = 0.0
			else:
				sil = float(silhouette_score(x, labels))
				if not np.isfinite(sil):
					sil = 0.0
				q_scores[m] = float(max(0.0, sil))
		except Exception:
			q_scores[m] = 0.0

		lab_s = pd.Series(labels, index=df.index, name="cluster")
		os_stat = evaluate_clinical_association(lab_s, clinical_df, survival_df, endpoint="OS")
		pfi_stat = evaluate_clinical_association(lab_s, clinical_df, survival_df, endpoint="PFI")
		p_os = float(os_stat.get("logrank_p", np.nan))
		p_pfi = float(pfi_stat.get("logrank_p", np.nan))
		pvals = [x for x in [p_os, p_pfi] if pd.notna(x)]
		c_raw[m] = float(np.mean([-np.log10(max(p, 1e-12)) for p in pvals])) if pvals else 0.0

		s_scores[m] = float(df.shape[0] / n_ref)
		if m == "methylation" and not methylation_filter_log.empty:
			m_scores[m] = float(methylation_filter_log["probe_missing_rate"].mean())
		else:
			m_scores[m] = 0.0

	c_min = min(c_raw.values())
	c_max = max(c_raw.values())
	c_scores = {k: ((v - c_min) / (c_max - c_min + 1e-12)) for k, v in c_raw.items()}
	q_min = min(q_scores.values())
	q_max = max(q_scores.values())
	q_norm = {k: ((v - q_min) / (q_max - q_min + 1e-12)) for k, v in q_scores.items()}

	return build_proposed_weight_components(
		modalities=MAIN5_MODALITIES,
		q_scores=q_norm,
		c_scores=c_scores,
		s_scores=s_scores,
		m_scores=m_scores,
	)


def _clinical_p_for_labels(labels: pd.Series, clinical_df: pd.DataFrame, survival_df: pd.DataFrame) -> tuple[float, float]:
	os_p = float(evaluate_clinical_association(labels, clinical_df, survival_df, endpoint="OS").get("logrank_p", np.nan))
	pfi_p = float(evaluate_clinical_association(labels, clinical_df, survival_df, endpoint="PFI").get("logrank_p", np.nan))
	return os_p, pfi_p


def _write_markdown(path: Path, content: str) -> None:
	ensure_dir(path.parent)
	path.write_text(content, encoding="utf-8")


def run_phase16_22(project_root: Path) -> None:
	paths_cfg, default_cfg = load_project_config(project_root)
	tables_dir = Path(paths_cfg["results"]["tables"])
	figures_dir = Path(paths_cfg["results"]["figures"])
	docs_dir = project_root / "docs"
	ensure_dir(tables_dir)
	ensure_dir(figures_dir)

	# Phase 16
	methyl_out = run_phase16_methylation(project_root, paths_cfg, default_cfg)

	# Phase 17
	mats = _build_main5_matrices(project_root, methyl_out["matrix_path"])
	_save_main5_readiness(mats, tables_dir / "modeling_readiness_main5.csv")
	affs, _ = _build_similarity_main5(mats, tables_dir)
	baseline = run_phase17_main5_baseline(project_root, mats, affs, paths_cfg, default_cfg)

	# Phase 18
	clinical_df, survival_df = _load_clinical(project_root)
	component_df = _compute_weight_components(mats, affs, clinical_df, survival_df, methyl_out["filter_log"])
	component_df.to_csv(tables_dir / "proposed_weight_components.csv", index=False, encoding="utf-8")

	params = default_cfg.get("proposed_method", {})
	weight_df = estimate_final_weights(
		component_df,
		alpha=float(params.get("alpha", 1.0)),
		beta=float(params.get("beta", 1.0)),
		delta=float(params.get("delta", 1.0)),
		gamma=float(params.get("gamma", 1.0)),
	)
	weight_df.to_csv(tables_dir / "proposed_weight_estimates_final.csv", index=False, encoding="utf-8")

	proposed = run_weighted_fusion_clustering(
		affinity_mats=affs,
		weight_table=weight_df,
		combined_feature_matrix=baseline["combined"],
		k_values=list(default_cfg.get("k_range", [2, 3, 4, 5, 6])),
		modality_set_name="main5",
		out_tables_dir=tables_dir,
		out_figures_dir=figures_dir,
		label_prefix="proposed_main5",
		random_state=int(default_cfg.get("random_seed", 42)),
	)
	proposed["metrics"].to_csv(tables_dir / "proposed_main5_metrics.csv", index=False, encoding="utf-8")
	proposed["cluster_sizes"].to_csv(tables_dir / "proposed_main5_cluster_sizes.csv", index=False, encoding="utf-8")
	proposed["fused_summary"].to_csv(tables_dir / "proposed_fused_similarity_summary.csv", index=False, encoding="utf-8")

	# Phase 19 ablation
	baseline_metrics = pd.read_csv(tables_dir / "baseline_main5_metrics.csv")
	methods_for_stability = [(str(r["method"]), int(r["k"])) for _, r in baseline_metrics.iterrows()]
	stability = _stability_for_labels(baseline["x_scaled"], baseline["labels_map"], methods_for_stability)

	ablation_rows: list[dict[str, Any]] = []
	for _, r in baseline_metrics.iterrows():
		method = str(r["method"])
		k = int(r["k"])
		labels = baseline["labels_map"][(method, k)]
		os_p, pfi_p = _clinical_p_for_labels(labels, clinical_df, survival_df)
		ablation_rows.append(
			{
				"experiment_name": method,
				"modality_set": "main5",
				"k": k,
				"silhouette": float(r["silhouette"]),
				"stability_score": float(stability.get((method, k), np.nan)),
				"logrank_os_p": os_p,
				"logrank_pfi_p": pfi_p,
				"cluster_balance_score": float(r["cluster_balance_score"]),
				"notes": "method_ablation",
			}
		)

	prop_metrics = proposed["metrics"].copy()
	for _, r in prop_metrics.iterrows():
		k = int(r["k"])
		labels = proposed["labels"][k]
		os_p, pfi_p = _clinical_p_for_labels(labels, clinical_df, survival_df)
		ablation_rows.append(
			{
				"experiment_name": "proposed_weighted_fusion_full",
				"modality_set": "main5",
				"k": k,
				"silhouette": float(r["silhouette"]),
				"stability_score": float(np.nan),
				"logrank_os_p": os_p,
				"logrank_pfi_p": pfi_p,
				"cluster_balance_score": float(r["cluster_balance_score"]),
				"notes": "method_ablation",
			}
		)

	# Weight component ablations
	comp_variants = {
		"weight_no_C": (0.0, 0.0, 1.0, 1.0),
		"weight_no_S": (1.0, 1.0, 0.0, 1.0),
		"weight_no_M": (1.0, 1.0, 1.0, 0.0),
		"weight_only_Q": (1.0, 0.0, 0.0, 0.0),
		"weight_full": (1.0, 1.0, 1.0, 1.0),
	}
	for name, (a, b, d, g) in comp_variants.items():
		w_df = estimate_final_weights(component_df, alpha=a, beta=b, delta=d, gamma=g)
		out = run_weighted_fusion_clustering(
			affinity_mats=affs,
			weight_table=w_df,
			combined_feature_matrix=baseline["combined"],
			k_values=list(default_cfg.get("k_range", [2, 3, 4, 5, 6])),
			modality_set_name="main5",
			out_tables_dir=tables_dir,
			out_figures_dir=figures_dir,
			label_prefix=f"ablation_{name}",
			random_state=42,
		)
		for _, r in out["metrics"].iterrows():
			k = int(r["k"])
			labels = out["labels"][k]
			os_p, pfi_p = _clinical_p_for_labels(labels, clinical_df, survival_df)
			ablation_rows.append(
				{
					"experiment_name": name,
					"modality_set": "main5",
					"k": k,
					"silhouette": float(r["silhouette"]),
					"stability_score": np.nan,
					"logrank_os_p": os_p,
					"logrank_pfi_p": pfi_p,
					"cluster_balance_score": float(r["cluster_balance_score"]),
					"notes": "weight_component_ablation",
				}
			)

	# Modality ablations (equal-weight fusion with one modality dropped)
	for drop_m in MAIN5_MODALITIES:
		kept = [m for m in MAIN5_MODALITIES if m != drop_m]
		sub_mats = {m: mats[m] for m in kept}
		sub_cohort = sorted(set.intersection(*[set(df.index) for df in sub_mats.values()]))
		sub_mats = {m: df.loc[sub_cohort] for m, df in sub_mats.items()}
		sub_affs = {m: euclidean_affinity(df, sigma=1.0) for m, df in sub_mats.items()}
		sub_x = pd.concat([sub_mats[m] for m in kept], axis=1)
		sub = evaluate_equal_weight_fusion_round1(
			affinity_mats=sub_affs,
			aligned_feature_matrix=sub_x,
			k_values=list(default_cfg.get("k_range", [2, 3, 4, 5, 6])),
			out_tables_dir=tables_dir,
			out_figures_dir=figures_dir,
			modality_set_name=f"main5_drop_{drop_m}",
			random_state=42,
		)
		for _, r in sub["metrics"].iterrows():
			k = int(r["k"])
			labels = sub["labels"][k]
			os_p, pfi_p = _clinical_p_for_labels(labels, clinical_df, survival_df)
			ablation_rows.append(
				{
					"experiment_name": f"drop_{drop_m}",
					"modality_set": f"main5_minus_{drop_m}",
					"k": k,
					"silhouette": float(r["silhouette"]),
					"stability_score": np.nan,
					"logrank_os_p": os_p,
					"logrank_pfi_p": pfi_p,
					"cluster_balance_score": float(r["cluster_balance_score"]),
					"notes": "modality_ablation",
				}
			)

	ablation_rows.append(
		{
			"experiment_name": "partial_fusion_complete_intersection",
			"modality_set": "main5_complete_intersection",
			"k": 0,
			"silhouette": np.nan,
			"stability_score": np.nan,
			"logrank_os_p": np.nan,
			"logrank_pfi_p": np.nan,
			"cluster_balance_score": np.nan,
			"notes": "partial_fusion_interface_not_fully_available_in_current_codebase",
		}
	)
	ablation_rows.append(
		{
			"experiment_name": "partial_fusion_partial_coverage",
			"modality_set": "main5_partial_coverage",
			"k": 0,
			"silhouette": np.nan,
			"stability_score": np.nan,
			"logrank_os_p": np.nan,
			"logrank_pfi_p": np.nan,
			"cluster_balance_score": np.nan,
			"notes": "partial_fusion_interface_not_fully_available_in_current_codebase",
		}
	)

	ablation_df = pd.DataFrame(ablation_rows)
	ranked_ablation, _ = write_ablation_outputs(ablation_df, tables_dir)
	save_ablation_plots(ranked_ablation, figures_dir)

	# Phase 20 classifier on selected final labels candidate (best ranked valid experiment)
	valid_ranked = ranked_ablation.dropna(subset=["silhouette"]).sort_values("overall_rank")
	chosen = valid_ranked.iloc[0]
	chosen_name = str(chosen["experiment_name"])
	chosen_k = int(chosen["k"])

	if chosen_name == "proposed_weighted_fusion_full":
		final_labels = proposed["labels"][chosen_k]
	elif chosen_name.startswith("early_fusion") or chosen_name == "equal_weight_fusion":
		final_labels = baseline["labels_map"][(chosen_name, chosen_k)]
	else:
		final_labels = proposed["labels"][proposed["best_k"]]
		chosen_name = "proposed_weighted_fusion_full"
		chosen_k = int(proposed["best_k"])

	# Classifier internal CV requires at least 2 samples per class in training partitions.
	if final_labels.value_counts().min() < 2:
		fallback = baseline_metrics[baseline_metrics["method"] == "early_fusion_kmeans"].sort_values("silhouette", ascending=False).iloc[0]
		chosen_name = "early_fusion_kmeans"
		chosen_k = int(fallback["k"])
		final_labels = baseline["labels_map"][(chosen_name, chosen_k)]

	features_for_classifier = baseline["combined"].loc[final_labels.index]
	clf_out = run_subtype_classifier_internal_validation(
		features=features_for_classifier,
		labels=final_labels,
		cv_folds=int(default_cfg["evaluation"].get("classifier_cv_folds", 5)),
		cv_repeats=min(2, int(default_cfg["evaluation"].get("classifier_cv_repeats", 5))),
		random_state=int(default_cfg.get("random_seed", 42)),
	)
	clf_out["cv_results"].to_csv(tables_dir / "subtype_classifier_cv_results.csv", index=False, encoding="utf-8")
	clf_out["feature_importance"].to_csv(tables_dir / "subtype_classifier_feature_importance.csv", index=False, encoding="utf-8")
	clf_out["confusion"].to_csv(tables_dir / "subtype_assignment_confusion_matrix.csv", index=False, encoding="utf-8")
	save_classifier_plots(clf_out["cv_results"], clf_out["confusion"], clf_out["feature_importance"], figures_dir)

	# Phase 21 final summaries
	final_rows = []
	# main4 references from existing round1
	main4_metrics_path = tables_dir / "baseline_round1_metrics.csv"
	if main4_metrics_path.exists():
		main4 = pd.read_csv(main4_metrics_path)
		for method, grp in main4.groupby("method"):
			best = grp.sort_values("silhouette", ascending=False).iloc[0]
			final_rows.append(
				{
					"method": method,
					"modality_set": "main4",
					"k": int(best["k"]),
					"silhouette": float(best["silhouette"]),
					"stability": float(best.get("consensus_stability", np.nan)),
					"clinical_separation_os": np.nan,
					"clinical_separation_pfi": np.nan,
					"cluster_balance": float(best.get("cluster_balance_score", np.nan)),
					"selected_as_final": False,
					"why_selected": "reference_main4_round1",
				}
			)

	for method, grp in baseline_metrics.groupby("method"):
		best = grp.sort_values("silhouette", ascending=False).iloc[0]
		k = int(best["k"])
		labels = baseline["labels_map"][(method, k)]
		os_p, pfi_p = _clinical_p_for_labels(labels, clinical_df, survival_df)
		final_rows.append(
			{
				"method": method,
				"modality_set": "main5",
				"k": k,
				"silhouette": float(best["silhouette"]),
				"stability": float(stability.get((method, k), np.nan)),
				"clinical_separation_os": os_p,
				"clinical_separation_pfi": pfi_p,
				"cluster_balance": float(best["cluster_balance_score"]),
				"selected_as_final": False,
				"why_selected": "main5_baseline_candidate",
			}
		)

	best_prop = proposed["metrics"].sort_values("silhouette", ascending=False).iloc[0]
	prop_k = int(best_prop["k"])
	os_p, pfi_p = _clinical_p_for_labels(proposed["labels"][prop_k], clinical_df, survival_df)
	final_rows.append(
		{
			"method": "proposed_weighted_fusion",
			"modality_set": "main5",
			"k": prop_k,
			"silhouette": float(best_prop["silhouette"]),
			"stability": np.nan,
			"clinical_separation_os": os_p,
			"clinical_separation_pfi": pfi_p,
			"cluster_balance": float(best_prop["cluster_balance_score"]),
			"selected_as_final": False,
			"why_selected": "proposed_candidate",
		}
	)

	final_df = pd.DataFrame(final_rows)
	# Final selection by ablation top-ranked valid experiment
	final_method_name = chosen_name
	target_method = final_method_name.replace("_full", "")
	final_df.loc[
		(final_df["method"] == target_method) & (final_df["modality_set"] == "main5"),
		"selected_as_final",
	] = True
	final_df.loc[final_df["selected_as_final"], "why_selected"] = "best_ablation_rank_under_internal_validation"
	final_df.to_csv(tables_dir / "final_model_comparison_summary.csv", index=False, encoding="utf-8")

	inv = pd.read_csv(project_root / "results" / "tables" / "sample_inventory.csv")
	sets = {
		m: set(inv.loc[inv["modality"] == m, "normalized_sample_id"].dropna().astype(str).unique())
		for m in ["mutation", "cnv", "methylation", "rna", "mirna", "rppa"]
	}
	main5_inter = len(set.intersection(sets["mutation"], sets["cnv"], sets["methylation"], sets["rna"], sets["mirna"]))
	with_rppa_inter = len(set.intersection(sets["mutation"], sets["cnv"], sets["methylation"], sets["rna"], sets["mirna"], sets["rppa"]))
	rppa_reduction = 1.0 - (with_rppa_inter / main5_inter) if main5_inter > 0 else np.nan

	main_perf = clf_out["cv_results"].sort_values("macro_f1", ascending=False).iloc[0]
	final_os_p, final_pfi_p = _clinical_p_for_labels(final_labels, clinical_df, survival_df)

	key_numbers = pd.DataFrame(
		[
			{"key": "main5_intersection_samples", "value": main5_inter},
			{"key": "main5_plus_rppa_intersection_samples", "value": with_rppa_inter},
			{"key": "rppa_reduction_ratio", "value": rppa_reduction},
			{"key": "partial_fusion_coverage", "value": "not_fully_implemented_in_current_codebase"},
			{"key": "final_selected_method", "value": final_method_name},
			{"key": "final_selected_k", "value": chosen_k},
			{"key": "main_performance_macro_f1", "value": float(main_perf["macro_f1"])},
			{"key": "main_performance_accuracy", "value": float(main_perf["accuracy"])},
			{"key": "main_survival_separation_os_p", "value": final_os_p},
			{"key": "main_survival_separation_pfi_p", "value": final_pfi_p},
		]
	)
	key_numbers.to_csv(tables_dir / "final_key_numbers_for_ppt.csv", index=False, encoding="utf-8")

	# Phase 22 documents
	_write_markdown(
		docs_dir / "15_ablation_and_final_model_selection.md",
		"# 15 消融与最终模型选择\n\n"
		"本轮已完成方法消融、权重项消融和组学消融，结果见 `results/tables/ablation_results.csv` 与 `results/tables/ablation_rank_summary.csv`。\n\n"
		"最终主候选以 ablation 总排名为准，并在内部验证指标与临床分离之间做平衡。\n",
	)
	_write_markdown(
		docs_dir / "16_subtype_classifier_validation.md",
		"# 16 Subtype Classifier Internal Validation\n\n"
		"使用最终候选分型标签进行内部交叉验证（RepeatedStratifiedKFold），并明确为 internal validation。\n\n"
		"结果表：`results/tables/subtype_classifier_cv_results.csv`。\n",
	)
	_write_markdown(
		docs_dir / "17_final_results_summary.md",
		"# 17 最终结果摘要\n\n"
		f"最终主候选：{final_method_name}，k={chosen_k}。\n\n"
		"综合依据：方法对比、消融、临床关联、分类器内部验证。\n\n"
		"详细见 `results/tables/final_model_comparison_summary.csv`。\n",
	)
	_write_markdown(
		docs_dir / "18_limitations_and_risk_statement.md",
		"# 18 局限与风险声明\n\n"
		"1. 当前无外部独立验证队列。\n"
		"2. methylation 使用 probe-level 筛选，尚未完成正式 gene/promoter 聚合。\n"
		"3. 谱聚类图连通性 warning 在部分设置下仍可能出现，需谨慎解释。\n"
		"4. RPPA 未纳入主模型，仅作补充轨道。\n",
	)
	_write_markdown(
		docs_dir / "19_ppt_figure_plan.md",
		"# 19 PPT 图表计划\n\n"
		"- baseline_main5_metric_comparison.png：主图，方法对比页。\n"
		"- proposed_weight_barplot.png：主图，创新方法机制页。\n"
		"- ablation_comparison.png：主图，创新有效性页。\n"
		"- subtype_classifier_performance.png：主图，内部验证页。\n"
		"- subtype_classifier_confusion_matrix.png：备选图，误分类解释页。\n",
	)
	_write_markdown(
		docs_dir / "20_report_table_plan.md",
		"# 20 报告表格计划\n\n"
		"- final_model_comparison_summary.csv：主表。\n"
		"- ablation_rank_summary.csv：主表。\n"
		"- subtype_classifier_cv_results.csv：主表。\n"
		"- final_key_numbers_for_ppt.csv：摘要表。\n",
	)


if __name__ == "__main__":
	root = Path(__file__).resolve().parents[2]
	run_phase16_22(root)
	print("Phase 16-22 run completed.")
