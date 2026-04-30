"""Internal validation subtype classifier utilities.

Input: subtype labels + selected features
Output: repeated CV metrics for subtype assignment
Purpose: internal validation only (not external test).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.svm import SVC

try:
    from xgboost import XGBClassifier
except Exception:  # pragma: no cover - optional dependency path
    XGBClassifier = None


def _build_model_pipelines(random_state: int = 42) -> dict[str, Pipeline]:
    pre = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            # Leak-safe feature reduction inside CV fold for faster, stabler training.
            ("select", SelectKBest(score_func=f_classif, k=800)),
        ]
    )
    pipelines: dict[str, Pipeline] = {
        "RandomForest": Pipeline(
            steps=[
                ("pre", pre),
                (
                    "model",
                    RandomForestClassifier(
                        n_estimators=300,
                        max_depth=None,
                        min_samples_leaf=2,
                        n_jobs=-1,
                        random_state=random_state,
                    ),
                ),
            ]
        ),
        "SVM": Pipeline(
            steps=[
                ("pre", pre),
                ("model", SVC(C=1.0, kernel="linear", probability=False, random_state=random_state)),
            ]
        ),
    }

    if XGBClassifier is not None:
        pipelines["XGBoost"] = Pipeline(
            steps=[
                ("pre", pre),
                (
                    "model",
                    XGBClassifier(
                        n_estimators=120,
                        max_depth=3,
                        learning_rate=0.05,
                        subsample=0.9,
                        colsample_bytree=0.8,
                        objective="multi:softprob",
                        eval_metric="mlogloss",
                        random_state=random_state,
                    ),
                ),
            ]
        )
    return pipelines


def run_subtype_classifier_internal_validation(
    features: pd.DataFrame,
    labels: pd.Series,
    cv_folds: int = 5,
    cv_repeats: int = 5,
    random_state: int = 42,
) -> dict[str, Any]:
    x = features.copy()
    y = labels.loc[x.index].copy()
    le = LabelEncoder()
    y_enc = le.fit_transform(y.astype(str).values)

    min_class_count = int(pd.Series(y_enc).value_counts().min())
    effective_folds = max(2, min(cv_folds, min_class_count))

    cv = RepeatedStratifiedKFold(n_splits=effective_folds, n_repeats=cv_repeats, random_state=random_state)
    pipelines = _build_model_pipelines(random_state=random_state)

    rows: list[dict[str, Any]] = []
    confusion_rows: list[dict[str, Any]] = []
    importances: list[dict[str, Any]] = []

    for model_name, pipe in pipelines.items():
        accs: list[float] = []
        macro_f1s: list[float] = []
        weighted_f1s: list[float] = []
        aucs: list[float] = []
        for train_idx, test_idx in cv.split(x, y_enc):
            x_train, x_test = x.iloc[train_idx], x.iloc[test_idx]
            y_train, y_test = y_enc[train_idx], y_enc[test_idx]
            if np.unique(y_train).size < 2:
                continue
            try:
                pipe.fit(x_train, y_train)
            except Exception:
                continue
            pred = pipe.predict(x_test)
            accs.append(float(accuracy_score(y_test, pred)))
            macro_f1s.append(float(f1_score(y_test, pred, average="macro")))
            weighted_f1s.append(float(f1_score(y_test, pred, average="weighted")))
            try:
                proba = pipe.predict_proba(x_test)
                aucs.append(float(roc_auc_score(y_test, proba, multi_class="ovr", average="macro")))
            except Exception:
                pass

        if not accs:
            rows.append(
                {
                    "model": model_name,
                    "feature_set": "main5_concat",
                    "cv_scheme": f"RepeatedStratifiedKFold({effective_folds}x{cv_repeats})",
                    "accuracy": np.nan,
                    "macro_f1": np.nan,
                    "weighted_f1": np.nan,
                    "auroc_ovr": np.nan,
                    "notes": "internal_validation_only;model_failed_or_single_class_folds",
                }
            )
            continue

        auroc = float(np.mean(aucs)) if aucs else np.nan

        rows.append(
            {
                "model": model_name,
                "feature_set": "main5_concat",
                "cv_scheme": f"RepeatedStratifiedKFold({effective_folds}x{cv_repeats})",
                "accuracy": float(np.mean(accs)),
                "macro_f1": float(np.mean(macro_f1s)),
                "weighted_f1": float(np.mean(weighted_f1s)),
                "auroc_ovr": auroc,
                "notes": "internal_validation_only",
            }
        )

        # Confusion matrix from a leakage-safe single stratified K-fold partition.
        cm_cv = StratifiedKFold(n_splits=effective_folds, shuffle=True, random_state=random_state)
        try:
            y_pred_once = cross_val_predict(pipe, x, y_enc, cv=cm_cv, method="predict", n_jobs=1)
            cm = confusion_matrix(y_enc, y_pred_once)
        except Exception:
            cm = np.zeros((len(np.unique(y_enc)), len(np.unique(y_enc))), dtype=int)
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                confusion_rows.append(
                    {
                        "model": model_name,
                        "true_class": int(i),
                        "pred_class": int(j),
                        "count": int(cm[i, j]),
                    }
                )

        try:
            pipe.fit(x, y_enc)
            pre = pipe.named_steps["pre"]
            model = pipe.named_steps["model"]
            if hasattr(model, "feature_importances_"):
                select = pre.named_steps["select"]
                selected_cols = x.columns[select.get_support()].tolist()
                feats = pd.Series(model.feature_importances_, index=selected_cols).sort_values(ascending=False)
                for rank, (fname, score) in enumerate(feats.head(200).items(), start=1):
                    importances.append(
                        {
                            "model": model_name,
                            "feature": fname,
                            "importance": float(score),
                            "rank": int(rank),
                        }
                    )
        except Exception:
            pass

    cv_df = pd.DataFrame(rows)
    cm_df = pd.DataFrame(confusion_rows)
    imp_df = pd.DataFrame(importances)
    return {
        "cv_results": cv_df,
        "confusion": cm_df,
        "feature_importance": imp_df,
        "label_encoder": le,
    }


def save_classifier_plots(
    cv_results: pd.DataFrame,
    confusion_df: pd.DataFrame,
    feature_importance_df: pd.DataFrame,
    figures_dir: Path,
) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")

    plt.figure(figsize=(8, 4))
    perf_long = cv_results.melt(
        id_vars=["model"],
        value_vars=["accuracy", "macro_f1", "weighted_f1", "auroc_ovr"],
        var_name="metric",
        value_name="value",
    )
    sns.barplot(data=perf_long, x="metric", y="value", hue="model")
    plt.ylim(0, 1)
    plt.title("Subtype Classifier Internal CV Performance")
    plt.tight_layout()
    plt.savefig(figures_dir / "subtype_classifier_performance.png", dpi=220)
    plt.close()

    if not confusion_df.empty:
        best_model = cv_results.sort_values("macro_f1", ascending=False).iloc[0]["model"]
        sub = confusion_df[confusion_df["model"] == best_model]
        pivot = sub.pivot(index="true_class", columns="pred_class", values="count")
        plt.figure(figsize=(6, 5))
        sns.heatmap(pivot, annot=True, fmt=".0f", cmap="Blues")
        plt.title(f"Confusion Matrix ({best_model})")
        plt.tight_layout()
        plt.savefig(figures_dir / "subtype_classifier_confusion_matrix.png", dpi=220)
        plt.close()

    if not feature_importance_df.empty:
        best_imp_model = (
            cv_results[cv_results["model"].isin(feature_importance_df["model"].unique())]
            .sort_values("macro_f1", ascending=False)
            .iloc[0]["model"]
        )
        top20 = feature_importance_df[feature_importance_df["model"] == best_imp_model].head(20)
        plt.figure(figsize=(8, 6))
        sns.barplot(data=top20, y="feature", x="importance", color="#1d3557")
        plt.title(f"Top20 Feature Importance ({best_imp_model})")
        plt.tight_layout()
        plt.savefig(figures_dir / "subtype_feature_importance_top20.png", dpi=220)
        plt.close()
