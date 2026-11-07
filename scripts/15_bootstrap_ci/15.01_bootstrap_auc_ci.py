"""Purpose: Bootstrap auc ci for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import os
import numpy as np
import pandas as pd

from sklearn.metrics import (
    roc_auc_score,
    accuracy_score,
    f1_score,
    confusion_matrix,
)

FRAG_CSV = "results/09_ulp_validation/ULP_FRAG_prediction_results_with_scores.csv"
CNA_CSV = "results/09_ulp_validation/ULP_CNA_prediction_results_with_scores.csv"
METH_CSV = "results/09_ulp_validation/ULP_METH_prediction_results_with_scores.csv"
MULTI_CSV = "results/09_ulp_validation/ULP_Multi_prediction_results_with_scores.csv"
GROUP_FILE = "metadata/ULP_groupfile.csv"

OUTPUT_DIR = "results/15_bootstrap_ci"
OUTPUT_CSV = os.path.join(
    OUTPUT_DIR,
    "ULP_model_metrics_AUROC_bootstrap_95CI.csv"
)

N_BOOTSTRAP = 5000
RANDOM_STATE = 2026

MODEL_FILES = {
    "FRAG": {
        "file": FRAG_CSV,
        "prob_col": "FRAG_Score_Prob",
        "label_col": "FRAG_Predicted_Label",
    },
    "CNA": {
        "file": CNA_CSV,
        "prob_col": "CNA_Score_Prob",
        "label_col": "CNA_Predicted_Label",
    },
    "METH": {
        "file": METH_CSV,
        "prob_col": "METH_Score_Prob",
        "label_col": "METH_Predicted_Label",
    },
    "Multi": {
        "file": MULTI_CSV,
        "prob_col": "Multi_Score_Prob",
        "label_col": "Multi_Predicted_Label",
    },
}


def calculate_metrics(y_true, y_prob, y_pred):
    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob, dtype=float)
    y_pred = np.asarray(y_pred).astype(int)

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

    metrics = {}

    if len(np.unique(y_true)) == 2:
        metrics["AUROC"] = roc_auc_score(y_true, y_prob)
    else:
        metrics["AUROC"] = np.nan

    metrics["Accuracy"] = accuracy_score(y_true, y_pred)
    metrics["Sensitivity"] = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    metrics["Specificity"] = tn / (tn + fp) if (tn + fp) > 0 else np.nan
    metrics["Precision"] = tp / (tp + fp) if (tp + fp) > 0 else np.nan
    metrics["F1 score"] = f1_score(y_true, y_pred, zero_division=0)

    return metrics


def bootstrap_ci(y_true, y_prob, y_pred, n_bootstrap=5000, random_state=2026):
    rng = np.random.default_rng(random_state)

    y_true = np.asarray(y_true).astype(int)
    y_prob = np.asarray(y_prob, dtype=float)
    y_pred = np.asarray(y_pred).astype(int)

    n = len(y_true)

    boot_values = {
        "AUROC": [],
        "Accuracy": [],
        "Sensitivity": [],
        "Specificity": [],
        "Precision": [],
        "F1 score": [],
    }

    for _ in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)

        y_t = y_true[idx]
        y_pb = y_prob[idx]
        y_pd = y_pred[idx]

        metrics = calculate_metrics(y_t, y_pb, y_pd)

        for metric, value in metrics.items():
            if not pd.isna(value):
                boot_values[metric].append(value)

    ci = {}

    for metric, values in boot_values.items():
        values = np.asarray(values, dtype=float)

        if len(values) == 0:
            ci[metric] = {
                "CI_lower": np.nan,
                "CI_upper": np.nan,
            }
        else:
            ci[metric] = {
                "CI_lower": np.percentile(values, 2.5),
                "CI_upper": np.percentile(values, 97.5),
            }

    return ci


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    group_df = pd.read_csv(GROUP_FILE)
    group_df["Sample"] = group_df["Sample"].astype(str).str.strip()
    group_df["group"] = group_df["group"].astype(str).str.strip().str.lower()
    group_df["split"] = group_df["split"].astype(str).str.strip().str.lower()
    group_df = group_df[["Sample", "group", "split"]]

    rows = []

    for model_name, info in MODEL_FILES.items():
        df = pd.read_csv(info["file"])

        df["Sample"] = df["Sample"].astype(str).str.strip()
        df = pd.merge(group_df, df, on="Sample", how="inner")

        prob_col = info["prob_col"]
        label_col = info["label_col"]

        df = df[df["group"].isin(["mpe", "bpe"])].copy()
        df["y_true"] = df["group"].map({
            "mpe": 1,
            "bpe": 0,
        })

        df["y_true"] = df["y_true"].astype(int)
        df["y_prob"] = df[prob_col].astype(float)
        df["y_pred"] = df[label_col].astype(int)

        split_order = sorted(df["split"].unique())

        for split_name in split_order:
            sub = df[df["split"] == split_name].copy()

            if sub.empty:
                continue

            y_true = sub["y_true"].values
            y_prob = sub["y_prob"].values
            y_pred = sub["y_pred"].values

            observed_metrics = calculate_metrics(
                y_true,
                y_prob,
                y_pred,
            )

            ci_results = bootstrap_ci(
                y_true,
                y_prob,
                y_pred,
                n_bootstrap=N_BOOTSTRAP,
                random_state=RANDOM_STATE,
            )

            for metric, estimate in observed_metrics.items():
                ci_l = ci_results[metric]["CI_lower"]
                ci_u = ci_results[metric]["CI_upper"]

                rows.append({
                    "Set": split_name,
                    "Model": model_name,
                    "Metric": metric,
                    "Estimate_95CI_percent": (
                        f"{estimate * 100:.1f}% "
                        f"({ci_l * 100:.1f}-{ci_u * 100:.1f}%)"
                    ),
                })

    result_df = pd.DataFrame(rows)
    result_df.to_csv(OUTPUT_CSV, index=False)



if __name__ == "__main__":
    main()
