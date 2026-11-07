"""Purpose: Reproduce the Cluster 2 ROC curves and Table S7 performance metrics.

Inputs: Public cluster assignments, Discovery/ULP model predictions, and optional
        approved cytology metadata containing Sample, CB, and LBC.
Outputs: Discovery, ULP, and combined ROC PDFs; AUROC and bootstrap metric tables.
Usage: python scripts/17_cluster2_analysis/17.01_cluster2_roc_and_metrics.py
"""

import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, roc_auc_score, roc_curve


CLUSTER_FILE = "results/17_cluster2_analysis/cluster_assignments.csv"
DISCOVERY_METADATA_FILE = "metadata/discovery_groupfile.csv"
ULP_METADATA_FILE = "metadata/ULP_groupfile.csv"
OUTPUT_DIR = "results/17_cluster2_analysis"
OUTPUT_AUC = os.path.join(OUTPUT_DIR, "cluster2_ROC_AUC.csv")
OUTPUT_METRICS = os.path.join(OUTPUT_DIR, "Table_S7_cluster2_metrics_recalculated.csv")

N_BOOTSTRAP = 5000
RANDOM_STATE = 2026
CYTOLOGY_POSITIVE_VALUE = 3

PREDICTION_FILES = {
    "Discovery": {
        "METH": "results/02_nanopore_methylation/discovery_METH_prediction_results_with_scores.csv",
        "CNA": "results/03_nanopore_cna/discovery_CNA_prediction_results_with_scores.csv",
        "FRAG": "results/04_nanopore_fragmentomics/discovery_FRAG_prediction_results_with_scores.csv",
        "Multi": "results/05_nanopore_multimodal/discovery_Multi_prediction_results_with_scores.csv",
    },
    "ULP": {
        "METH": "results/09_ulp_validation/ULP_METH_prediction_results_with_scores.csv",
        "CNA": "results/09_ulp_validation/ULP_CNA_prediction_results_with_scores.csv",
        "FRAG": "results/09_ulp_validation/ULP_FRAG_prediction_results_with_scores.csv",
        "Multi": "results/09_ulp_validation/ULP_Multi_prediction_results_with_scores.csv",
    },
}

MODEL_COLUMNS = {
    "METH": ("METH_Score_Prob", "METH_Predicted_Label"),
    "CNA": ("CNA_Score_Prob", "CNA_Predicted_Label"),
    "FRAG": ("FRAG_Score_Prob", "FRAG_Predicted_Label"),
    "Multi": ("Multi_Score_Prob", "Multi_Predicted_Label"),
}

MODEL_COLORS = {
    "METH": "#009E73",
    "CNA": "#D55E00",
    "FRAG": "#0072B2",
    "Multi": "#CC0F0C",
}

PERCENT_METRICS = [
    "Accuracy",
    "Sensitivity",
    "Specificity",
    "Precision",
    "PPV",
    "NPV",
    "F1",
]
CI_METRICS = PERCENT_METRICS + ["AUROC"]

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "Arial"


def load_cluster2_dataset(cohort: str) -> pd.DataFrame:
    assignments = pd.read_csv(CLUSTER_FILE)
    assignments["Sample"] = assignments["Sample"].astype(str).str.strip()
    assignments["group"] = assignments["group"].astype(str).str.strip().str.lower()
    assignments["clusterlabel"] = (
        assignments["clusterlabel"].astype(str).str.strip().str.lower()
    )

    data = assignments[
        (assignments["Cohort"] == cohort)
        & (
            (assignments["group"] == "bpe")
            | (assignments["clusterlabel"] == "cluster2")
        )
    ].copy()
    data["y_true"] = data["group"].map({"bpe": 0, "mpe": 1})

    for model, prediction_file in PREDICTION_FILES[cohort].items():
        probability_column, label_column = MODEL_COLUMNS[model]
        prediction = pd.read_csv(prediction_file)[
            ["Sample", probability_column, label_column]
        ].copy()
        prediction["Sample"] = prediction["Sample"].astype(str).str.strip()
        prediction = prediction.rename(
            columns={
                probability_column: f"{model}_probability",
                label_column: f"{model}_prediction",
            }
        )
        data = data.merge(prediction, on="Sample", how="inner")

    return data


def load_cytology() -> pd.DataFrame:
    metadata_files = {
        "Discovery": DISCOVERY_METADATA_FILE,
        "ULP": ULP_METADATA_FILE,
    }
    tables = []

    for cohort, metadata_file in metadata_files.items():
        metadata = pd.read_csv(metadata_file)[["Sample", "CB", "LBC"]].copy()
        metadata["Sample"] = metadata["Sample"].astype(str).str.strip()
        cb = pd.to_numeric(metadata["CB"], errors="coerce")
        lbc = pd.to_numeric(metadata["LBC"], errors="coerce")
        has_result = cb.notna() | lbc.notna()
        metadata["Cytology"] = np.where(
            has_result,
            ((cb == CYTOLOGY_POSITIVE_VALUE) | (lbc == CYTOLOGY_POSITIVE_VALUE)).astype(int),
            np.nan,
        )
        metadata["Cohort"] = cohort
        tables.append(metadata[["Sample", "Cohort", "Cytology"]])

    return pd.concat(tables, ignore_index=True)


def calculate_binary_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_probability: np.ndarray | None = None,
) -> dict:
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    accuracy = (tp + tn) / (tp + tn + fp + fn)
    sensitivity = tp / (tp + fn) if tp + fn else np.nan
    specificity = tn / (tn + fp) if tn + fp else np.nan
    precision = tp / (tp + fp) if tp + fp else np.nan
    npv = tn / (tn + fn) if tn + fn else np.nan
    f1 = (
        2 * precision * sensitivity / (precision + sensitivity)
        if pd.notna(precision) and pd.notna(sensitivity) and precision + sensitivity
        else np.nan
    )
    auroc = (
        roc_auc_score(y_true, y_probability)
        if y_probability is not None and len(np.unique(y_true)) == 2
        else np.nan
    )

    return {
        "N": len(y_true),
        "BPE_N": int((y_true == 0).sum()),
        "MPE_N": int((y_true == 1).sum()),
        "TP": int(tp),
        "FP": int(fp),
        "TN": int(tn),
        "FN": int(fn),
        "Accuracy": accuracy,
        "Sensitivity": sensitivity,
        "Specificity": specificity,
        "Precision": precision,
        "PPV": precision,
        "NPV": npv,
        "F1": f1,
        "AUROC": auroc,
    }


def bootstrap_ci(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_probability: np.ndarray | None,
) -> dict:
    rng = np.random.default_rng(RANDOM_STATE)
    values = {metric: [] for metric in CI_METRICS}

    for _ in range(N_BOOTSTRAP):
        indices = rng.choice(len(y_true), size=len(y_true), replace=True)
        bootstrap_probability = (
            y_probability[indices] if y_probability is not None else None
        )
        metrics = calculate_binary_metrics(
            y_true[indices],
            y_pred[indices],
            bootstrap_probability,
        )
        for metric in CI_METRICS:
            if pd.notna(metrics[metric]):
                values[metric].append(metrics[metric])

    confidence_intervals = {}
    for metric, metric_values in values.items():
        if metric_values:
            lower, upper = np.percentile(metric_values, [2.5, 97.5])
            confidence_intervals[metric] = (lower, upper, len(metric_values))
        else:
            confidence_intervals[metric] = (np.nan, np.nan, 0)

    return confidence_intervals


def format_ci(metric: str, estimate: float, lower: float, upper: float) -> str:
    if pd.isna(estimate) or pd.isna(lower) or pd.isna(upper):
        return "NA"
    if metric == "AUROC":
        return f"{estimate:.3f} ({lower:.3f}-{upper:.3f})"
    return f"{estimate * 100:.1f}% ({lower * 100:.1f}-{upper * 100:.1f}%)"


def calculate_dataset_metrics(data: pd.DataFrame, dataset: str) -> list[dict]:
    rows = []
    methods = list(MODEL_COLUMNS)
    if data["Cytology"].notna().any():
        methods.append("Cytology")

    for method in methods:
        if method == "Cytology":
            selected = data[["y_true", "Cytology"]].dropna().copy()
            y_probability = None
            probability_column = ""
            prediction_column = "Cytology"
        else:
            prediction_column = f"{method}_prediction"
            probability_column = f"{method}_probability"
            selected = data[
                ["y_true", prediction_column, probability_column]
            ].dropna().copy()
            y_probability = selected[probability_column].astype(float).to_numpy()

        y_true = selected["y_true"].astype(int).to_numpy()
        y_pred = selected[prediction_column].astype(int).to_numpy()
        metrics = calculate_binary_metrics(y_true, y_pred, y_probability)
        intervals = bootstrap_ci(y_true, y_pred, y_probability)

        row = {
            "Dataset": dataset,
            "Method": method,
            "Probability_Column": MODEL_COLUMNS[method][0] if method != "Cytology" else "",
            **metrics,
        }
        for metric in CI_METRICS:
            lower, upper, valid_n = intervals[metric]
            row[f"{metric}_CI_lower"] = lower
            row[f"{metric}_CI_upper"] = upper
            row[f"{metric}_95CI"] = format_ci(metric, metrics[metric], lower, upper)
            row[f"{metric}_Bootstrap_valid_n"] = valid_n
        rows.append(row)

    return rows


def plot_roc(data: pd.DataFrame, dataset: str) -> list[dict]:
    figure, axis = plt.subplots(figsize=(8, 8))
    auc_rows = []

    for model in ["CNA", "METH", "FRAG", "Multi"]:
        probability_column = f"{model}_probability"
        selected = data[["y_true", probability_column]].dropna().copy()
        y_true = selected["y_true"].astype(int)
        y_probability = selected[probability_column].astype(float)
        fpr, tpr, _ = roc_curve(y_true, y_probability)
        auroc = roc_auc_score(y_true, y_probability)
        axis.plot(
            1 - fpr,
            tpr,
            color=MODEL_COLORS[model],
            linewidth=2.5,
            label=f"{model} (AUC={auroc:.3f})",
        )
        auc_rows.append(
            {
                "Dataset": dataset,
                "Model": model,
                "Probability_Column": MODEL_COLUMNS[model][0],
                "N": len(selected),
                "BPE_N": int((y_true == 0).sum()),
                "MPE_N": int((y_true == 1).sum()),
                "AUROC": auroc,
            }
        )

    axis.plot([1, 0], [0, 1], color="gray", linewidth=1.5, linestyle="--")
    axis.set_xlim(1.05, -0.05)
    axis.set_ylim(-0.05, 1.05)
    axis.set_xticks([0, 0.25, 0.5, 0.75, 1])
    axis.set_yticks([0, 0.25, 0.5, 0.75, 1])
    axis.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
    axis.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
    axis.set_xlabel("Specificity", fontsize=14)
    axis.set_ylabel("Sensitivity", fontsize=14)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.legend(loc="lower right", frameon=False)
    figure.tight_layout()
    figure.savefig(
        os.path.join(OUTPUT_DIR, f"Figure_S5_{dataset}_cluster2_ROC.pdf"),
        dpi=300,
    )
    plt.close(figure)

    return auc_rows


def main() -> None:
    discovery = load_cluster2_dataset("Discovery")
    ulp = load_cluster2_dataset("ULP")
    cytology = load_cytology()

    discovery = discovery.merge(cytology, on=["Sample", "Cohort"], how="left")
    ulp = ulp.merge(cytology, on=["Sample", "Cohort"], how="left")
    combined = pd.concat([discovery, ulp], ignore_index=True)

    datasets = {
        "Discovery": discovery,
        "ULP": ulp,
        "Combined": combined,
    }
    auc_rows = []
    metric_rows = []

    for dataset, data in datasets.items():
        auc_rows.extend(plot_roc(data, dataset))
        metric_rows.extend(calculate_dataset_metrics(data, dataset))

    pd.DataFrame(auc_rows).to_csv(OUTPUT_AUC, index=False)

    metrics = pd.DataFrame(metric_rows)
    base_columns = [
        "Dataset",
        "Method",
        "Probability_Column",
        "N",
        "BPE_N",
        "MPE_N",
        "TP",
        "FP",
        "TN",
        "FN",
    ]
    metric_columns = []
    for metric in CI_METRICS:
        metric_columns.extend(
            [
                metric,
                f"{metric}_CI_lower",
                f"{metric}_CI_upper",
                f"{metric}_95CI",
                f"{metric}_Bootstrap_valid_n",
            ]
        )
    metrics = metrics[base_columns + metric_columns]

    scale_columns = []
    for metric in PERCENT_METRICS:
        scale_columns.extend([metric, f"{metric}_CI_lower", f"{metric}_CI_upper"])
    metrics[scale_columns] = metrics[scale_columns] * 100
    metrics.to_csv(OUTPUT_METRICS, index=False)

    print(f"Saved: {OUTPUT_AUC}")
    print(f"Saved: {OUTPUT_METRICS}")


if __name__ == "__main__":
    main()
