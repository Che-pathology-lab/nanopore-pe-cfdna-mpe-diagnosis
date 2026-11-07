"""Purpose: Overall head to head for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from sklearn.metrics import confusion_matrix
from statsmodels.stats.contingency_tables import mcnemar
from statsmodels.stats.multitest import multipletests


mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "Arial"


MODEL_FILE = "results/09_ulp_validation/ULP_Multi_prediction_results_with_scores.csv"
GROUP_FILE = "metadata/ULP_groupfile.csv"
IHC_FILE = "results/10_cytology_comparison/ihc_result.csv"

OUTPUT_DIR = "results/10_cytology_comparison"

CYTO_MODEL_PDF = os.path.join(
    OUTPUT_DIR,
    "ULP_cytology_vs_multimodal_model_metrics.pdf",
)

COMBINED_GUIDED_PDF = os.path.join(
    OUTPUT_DIR,
    "ULP_combined_vs_model_guided_cytology_metrics.pdf",
)


FIGSIZE = (5.2, 5.2)
YLIM = (0, 142)

METRICS_ORDER = ["Accuracy", "Sensitivity", "Specificity"]

COMMON_COLORS = {
    "Cell Block": "#aec7e8",
    "Liquid-based Cytology": "#5d99c4",
    "Combined Cytology": "#0D75BB",
    "Multimodal Model": "#d62728",
    "Model-guided Cytology": "#D99A00",
}


def calculate_metrics(y_true, y_pred):
    y_true = np.asarray(y_true).astype(int)
    y_pred = np.asarray(y_pred).astype(int)

    if len(y_true) == 0:
        return {
            "Accuracy": np.nan,
            "Sensitivity": np.nan,
            "Specificity": np.nan,
        }

    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    tn, fp, fn, tp = cm.ravel()

    accuracy = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) > 0 else np.nan
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    specificity = tn / (tn + fp) if (tn + fp) > 0 else np.nan

    return {
        "Accuracy": accuracy,
        "Sensitivity": sensitivity,
        "Specificity": specificity,
    }


def perform_mcnemar(y_true, y_pred_a, y_pred_b):
    n01 = np.sum((y_pred_a == y_true) & (y_pred_b != y_true))
    n10 = np.sum((y_pred_a != y_true) & (y_pred_b == y_true))
    table = [[0, n01], [n10, 0]]

    if (n01 + n10) < 25:
        result = mcnemar(table, exact=True)
    else:
        result = mcnemar(table, exact=False, correction=True)

    return result.pvalue


def get_significance_label(p):
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def plot_metric_bars(metrics_df, methods_order, compare_pairs, significance_map, out_pdf):
    os.makedirs(os.path.dirname(out_pdf), exist_ok=True)

    fig, ax = plt.subplots(figsize=FIGSIZE)

    x = np.arange(len(METRICS_ORDER))

    if len(methods_order) > 2:
        bar_width = 0.15
        bar_step = 0.19
    else:
        bar_width = 0.23
        bar_step = 0.30

    offsets = (
        np.arange(len(methods_order)) - (len(methods_order) - 1) / 2
    ) * bar_step

    bar_centers = {}

    for i, method in enumerate(methods_order):
        values = metrics_df.loc[method, METRICS_ORDER].values.astype(float)

        bars = ax.bar(
            x + offsets[i],
            values,
            width=bar_width,
            color=COMMON_COLORS[method],
            edgecolor="white",
            linewidth=0.9,
            label=method,
            zorder=2,
        )

        for metric, bar, value in zip(METRICS_ORDER, bars, values):
            center = bar.get_x() + bar.get_width() / 2
            bar_centers[(metric, method)] = center

            if pd.notna(value) and value > 0:
                ax.text(
                    center,
                    value + 1.8,
                    f"{value:.1f}",
                    ha="center",
                    va="bottom",
                    fontsize=6.8,
                    fontweight="bold",
                    color="#222222",
                    zorder=6,
                )

    bracket_drop = 1.6
    bracket_step = 7.0

    for metric in METRICS_ORDER:
        metric_pairs = []

        for m1, m2 in compare_pairs:
            sig = significance_map.get((metric, m1, m2), "")

            if sig and sig != "ns":
                x1 = bar_centers[(metric, m1)]
                x2 = bar_centers[(metric, m2)]
                span = abs(x2 - x1)
                metric_pairs.append((span, x1, x2, sig))

        metric_pairs = sorted(metric_pairs, key=lambda z: z[0])
        base_y = 112

        for level, (_, x1, x2, sig) in enumerate(metric_pairs):
            y = base_y + level * bracket_step

            ax.plot(
                [x1, x1, x2, x2],
                [y - bracket_drop, y, y, y - bracket_drop],
                lw=1.0,
                color="black",
                zorder=5,
            )

            ax.text(
                (x1 + x2) / 2,
                y + 0.45,
                sig,
                ha="center",
                va="bottom",
                fontsize=10.5,
                fontweight="bold",
                color="black",
                zorder=6,
            )

    ax.set_xlim(-0.55, len(METRICS_ORDER) - 0.45)
    ax.set_ylim(*YLIM)
    ax.set_yticks(np.arange(0, 101, 20))

    ax.set_xticks(x)
    ax.set_xticklabels(METRICS_ORDER, fontsize=9.5, fontweight="bold")

    ax.set_ylabel("Percentage (%)", fontsize=10.5, fontweight="bold")
    ax.set_xlabel("")

    ax.tick_params(axis="y", labelsize=10, width=1.2, length=4)
    ax.tick_params(axis="x", width=1.2, length=0)

    ax.yaxis.grid(
        True,
        linestyle="-",
        linewidth=0.5,
        color="#E8E8E8",
        zorder=0,
    )
    ax.xaxis.grid(False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["bottom"].set_zorder(10)
    ax.spines["left"].set_zorder(10)
    ax.spines["bottom"].set_linewidth(1.35)
    ax.spines["left"].set_linewidth(1.35)

    legend = ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.19),
        ncol=2,
        frameon=False,
        fontsize=7.8,
        handlelength=1.15,
        columnspacing=1.1,
    )

    for text in legend.get_texts():
        text.set_fontweight("bold")

    fig.subplots_adjust(left=0.16, right=0.98, bottom=0.13, top=0.80)
    fig.savefig(out_pdf, dpi=300)


def load_ulp_data():
    df_model = pd.read_csv(MODEL_FILE)
    df_group = pd.read_csv(GROUP_FILE)
    df_ihc = pd.read_csv(IHC_FILE)

    for df in [df_model, df_group, df_ihc]:
        df["Sample"] = df["Sample"].astype(str).str.strip()

    df_group["split"] = df_group["split"].astype(str).str.strip().str.lower()
    df_group["group"] = df_group["group"].astype(str).str.strip().str.lower()

    df_merge = pd.merge(
        df_group[["Sample", "group", "split", "CB", "LBC", "Combined_Cyto"]],
        df_model[["Sample", "Multi_Predicted_Label"]],
        on="Sample",
        how="inner",
    )

    df_merge = df_merge[df_merge["split"] == "ulp"].copy()
    df_merge = df_merge[df_merge["group"].isin(["mpe", "bpe"])].copy()

    df_merge["y_true"] = df_merge["group"].map({
        "mpe": 1,
        "bpe": 0,
    })

    df_merge["CB_Pred"] = pd.to_numeric(df_merge["CB"], errors="coerce")
    df_merge["LBC_Pred"] = pd.to_numeric(df_merge["LBC"], errors="coerce")
    df_merge["Combined_Pred"] = pd.to_numeric(df_merge["Combined_Cyto"], errors="coerce")
    df_merge["Multi_Predicted_Label"] = pd.to_numeric(
        df_merge["Multi_Predicted_Label"],
        errors="coerce",
    )

    df_ihc["ihc"] = pd.to_numeric(df_ihc["ihc"], errors="coerce")

    df_merge = df_merge.merge(
        df_ihc[["Sample", "ihc"]],
        on="Sample",
        how="left",
    )

    ihc_target_mask = (
        (df_merge["Multi_Predicted_Label"].astype(float) == 1)
        & (df_merge["Combined_Pred"].astype(float) == 0)
    )

    bpe_ihc_default_mask = (
        ihc_target_mask
        & (df_merge["y_true"] == 0)
        & (df_merge["ihc"].isna())
    )

    df_merge.loc[bpe_ihc_default_mask, "ihc"] = 0

    df_merge["Cytology_IHC_Pred"] = df_merge["Combined_Pred"].copy()
    df_merge.loc[ihc_target_mask & (df_merge["ihc"] == 1), "Cytology_IHC_Pred"] = 1
    df_merge.loc[ihc_target_mask & (df_merge["ihc"] == 0), "Cytology_IHC_Pred"] = 0

    return df_merge


def build_cyto_model_plot(df_merge):
    methods_order = [
        "Cell Block",
        "Liquid-based Cytology",
        "Combined Cytology",
        "Multimodal Model",
    ]

    pred_cols = {
        "Cell Block": "CB_Pred",
        "Liquid-based Cytology": "LBC_Pred",
        "Combined Cytology": "Combined_Pred",
        "Multimodal Model": "Multi_Predicted_Label",
    }

    results = {}

    for method, pred_col in pred_cols.items():
        mask = ~df_merge[pred_col].isna() & ~df_merge["y_true"].isna()
        y_true = df_merge.loc[mask, "y_true"].astype(int).values
        y_pred = df_merge.loc[mask, pred_col].astype(int).values
        results[method] = calculate_metrics(y_true, y_pred)

    metrics_df = pd.DataFrame(results).T * 100
    metrics_df = metrics_df.reindex(methods_order)
    metrics_df = metrics_df.reindex(columns=METRICS_ORDER)

    methods_compare = [
        "Cell Block",
        "Liquid-based Cytology",
        "Combined Cytology",
    ]

    raw_p_values = {
        "Accuracy": {},
        "Sensitivity": {},
        "Specificity": {},
    }

    for method in methods_compare:
        pred_col = pred_cols[method]

        mask = (
            ~df_merge["y_true"].isna()
            & ~df_merge["Multi_Predicted_Label"].isna()
            & ~df_merge[pred_col].isna()
        )

        y_true = df_merge.loc[mask, "y_true"].astype(int).values
        y_model = df_merge.loc[mask, "Multi_Predicted_Label"].astype(int).values
        y_cyto = df_merge.loc[mask, pred_col].astype(int).values

        for metric in METRICS_ORDER:
            if metric == "Sensitivity":
                metric_mask = y_true == 1
            elif metric == "Specificity":
                metric_mask = y_true == 0
            else:
                metric_mask = np.ones_like(y_true, dtype=bool)

            if np.sum(metric_mask) > 0:
                raw_p_values[metric][method] = perform_mcnemar(
                    y_true[metric_mask],
                    y_model[metric_mask],
                    y_cyto[metric_mask],
                )
            else:
                raw_p_values[metric][method] = np.nan

    significance_map = {}

    for metric in METRICS_ORDER:
        pvals = [raw_p_values[metric].get(m, np.nan) for m in methods_compare]
        valid_indices = [i for i, p in enumerate(pvals) if not pd.isna(p)]
        valid_pvals = [pvals[i] for i in valid_indices]

        if valid_pvals:
            _, pvals_corrected, _, _ = multipletests(
                valid_pvals,
                alpha=0.05,
                method="holm",
            )

            for i, p_corr in zip(valid_indices, pvals_corrected):
                method = methods_compare[i]
                significance_map[(metric, method, "Multimodal Model")] = (
                    get_significance_label(p_corr)
                )

    compare_pairs = [
        ("Combined Cytology", "Multimodal Model"),
        ("Liquid-based Cytology", "Multimodal Model"),
        ("Cell Block", "Multimodal Model"),
    ]

    plot_metric_bars(
        metrics_df=metrics_df,
        methods_order=methods_order,
        compare_pairs=compare_pairs,
        significance_map=significance_map,
        out_pdf=CYTO_MODEL_PDF,
    )


def build_combined_guided_plot(df_merge):
    methods_order = [
        "Combined Cytology",
        "Model-guided Cytology",
    ]

    pred_cols = {
        "Combined Cytology": "Combined_Pred",
        "Model-guided Cytology": "Cytology_IHC_Pred",
    }

    results = {}

    for method, pred_col in pred_cols.items():
        mask = ~df_merge[pred_col].isna() & ~df_merge["y_true"].isna()
        y_true = df_merge.loc[mask, "y_true"].astype(int).values
        y_pred = df_merge.loc[mask, pred_col].astype(int).values
        results[method] = calculate_metrics(y_true, y_pred)

    metrics_df = pd.DataFrame(results).T * 100
    metrics_df = metrics_df.reindex(methods_order)
    metrics_df = metrics_df.reindex(columns=METRICS_ORDER)

    mask = (
        ~df_merge["Combined_Pred"].isna()
        & ~df_merge["Cytology_IHC_Pred"].isna()
        & ~df_merge["y_true"].isna()
    )

    y_true = df_merge.loc[mask, "y_true"].astype(int).values
    y_combined = df_merge.loc[mask, "Combined_Pred"].astype(int).values
    y_guided = df_merge.loc[mask, "Cytology_IHC_Pred"].astype(int).values

    significance_map = {}

    for metric in METRICS_ORDER:
        if metric == "Sensitivity":
            metric_mask = y_true == 1
        elif metric == "Specificity":
            metric_mask = y_true == 0
        else:
            metric_mask = np.ones_like(y_true, dtype=bool)

        if np.sum(metric_mask) > 0:
            p_value = perform_mcnemar(
                y_true[metric_mask],
                y_combined[metric_mask],
                y_guided[metric_mask],
            )
        else:
            p_value = np.nan

        significance_map[
            (metric, "Combined Cytology", "Model-guided Cytology")
        ] = get_significance_label(p_value)

    compare_pairs = [
        ("Combined Cytology", "Model-guided Cytology"),
    ]

    plot_metric_bars(
        metrics_df=metrics_df,
        methods_order=methods_order,
        compare_pairs=compare_pairs,
        significance_map=significance_map,
        out_pdf=COMBINED_GUIDED_PDF,
    )


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df_merge = load_ulp_data()

    build_cyto_model_plot(df_merge)
    build_combined_guided_plot(df_merge)


if __name__ == "__main__":
    main()