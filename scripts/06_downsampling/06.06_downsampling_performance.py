"""Plot Multi-model performance at original, 0.1x, and 0.01x depth."""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import wilcoxon
from sklearn.metrics import (
    accuracy_score,
    auc,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
GROUP_FILE = "metadata/discovery_groupfile.csv"
OUTPUT_DIR = "results/06_downsampling/performance"

PRED_FILES = {
    "Original": "results/05_nanopore_multimodal/discovery_Multi_prediction_results_with_scores.csv",
    "0.1x": "results/06_downsampling/Multi/downsampling_0.1x_Multi_prediction_results_with_scores.csv",
    "0.01x": "results/06_downsampling/Multi/downsampling_0.01x_Multi_prediction_results_with_scores.csv",
}

DEPTH_ORDER = ["Original", "0.1x", "0.01x"]
TEST_PAIRS = [("Original", "0.1x"), ("Original", "0.01x")]
DEPTH_COLORS = {"Original": "#0F2347", "0.1x": "#2F78B7", "0.01x": "#E8B46B"}
PROB_COLUMN = "Multi_Score_Prob"
PRED_COLUMN = "Multi_Predicted_Label"
METRIC_ORDER = ["AUC", "Accuracy", "Sensitivity", "Specificity", "PPV", "NPV", "F1-Score"]

mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "Arial",
    "axes.linewidth": 2.0,
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "axes.titlesize": 15,
    "axes.labelsize": 14,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 12,
    "figure.dpi": 120,
    "savefig.dpi": 600,
})


def save_pdf(fig, filename):
    fig.savefig(OUTPUT_DIR / filename, bbox_inches="tight", dpi=300)
    plt.close(fig)


def load_data():
    meta = pd.read_csv(GROUP_FILE)[["Sample", "group", "split"]].copy()
    meta["Sample"] = meta["Sample"].astype(str).str.strip()
    meta["group"] = meta["group"].astype(str).str.lower().str.strip()
    meta["split"] = meta["split"].astype(str).str.lower().str.strip()

    all_data = []
    for depth, file_path in PRED_FILES.items():
        df = pd.read_csv(file_path)
        df["Sample"] = df["Sample"].astype(str).str.strip()
        df["Multi_Prob"] = df[PROB_COLUMN].astype(float)
        df["Multi_Pred"] = df[PRED_COLUMN].astype(int)

        if depth == "Original":
            df["Base_Sample"] = df["Sample"]
            df["Replicate"] = 1
        else:
            df["Base_Sample"] = df["Original_Sample"].astype(str).str.strip()

        df = df.drop(columns=df.columns.intersection(["group", "split"]))
        df = df.merge(meta, left_on="Base_Sample", right_on="Sample", how="left", suffixes=("_pred", ""))
        df["Depth"] = depth
        all_data.append(df)

    df_all = pd.concat(all_data, ignore_index=True)
    df_all["group"] = df_all["group"].astype(str).str.lower().str.strip()
    df_all["split"] = df_all["split"].astype(str).str.lower().str.strip()
    df_all = df_all[(df_all["split"] == "test") & df_all["group"].isin(["mpe", "bpe"])].copy()
    df_all["y_true"] = (df_all["group"] == "mpe").astype(int)
    df_all["Depth"] = pd.Categorical(df_all["Depth"], DEPTH_ORDER, ordered=True)
    return df_all


def calculate_results(df_all):
    metric_rows = []
    mean_fpr = np.linspace(0, 1, 1000)
    roc_data = {depth: {"tprs": [], "tpr_at_fpr0": []} for depth in DEPTH_ORDER}

    for (depth, replicate), group_df in df_all.groupby(["Depth", "Replicate"], observed=True):
        y_true = group_df["y_true"].astype(int).to_numpy()
        y_prob = group_df["Multi_Prob"].astype(float).to_numpy()
        y_pred = group_df["Multi_Pred"].astype(int).to_numpy()

        auc_value = roc_auc_score(y_true, y_prob)
        fpr, tpr, _ = roc_curve(y_true, y_prob)
        interpolated_tpr = np.interp(mean_fpr, fpr, tpr)
        interpolated_tpr[0] = 0.0
        roc_data[depth]["tprs"].append(interpolated_tpr)
        roc_data[depth]["tpr_at_fpr0"].append(tpr[fpr == 0][-1])

        metrics = {
            "AUC": auc_value,
            "Accuracy": accuracy_score(y_true, y_pred),
            "Sensitivity": recall_score(y_true, y_pred, pos_label=1, zero_division=0),
            "Specificity": recall_score(y_true, y_pred, pos_label=0, zero_division=0),
            "PPV": precision_score(y_true, y_pred, pos_label=1, zero_division=0),
            "NPV": precision_score(y_true, y_pred, pos_label=0, zero_division=0),
            "F1-Score": f1_score(y_true, y_pred, pos_label=1, zero_division=0),
        }
        for metric, score in metrics.items():
            metric_rows.append({"Depth": depth, "Replicate": replicate, "Metric": metric, "Score": score})

    metrics_df = pd.DataFrame(metric_rows)
    metrics_df["Depth"] = pd.Categorical(metrics_df["Depth"], DEPTH_ORDER, ordered=True)
    metrics_df["Metric"] = pd.Categorical(metrics_df["Metric"], METRIC_ORDER, ordered=True)
    return metrics_df, roc_data, mean_fpr


def significance_label(p_value):
    if pd.isna(p_value):
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def perform_wilcoxon_tests(metrics_df):
    rows = []
    for metric in METRIC_ORDER:
        metric_df = metrics_df[metrics_df["Metric"] == metric]
        original_score = metric_df.loc[metric_df["Depth"] == "Original", "Score"].to_numpy(dtype=float)[0]

        for _, downsample_depth in TEST_PAIRS:
            downsample_scores = metric_df.loc[
                metric_df["Depth"] == downsample_depth, "Score"
            ].to_numpy(dtype=float)
            differences = downsample_scores - original_score
            nonzero_differences = differences[~np.isclose(differences, 0.0, atol=1e-12)]
            p_value = (
                1.0 if len(nonzero_differences) == 0
                else wilcoxon(nonzero_differences, alternative="two-sided", method="auto").pvalue
            )

            rows.append({
                "Metric": metric,
                "Depth_A": "Original",
                "Depth_B": downsample_depth,
                "N_Downsampling_Repeats": len(downsample_scores),
                "Original_Score": original_score,
                "Downsampling_Mean": np.mean(downsample_scores),
                "Downsampling_SD": np.std(downsample_scores, ddof=1) if len(downsample_scores) > 1 else 0.0,
                "Mean_Difference_Low_minus_Original": np.mean(differences),
                "N_Lower_Than_Original": int(np.sum(differences < 0)),
                "N_Equal_To_Original": int(np.sum(np.isclose(differences, 0.0, atol=1e-12))),
                "N_Higher_Than_Original": int(np.sum(differences > 0)),
                "P_Raw_Two_Sided": p_value,
            })

    results_df = pd.DataFrame(rows)
    results_df["P_Adj_Holm_By_Metric"] = np.nan

    for metric in METRIC_ORDER:
        mask = (results_df["Metric"] == metric) & results_df["P_Raw_Two_Sided"].notna()
        adjusted_p = multipletests(results_df.loc[mask, "P_Raw_Two_Sided"], method="holm")[1]
        results_df.loc[mask, "P_Adj_Holm_By_Metric"] = adjusted_p

    results_df["Significance"] = results_df["P_Adj_Holm_By_Metric"].apply(significance_label)
    mean_difference = results_df["Mean_Difference_Low_minus_Original"]
    results_df["Change_Direction"] = np.select(
        [mean_difference < 0, mean_difference > 0], ["Decrease", "Increase"], default="No change"
    )
    results_df["Significant_Change"] = results_df["P_Adj_Holm_By_Metric"] < 0.05
    results_df["Significant_Decrease"] = (
        results_df["Significant_Change"] & (results_df["Change_Direction"] == "Decrease")
    )
    return results_df


def add_significance_brackets(
    ax, metric, significance_df, x_positions, base_y, step_y, bracket_drop, show_ns=True
):
    metric_df = significance_df[significance_df["Metric"] == metric]
    level = 0

    for depth_a, depth_b in TEST_PAIRS:
        row = metric_df[(metric_df["Depth_A"] == depth_a) & (metric_df["Depth_B"] == depth_b)]
        label = row.iloc[0]["Significance"]
        if not label or (label == "ns" and not show_ns):
            continue

        x1, x2 = x_positions[depth_a], x_positions[depth_b]
        y = base_y + level * step_y
        ax.plot(
            [x1, x1, x2, x2], [y - bracket_drop, y, y, y - bracket_drop],
            color="black", linewidth=0.9, clip_on=False, zorder=10,
        )
        ax.text(
            (x1 + x2) / 2, y + bracket_drop * 0.35, label,
            ha="center", va="bottom", fontsize=7.5, fontweight="bold", clip_on=False, zorder=11,
        )
        level += 1


def style_open_axes(ax, linewidth=1.5):
    ax.spines["bottom"].set_linewidth(linewidth)
    ax.spines["left"].set_linewidth(linewidth)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_auc(metrics_df, significance_df):
    auc_df = metrics_df[metrics_df["Metric"] == "AUC"].copy()
    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    common = {
        "data": auc_df, "x": "Depth", "y": "Score", "hue": "Depth", "order": DEPTH_ORDER,
        "hue_order": DEPTH_ORDER, "palette": DEPTH_COLORS, "legend": False, "ax": ax,
    }
    sns.boxplot(**common, width=0.52, linewidth=1.2, fliersize=0)
    sns.stripplot(
        **common, alpha=0.85, size=4.2, jitter=0.07,
        edgecolor="black", linewidth=0.3, zorder=3,
    )

    add_significance_brackets(
        ax, "AUC", significance_df, {d: i for i, d in enumerate(DEPTH_ORDER)},
        base_y=1.05, step_y=0.06, bracket_drop=0.014,
    )
    ax.set(xlim=(-0.5, 2.5), ylim=(0.4, 1.16), xlabel="Sequencing Depth", ylabel="AUROC")
    ax.xaxis.label.set(fontsize=11, fontweight="bold")
    ax.yaxis.label.set(fontsize=11, fontweight="bold")
    ax.tick_params(axis="x", rotation=15, labelsize=9, width=1.5, length=4.5)
    ax.tick_params(axis="y", labelsize=9, width=1.5, length=4.5)
    for label in ax.get_xticklabels():
        label.set_fontweight("bold")
    style_open_axes(ax)
    fig.tight_layout()
    save_pdf(fig, "Figure1_Test_AUROC_Boxplot_TwoSided_Wilcoxon.pdf")


def plot_roc(roc_data, mean_fpr):
    fig, ax = plt.subplots(figsize=(4.2, 4.2))

    for depth in DEPTH_ORDER:
        tprs = np.asarray(roc_data[depth]["tprs"])
        mean_tpr = tprs.mean(axis=0)
        mean_tpr[-1] = 1.0
        std_tpr = tprs.std(axis=0)
        tpr_at_zero = np.mean(roc_data[depth]["tpr_at_fpr0"])
        std_at_zero = np.std(roc_data[depth]["tpr_at_fpr0"])

        plot_fpr = np.insert(mean_fpr, 1, 0.0)
        plot_tpr = np.insert(mean_tpr, 1, tpr_at_zero)
        plot_std = np.insert(std_tpr, 1, std_at_zero)
        plot_tpr[0] = plot_std[0] = 0.0
        specificity = 1 - plot_fpr
        mean_auc = auc(mean_fpr, mean_tpr)

        ax.plot(
            specificity, plot_tpr, color=DEPTH_COLORS[depth], linewidth=2.0,
            label=f"{depth} (AUC={mean_auc:.3f})",
        )
        if depth != "Original":
            ax.fill_between(
                specificity, np.clip(plot_tpr - plot_std, 0, 1), np.clip(plot_tpr + plot_std, 0, 1),
                color=DEPTH_COLORS[depth], alpha=0.13, linewidth=0,
            )

    ax.plot([1, 0], [0, 1], color="gray", linewidth=1.2, linestyle="--", alpha=0.7)
    ticks = [0, 0.25, 0.50, 0.75, 1.00]
    ax.set(xlim=(1.05, -0.05), ylim=(-0.05, 1.05), xticks=ticks, yticks=ticks,
           xlabel="Specificity", ylabel="Sensitivity")
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
    ax.xaxis.label.set(fontsize=11, fontweight="bold")
    ax.yaxis.label.set(fontsize=11, fontweight="bold")
    ax.tick_params(width=1.5, length=4.5, labelsize=9)
    ax.legend(loc="lower right", fontsize=8.5, frameon=False, handlelength=1.4,
              handletextpad=0.4, labelspacing=0.35)
    style_open_axes(ax)
    fig.tight_layout()
    save_pdf(fig, "Figure2_Test_ROC_Curves.pdf")


def plot_clinical_metrics(metrics_df, significance_df):
    clinical_order = ["Accuracy", "Sensitivity", "Specificity", "PPV", "NPV", "F1-Score"]
    plot_df = metrics_df[metrics_df["Metric"].isin(clinical_order)].copy()
    summary = plot_df.groupby(["Metric", "Depth"], observed=True)["Score"].agg(["mean", "std"]).reset_index()

    sns.set_theme(style="ticks", rc={"font.family": "Arial"})
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.0), sharey=True)
    x_positions = np.arange(len(DEPTH_ORDER)) * 0.70

    for ax, metric in zip(axes.flat, clinical_order):
        metric_summary = summary[summary["Metric"] == metric].set_index("Depth").reindex(DEPTH_ORDER)
        means = metric_summary["mean"].to_numpy()
        standard_deviations = metric_summary["std"].fillna(0).to_numpy()
        ax.bar(
            x_positions, means, yerr=standard_deviations, width=0.42,
            color=[DEPTH_COLORS[d] for d in DEPTH_ORDER], edgecolor="black", linewidth=0.55, zorder=2,
            error_kw={"elinewidth": 0.9, "ecolor": "black", "capsize": 2.5, "capthick": 0.9},
        )
        add_significance_brackets(
            ax, metric, significance_df, dict(zip(DEPTH_ORDER, x_positions)),
            base_y=1.08, step_y=0.07, bracket_drop=0.018,
        )
        ax.set_title(metric, fontsize=9, fontweight="bold")
        ax.set(xlim=(x_positions[0] - 0.35, x_positions[-1] + 0.35), ylim=(0, 1.23),
               yticks=np.arange(0, 1.01, 0.25), xticks=x_positions)
        ax.set_xticklabels(DEPTH_ORDER, rotation=35, ha="right", fontsize=8, fontweight="bold")
        ax.tick_params(axis="y", labelsize=9, width=0.8, length=3)
        style_open_axes(ax, linewidth=1.0)

    axes[0, 0].set_ylabel("Metric value", fontsize=9, fontweight="bold")
    axes[1, 0].set_ylabel("Metric value", fontsize=9, fontweight="bold")
    fig.subplots_adjust(wspace=0.18, hspace=0.35)
    save_pdf(fig, "Figure3_Test_Clinical_Metrics_Facet_Barplot_TwoSided_Wilcoxon.pdf")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df_all = load_data()
    metrics_df, roc_data, mean_fpr = calculate_results(df_all)
    significance_df = perform_wilcoxon_tests(metrics_df)

    plot_auc(metrics_df, significance_df)
    plot_roc(roc_data, mean_fpr)
    plot_clinical_metrics(metrics_df, significance_df)


if __name__ == "__main__":
    main()
