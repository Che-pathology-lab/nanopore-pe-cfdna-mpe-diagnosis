"""Plot patient-level recurrence and burden of downsampling errors."""

import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator


GROUP_FILE = "metadata/discovery_groupfile.csv"
OUTPUT_DIR = "results/06_downsampling/error_recurrence"

PRED_FILES = {
    "Original": "results/05_nanopore_multimodal/discovery_Multi_prediction_results_with_scores.csv",
    "0.1x": "results/06_downsampling/Multi/downsampling_0.1x_Multi_prediction_results_with_scores.csv",
    "0.01x": "results/06_downsampling/Multi/downsampling_0.01x_Multi_prediction_results_with_scores.csv",
}

PRED_COLUMN = "Multi_Predicted_Label"
PATIENT_FIGURE_PDF = f"{OUTPUT_DIR}/FigureA_downsampling_patient_level_misclassification.pdf"
BURDEN_FIGURE_PDF = f"{OUTPUT_DIR}/FigureB_downsampling_error_burden_across_repeats.pdf"

DEPTH_ORDER = ["0.1x", "0.01x"]
DEPTH_COLORS = {"0.1x": "#2F78B7", "0.01x": "#E8B46B"}
GROUP_COLORS = {"BPE": "#0D75BB", "MPE": "#d62728"}
NEW_FP_COLOR = GROUP_COLORS["BPE"]
NEW_FN_COLOR = GROUP_COLORS["MPE"]

FIGSIZE = (4.2, 4.2)
AXIS_LINEWIDTH = 1.2
LABEL_FONTSIZE = 10
TITLE_FONTSIZE = 10
TICK_FONTSIZE = 8
LEGEND_FONTSIZE = 8
TICK_WIDTH = 1.0
TICK_LENGTH = 4.0

mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "Arial",
    "font.size": TICK_FONTSIZE,
    "axes.linewidth": AXIS_LINEWIDTH,
    "axes.labelweight": "bold",
    "axes.titleweight": "bold",
    "axes.titlesize": TITLE_FONTSIZE,
    "axes.labelsize": LABEL_FONTSIZE,
    "xtick.labelsize": TICK_FONTSIZE,
    "ytick.labelsize": TICK_FONTSIZE,
    "legend.fontsize": LEGEND_FONTSIZE,
    "figure.dpi": 120,
    "savefig.dpi": 600,
})


def load_original_prediction():
    df = pd.read_csv(PRED_FILES["Original"])
    df["Sample"] = df["Sample"].astype(str).str.strip()
    return pd.DataFrame({
        "Base_Sample": df["Sample"],
        "Original_Predicted_Label": df[PRED_COLUMN].astype(int),
    })


def load_downsampling_prediction(file_path, depth):
    df = pd.read_csv(file_path)
    return pd.DataFrame({
        "Base_Sample": df["Original_Sample"].astype(str).str.strip(),
        "Replicate": df["Replicate"].astype(str).str.strip(),
        "Depth": depth,
        "Predicted_Label": df[PRED_COLUMN].astype(int),
    })


def prepare_event_data():
    metadata = pd.read_csv(GROUP_FILE)
    metadata["Sample"] = metadata["Sample"].astype(str).str.strip()
    metadata["group"] = metadata["group"].astype(str).str.lower().str.strip()
    metadata["split"] = metadata["split"].astype(str).str.lower().str.strip()
    metadata = metadata[
        (metadata["split"] == "test") & metadata["group"].isin(["mpe", "bpe"])
    ].copy()
    metadata["Base_Sample"] = metadata["Sample"]
    metadata["y_true"] = (metadata["group"] == "mpe").astype(int)

    original = load_original_prediction()
    baseline = metadata.merge(
        original[["Base_Sample", "Original_Predicted_Label"]], on="Base_Sample", how="inner"
    )

    event_tables = []
    for depth in DEPTH_ORDER:
        events = load_downsampling_prediction(PRED_FILES[depth], depth).merge(
            baseline[["Base_Sample", "y_true", "group", "Original_Predicted_Label"]],
            on="Base_Sample",
            how="inner",
        )
        y_true = events["y_true"]
        original_pred = events["Original_Predicted_Label"]
        downsample_pred = events["Predicted_Label"]
        events["Original_Error"] = original_pred.ne(y_true)
        events["Downsampling_Error"] = downsample_pred.ne(y_true)
        event_tables.append(events)

    return pd.concat(event_tables, ignore_index=True)


def make_replicate_summary(events):
    summary = (
        events.groupby(["Depth", "Replicate"], observed=True)
        .agg(Original_Errors=("Original_Error", "sum"), Downsampling_Errors=("Downsampling_Error", "sum"))
        .reset_index()
    )
    summary["Net_Additional_Errors"] = summary["Downsampling_Errors"] - summary["Original_Errors"]
    return summary.sort_values(["Depth", "Replicate"])


def make_sample_summary(events):
    return (
        events.groupby(["Depth", "Base_Sample"], observed=True)
        .agg(
            group=("group", "first"),
            N_Repeats=("Replicate", "size"),
            Original_Error=("Original_Error", "first"),
            N_Downsampling_Error=("Downsampling_Error", "sum"),
        )
        .reset_index()
    )


def add_frequency_marker(ax, x, y, count, total, color, fontsize=None):
    if count == 0:
        ax.scatter(
            x, y, s=57, marker="o", facecolor="#E8E8E8", edgecolor="#B8B8B8",
            linewidth=0.6, zorder=3,
        )
        return

    ax.scatter(
        x, y, s=165 + 41 * count, marker="o", facecolor=color,
        edgecolor="white", linewidth=0.8, zorder=4,
    )
    ax.text(
        x, y, f"{100 * count / total:.0f}", ha="center", va="center",
        fontsize=fontsize or (6.7 if count < 10 else 6.2),
        fontweight="bold", color="white", zorder=5,
    )


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)
    ax.spines["left"].set_linewidth(AXIS_LINEWIDTH)
    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE, width=TICK_WIDTH, length=TICK_LENGTH)


def plot_error_summary(sample_summary, replicate_summary):
    original_error = (
        sample_summary.groupby("Base_Sample", observed=True)["Original_Error"].first().astype(bool)
    )
    downsampling_error = sample_summary.groupby("Base_Sample", observed=True)["N_Downsampling_Error"].sum()
    affected_samples = original_error.index[original_error | downsampling_error.gt(0)]
    if len(affected_samples) == 0:
        return

    sample_plot = sample_summary[
        sample_summary["Base_Sample"].isin(affected_samples) & sample_summary["Depth"].isin(DEPTH_ORDER)
    ].copy()
    error_n = sample_plot.pivot(
        index="Base_Sample", columns="Depth", values="N_Downsampling_Error"
    ).reindex(columns=DEPTH_ORDER, fill_value=0).fillna(0)
    depth_repeat_n = sample_plot.groupby("Depth", observed=True)["N_Repeats"].max().astype(int).to_dict()
    sample_group = sample_plot.drop_duplicates("Base_Sample").set_index("Base_Sample")["group"]
    sample_order = (
        error_n.assign(
            Original_Error=original_error.reindex(error_n.index).astype(int),
            Total=error_n.sum(axis=1),
            Sample=error_n.index,
        )
        .sort_values(
            ["Original_Error", "Total", "0.01x", "0.1x", "Sample"],
            ascending=[False, False, False, False, True],
        )
        .index.tolist()
    )

    sns.set_theme(style="ticks", rc={"font.family": "Arial"})
    figure_height = max(FIGSIZE[1], 0.42 * len(sample_order) + 1.55)
    shared_figsize = (FIGSIZE[0], figure_height)
    fig_sample, ax_sample = plt.subplots(figsize=shared_figsize, facecolor="white")
    fig_burden, ax_burden = plt.subplots(figsize=shared_figsize, facecolor="white")

    y_positions = np.arange(len(sample_order))
    y_labels = [
        f"{re.sub(r'_ONT_DISC$', '', str(sample))} "
        f"({'BPE' if sample_group.loc[sample] == 'bpe' else 'MPE'})"
        for sample in sample_order
    ]
    depth_positions = {"0.1x": 1, "0.01x": 2}

    for y, sample in zip(y_positions, sample_order):
        color = NEW_FP_COLOR if sample_group.loc[sample] == "bpe" else NEW_FN_COLOR
        original_count = 10 if original_error.loc[sample] else 0
        add_frequency_marker(ax_sample, 0, y, original_count, 10, color, fontsize=6.4)

        for depth in DEPTH_ORDER:
            count = int(error_n.loc[sample, depth])
            add_frequency_marker(
                ax_sample, depth_positions[depth], y, count, depth_repeat_n[depth], color
            )

    ax_sample.set_yticks(y_positions, labels=y_labels, rotation=45, ha="right", va="center",
                         rotation_mode="anchor", fontsize=7)
    ax_sample.invert_yaxis()
    ax_sample.set_xlim(-0.55, 2.55)
    ax_sample.set_xticks(
        [0, 1, 2],
        labels=[
            "Original\n(n=1)",
            f"0.1×\n(n={depth_repeat_n['0.1x']})",
            f"0.01×\n(n={depth_repeat_n['0.01x']})",
        ],
    )
    ax_sample.set_xlabel("Sequencing depth", fontsize=LABEL_FONTSIZE, fontweight="bold")
    ax_sample.set_ylabel("Sample", fontsize=LABEL_FONTSIZE, fontweight="bold")
    ax_sample.grid(axis="x", color="#E0E0E0", linewidth=0.65)
    ax_sample.grid(axis="y", visible=False)

    n_original_errors = int(original_error.reindex(sample_order).sum())
    if 0 < n_original_errors < len(sample_order):
        ax_sample.axhline(
            n_original_errors - 0.5, color="#8C8C8C", linewidth=0.8,
            linestyle=(0, (3, 2)), zorder=1,
        )

    for depth_index, depth in enumerate(DEPTH_ORDER):
        depth_data = replicate_summary[replicate_summary["Depth"] == depth].copy()
        depth_data["_sort"] = depth_data["Replicate"].str.extract(r"(\d+)", expand=False).astype(int)
        depth_data = depth_data.sort_values("_sort")
        values = depth_data["Net_Additional_Errors"].to_numpy(dtype=float)
        jitter = np.linspace(-0.11, 0.11, len(values))
        color = DEPTH_COLORS[depth]

        ax_burden.scatter(
            depth_index + jitter, values, s=36, color=color, edgecolor="white",
            linewidth=0.55, alpha=0.90, zorder=3,
        )
        mean_value = values.mean()
        sd_value = values.std(ddof=1) if len(values) > 1 else 0
        ax_burden.errorbar(
            depth_index, mean_value, yerr=sd_value, fmt="D", markersize=5.5,
            color=color, markeredgecolor="#202020", markeredgewidth=0.6,
            capsize=4, linewidth=1.3, zorder=5,
        )
        ax_burden.annotate(
            f"Mean = {mean_value:.1f}", (depth_index, mean_value + sd_value), xytext=(0, 8),
            textcoords="offset points", ha="center", va="bottom", fontsize=8,
            fontweight="bold", color=color,
        )

    ax_burden.axhline(0, color="#4D4D4D", linewidth=0.9, linestyle="--")
    ax_burden.set_xticks([0, 1], labels=["0.1×", "0.01×"])
    ax_burden.set_xlim(-0.45, 1.45)
    ax_burden.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax_burden.set_xlabel("Sequencing depth", fontsize=LABEL_FONTSIZE, fontweight="bold")
    ax_burden.set_ylabel(
        "Change in the number of misclassified samples\n(downsampled - original)",
        fontsize=LABEL_FONTSIZE,
        fontweight="bold",
    )
    ax_burden.grid(axis="y", color="#D9D9D9", linewidth=0.65)
    ax_burden.grid(axis="x", visible=False)

    patient_legend = [
        Line2D(
            [0], [0], marker="o", color="none", markerfacecolor="#E8E8E8",
            markeredgecolor="#B8B8B8", markersize=6, label="No error",
        ),
        Patch(facecolor=NEW_FP_COLOR, edgecolor="none", label="False positive (BPE)"),
        Patch(facecolor=NEW_FN_COLOR, edgecolor="none", label="False negative (MPE)"),
    ]
    ax_sample.legend(
        handles=patient_legend, title="Misclassification frequency (%)", loc="upper left",
        bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=LEGEND_FONTSIZE,
        title_fontsize=LEGEND_FONTSIZE, borderaxespad=0, handletextpad=0.55, labelspacing=0.45,
    )

    burden_legend = [
        Line2D(
            [0], [0], marker="o", color="none", markerfacecolor=DEPTH_COLORS[depth],
            markeredgecolor="white", markersize=6, label=depth.replace("x", "×"),
        )
        for depth in DEPTH_ORDER
    ]
    ax_burden.legend(
        handles=burden_legend, loc="upper left", bbox_to_anchor=(1.02, 1.0),
        frameon=False, fontsize=LEGEND_FONTSIZE, borderaxespad=0,
        handletextpad=0.55, labelspacing=0.45,
    )

    style_axes(ax_sample)
    style_axes(ax_burden)
    ax_sample.tick_params(axis="x", width=TICK_WIDTH, length=TICK_LENGTH, pad=6)
    ax_sample.tick_params(axis="y", labelsize=7, width=TICK_WIDTH, length=TICK_LENGTH, pad=5)

    fig_sample.subplots_adjust(left=0.29, right=0.97, bottom=min(0.24, 1.0 / figure_height), top=0.97)
    fig_burden.subplots_adjust(left=0.22, right=0.97, bottom=0.16, top=0.97)
    fig_sample.savefig(PATIENT_FIGURE_PDF, dpi=300, bbox_inches="tight", facecolor="white")
    fig_burden.savefig(BURDEN_FIGURE_PDF, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig_sample)
    plt.close(fig_burden)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    events = prepare_event_data()
    replicate_summary = make_replicate_summary(events)
    sample_summary = make_sample_summary(events)
    plot_error_summary(sample_summary, replicate_summary)


if __name__ == "__main__":
    main()
