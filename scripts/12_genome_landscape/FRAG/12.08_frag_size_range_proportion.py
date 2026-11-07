"""Purpose: Frag size range proportion for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")

FEATURE_PATH = "results/12_genome_landscape/FRAG/ONT_lengths_matrix.csv"
LABEL_PATH = "metadata/discovery_groupfile.csv"
OUTPUT_DIR = "results/12_genome_landscape/FRAG"
OUT_PDF = os.path.join(OUTPUT_DIR, "fragment_range_train.pdf")

FRAGMENT_RANGES = [
    ("short_mono", 100, 150, "100-150 bp\nShort"),
    ("canonical_mono", 150, 170, "150-170 bp\nCanonical mono"),
    ("long_mono", 170, 220, "170-220 bp\nLong mono"),
    ("mono_di_transition", 220, 270, "220-270 bp\nMono-di transition"),
    ("di_nucleosomal", 270, 320, "270-320 bp\nDi-nucleosomal"),
]

COLOR_DICT = {
    "BPE": "#0D75BB",
    "MPE": "#d62728",
}


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


def load_data():
    df_features = pd.read_csv(FEATURE_PATH, index_col=0)
    df_labels = pd.read_csv(LABEL_PATH)

    df_features.index = df_features.index.astype(str).str.strip()
    df_labels["Sample"] = df_labels["Sample"].astype(str).str.strip()
    df_labels["split"] = df_labels["split"].astype(str).str.strip().str.lower()
    df_labels["group"] = df_labels["group"].astype(str).str.strip().str.lower()

    df_labels = df_labels[
        df_labels["split"].isin(["train", "panelofnormal"])
        & df_labels["group"].isin(["mpe", "bpe"])
    ].copy()

    df_labels["group_label"] = df_labels["group"].map({
        "mpe": "MPE",
        "bpe": "BPE",
    })

    df = df_labels.merge(
        df_features,
        left_on="Sample",
        right_index=True,
        how="inner",
    )

    return df, df_features.columns


def calculate_range_features(df, feature_cols):
    length_cols = {}

    for col in feature_cols:
        try:
            length_cols[col] = int(float(col))
        except ValueError:
            continue

    for metric, start, end, _ in FRAGMENT_RANGES:
        cols = [
            col for col, length in length_cols.items()
            if start <= length < end
        ]

        if len(cols) == 0:
            raise ValueError(f"No columns found for {metric}")

        df[metric] = df[cols].apply(pd.to_numeric, errors="coerce").sum(axis=1)

    return df


def run_statistics(df):
    rows = []

    for metric, _, _, label in FRAGMENT_RANGES:
        bpe = df.loc[df["group_label"] == "BPE", metric].dropna()
        mpe = df.loc[df["group_label"] == "MPE", metric].dropna()

        if len(bpe) > 0 and len(mpe) > 0:
            _, p_value = mannwhitneyu(mpe, bpe, alternative="two-sided")
        else:
            p_value = np.nan

        rows.append({
            "metric": metric,
            "label": label,
            "p_value": p_value,
        })

    return pd.DataFrame(rows)


def add_box_overlay(ax, plot_df, metrics):
    dodge_offset = 0.21
    box_width = 0.055

    for i, metric in enumerate(metrics):
        for group, offset in zip(["BPE", "MPE"], [-dodge_offset, dodge_offset]):
            values = plot_df[
                (plot_df["metric"] == metric)
                & (plot_df["group_label"] == group)
            ]["value"].dropna()

            if values.empty:
                continue

            pos = i + offset

            q1 = values.quantile(0.25)
            q2 = values.quantile(0.50)
            q3 = values.quantile(0.75)
            iqr = q3 - q1

            lower = values[values >= q1 - 1.5 * iqr].min()
            upper = values[values <= q3 + 1.5 * iqr].max()

            ax.add_patch(
                plt.Rectangle(
                    (pos - box_width / 2, q1),
                    box_width,
                    q3 - q1,
                    facecolor="white",
                    edgecolor="black",
                    linewidth=0.8,
                    alpha=0.75,
                    zorder=4,
                )
            )

            ax.plot(
                [pos - box_width / 2, pos + box_width / 2],
                [q2, q2],
                color="black",
                linewidth=0.9,
                zorder=5,
            )

            ax.plot([pos, pos], [lower, q1], color="black", linewidth=0.8, zorder=4)
            ax.plot([pos, pos], [q3, upper], color="black", linewidth=0.8, zorder=4)

            cap_width = box_width * 0.75

            ax.plot(
                [pos - cap_width / 2, pos + cap_width / 2],
                [lower, lower],
                color="black",
                linewidth=0.8,
                zorder=4,
            )

            ax.plot(
                [pos - cap_width / 2, pos + cap_width / 2],
                [upper, upper],
                color="black",
                linewidth=0.8,
                zorder=4,
            )


def add_significance(ax, plot_df, stats_df, metrics):
    dodge_offset = 0.21

    y_min = plot_df["value"].min()
    y_max = plot_df["value"].max()
    y_range = y_max - y_min

    if y_range == 0:
        y_range = y_max if y_max > 0 else 1

    ax.set_ylim(max(0, y_min - y_range * 0.05), y_max + y_range * 0.30)

    for i, metric in enumerate(metrics):
        p_value = stats_df.loc[stats_df["metric"] == metric, "p_value"].iloc[0]
        sig_label = get_significance_label(p_value)

        if sig_label in ["ns", ""]:
            continue

        local_y_max = plot_df.loc[plot_df["metric"] == metric, "value"].max()

        x1 = i - dodge_offset
        x2 = i + dodge_offset

        line_y = local_y_max + y_range * 0.06
        line_h = y_range * 0.025
        text_y = line_y + line_h + y_range * 0.005

        ax.plot(
            [x1, x1, x2, x2],
            [line_y, line_y + line_h, line_y + line_h, line_y],
            lw=1.2,
            color="black",
            zorder=6,
        )

        ax.text(
            i,
            text_y,
            sig_label,
            ha="center",
            va="bottom",
            fontsize=14,
            fontweight="bold",
            color="black",
        )


def plot_fragment_range_violin(df, stats_df):
    sns.set_theme(style="ticks")

    metrics = [x[0] for x in FRAGMENT_RANGES]
    labels = [x[3] for x in FRAGMENT_RANGES]
    label_map = dict(zip(metrics, labels))

    plot_df = df[["group_label"] + metrics].melt(
        id_vars="group_label",
        value_vars=metrics,
        var_name="metric",
        value_name="value",
    )

    plot_df["fragment_range"] = plot_df["metric"].map(label_map)

    plt.rcParams.update({
        "font.family": "Arial",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(14, 5.8))

    sns.violinplot(
        data=plot_df,
        x="fragment_range",
        y="value",
        hue="group_label",
        order=labels,
        hue_order=["BPE", "MPE"],
        palette=COLOR_DICT,
        linewidth=1.2,
        inner=None,
        saturation=0.9,
        cut=2,
        dodge=True,
        width=0.85,
        ax=ax,
    )

    add_box_overlay(ax, plot_df, metrics)

    handles, legend_labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[:2],
        legend_labels[:2],
        title="",
        frameon=False,
        loc="upper right",
        fontsize=11,
    )

    add_significance(ax, plot_df, stats_df, metrics)

    ax.set_xlabel("")
    ax.set_ylabel("Fragment proportion", fontsize=12, fontweight="bold")

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(
        labels,
        rotation=0,
        ha="center",
        fontsize=10,
        fontweight="bold",
    )

    ax.spines["bottom"].set_linewidth(1.6)
    ax.spines["left"].set_linewidth(1.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(OUT_PDF, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df, feature_cols = load_data()
    df = calculate_range_features(df, feature_cols)
    stats_df = run_statistics(df)

    plot_fragment_range_violin(df, stats_df)



if __name__ == "__main__":
    main()
