"""Purpose: 1mer endmotif plot for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import os
import warnings
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")

FEATURE_PATH = "results/12_genome_landscape/FRAG/motif_1mer.csv"
GROUP_PATH = "metadata/discovery_groupfile.csv"
OUTPUT_DIR = "results/12_genome_landscape/FRAG"
OUT_PDF = os.path.join(OUTPUT_DIR, "1mer_EM.pdf")

MOTIF_FEATURES = [
    ("A", "A"),
    ("C", "C"),
    ("G", "G"),
    ("T", "T"),
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
    df_features = pd.read_csv(FEATURE_PATH)
    df_labels = pd.read_csv(GROUP_PATH)

    df_features["Sample"] = df_features["Sample"].astype(str).str.strip()
    df_labels["Sample"] = df_labels["Sample"].astype(str).str.strip()
    df_labels["split"] = df_labels["split"].astype(str).str.strip().str.lower()
    df_labels["group"] = df_labels["group"].astype(str).str.strip().str.lower()

    motif_cols = [feature for feature, _ in MOTIF_FEATURES]

    df_labels = df_labels[
        df_labels["split"].isin(["train", "panelofnormal"])
        & df_labels["group"].isin(["mpe", "bpe"])
    ].copy()

    df_labels["group_label"] = df_labels["group"].map({
        "mpe": "MPE",
        "bpe": "BPE",
    })

    missing_cols = [c for c in ["Sample"] + motif_cols if c not in df_features.columns]
    if missing_cols:
        raise ValueError(f"Missing columns in motif file: {missing_cols}")

    for col in motif_cols:
        df_features[col] = pd.to_numeric(df_features[col], errors="coerce")

    return df_labels.merge(
        df_features[["Sample"] + motif_cols],
        on="Sample",
        how="inner",
    )


def run_statistics(df):
    rows = []

    for metric, label in MOTIF_FEATURES:
        bpe = df.loc[df["group_label"] == "BPE", metric].dropna()
        mpe = df.loc[df["group_label"] == "MPE", metric].dropna()

        if len(bpe) > 0 and len(mpe) > 0:
            _, p_value = mannwhitneyu(mpe, bpe, alternative="two-sided")
        else:
            p_value = pd.NA

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


def add_significance(ax, stats_df, metrics):
    dodge_offset = 0.21

    for i, metric in enumerate(metrics):
        p_value = stats_df.loc[stats_df["metric"] == metric, "p_value"].iloc[0]
        sig_label = get_significance_label(p_value)

        if sig_label in ["ns", ""]:
            continue

        x1 = i - dodge_offset
        x2 = i + dodge_offset

        line_y = 0.535
        line_h = 0.018
        text_y = line_y + line_h + 0.006

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


def plot_motif_violin(df, stats_df):
    sns.set_theme(style="ticks")

    metrics = [feature for feature, _ in MOTIF_FEATURES]
    labels = [label for _, label in MOTIF_FEATURES]
    label_map = dict(zip(metrics, labels))

    plot_df = df[["group_label"] + metrics].melt(
        id_vars="group_label",
        value_vars=metrics,
        var_name="metric",
        value_name="value",
    )

    plot_df["motif"] = plot_df["metric"].map(label_map)

    plt.rcParams.update({
        "font.family": "Arial",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(8.5, 5.8))

    sns.violinplot(
        data=plot_df,
        x="motif",
        y="value",
        hue="group_label",
        order=labels,
        hue_order=["BPE", "MPE"],
        palette=COLOR_DICT,
        linewidth=1.2,
        inner=None,
        saturation=1,
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
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        fontsize=11,
    )

    ax.set_ylim(0, 0.6)
    add_significance(ax, stats_df, metrics)

    ax.set_xlabel("")
    ax.set_ylabel("1-mer End Motif Frequency", fontsize=12, fontweight="bold")

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(
        labels,
        rotation=0,
        ha="center",
        fontsize=12,
        fontweight="bold",
    )

    ax.spines["bottom"].set_linewidth(1.6)
    ax.spines["left"].set_linewidth(1.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout(rect=[0, 0, 0.86, 1])
    fig.savefig(OUT_PDF, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    df = load_data()
    stats_df = run_statistics(df)
    plot_motif_violin(df, stats_df)



if __name__ == "__main__":
    main()
