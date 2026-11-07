"""Purpose: Cna bins count for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import mannwhitneyu


INPUT_FILE = "results/03_nanopore_cna/discovery_Corrected_Copy_Number.csv"
META_FILE = "metadata/discovery_groupfile.csv"
OUT_PDF = "results/12_genome_landscape/CNA/Train_CNA_bins_count.pdf"

FIGSIZE = (5, 6)
PALETTE = {
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


def style_axis(ax):
    ax.spines["bottom"].set_linewidth(2.0)
    ax.spines["left"].set_linewidth(2.0)
    ax.tick_params(axis="y", labelsize=12, width=1.5, length=5)
    ax.tick_params(axis="x", labelsize=13, width=1.5, length=5)

    for label in ax.get_xticklabels():
        label.set_fontweight("bold")


def build_gain_loss_table(cna_file):
    df = pd.read_csv(cna_file)
    sample_columns = df.columns[3:]

    rows = []
    for col in sample_columns:
        sample_name = col.replace(".Corrected_Copy_Number", "")
        values = pd.to_numeric(df[col], errors="coerce")

        rows.append({
            "Sample": sample_name,
            "gain": (values > 2).sum(),
            "loss": (values < 2).sum(),
        })

    return pd.DataFrame(rows)


def main():
    plt.rcParams.update({
        "font.sans-serif": ["Arial"],
        "font.family": "Arial",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 2.0,
    })

    sns.set_theme(style="ticks")

    df_cna = build_gain_loss_table(INPUT_FILE)
    df_meta = pd.read_csv(META_FILE)

    df = pd.merge(df_cna, df_meta, on="Sample", how="inner")
    df = df[
        (df["split"] == "train")
        & (df["group"].isin(["mpe", "bpe"]))
    ].copy()

    df_melt = pd.melt(
        df,
        id_vars=["Sample", "group"],
        value_vars=["gain", "loss"],
        var_name="CNA_Type",
        value_name="Count",
    )

    df_melt["group_label"] = df_melt["group"].map({
        "bpe": "BPE",
        "mpe": "MPE",
    })

    fig, ax = plt.subplots(figsize=FIGSIZE)

    sns.violinplot(
        data=df_melt,
        x="CNA_Type",
        y="Count",
        hue="group_label",
        order=["gain", "loss"],
        hue_order=["BPE", "MPE"],
        palette=PALETTE,
        linewidth=1.5,
        ax=ax,
        inner=None,
        saturation=0.9,
        cut=1,
    )

    if ax.get_legend() is not None:
        ax.get_legend().remove()

    positions = [-0.2, 0.2, 0.8, 1.2]
    box_data = [
        df_melt[(df_melt["CNA_Type"] == "gain") & (df_melt["group_label"] == "BPE")]["Count"].dropna(),
        df_melt[(df_melt["CNA_Type"] == "gain") & (df_melt["group_label"] == "MPE")]["Count"].dropna(),
        df_melt[(df_melt["CNA_Type"] == "loss") & (df_melt["group_label"] == "BPE")]["Count"].dropna(),
        df_melt[(df_melt["CNA_Type"] == "loss") & (df_melt["group_label"] == "MPE")]["Count"].dropna(),
    ]

    ax.boxplot(
        box_data,
        positions=positions,
        widths=0.10,
        patch_artist=True,
        showfliers=False,
        boxprops=dict(facecolor="white", edgecolor="black", linewidth=1, alpha=0.75),
        medianprops=dict(color="black", linewidth=1),
        whiskerprops=dict(color="black", linewidth=1),
        capprops=dict(color="black", linewidth=1),
        zorder=2,
    )

    sig_y_positions = {"gain": 455, "loss": 455}

    for i, metric in enumerate(["gain", "loss"]):
        bpe_data = df[df["group"] == "bpe"][metric].dropna()
        mpe_data = df[df["group"] == "mpe"][metric].dropna()

        _, p_value = mannwhitneyu(mpe_data, bpe_data, alternative="two-sided")
        sig_label = get_significance_label(p_value)

        if sig_label not in ["ns", ""]:
            y = sig_y_positions[metric]
            x1 = i - 0.2
            x2 = i + 0.2

            ax.plot(
                [x1, x1, x2, x2],
                [y - 15, y, y, y - 15],
                lw=1.5,
                color="black",
                zorder=5,
            )

            ax.text(
                (x1 + x2) / 2,
                y + 8,
                sig_label,
                ha="center",
                va="bottom",
                fontsize=15,
                fontweight="bold",
                color="black",
                zorder=6,
            )

    ax.set_ylim(0, 500)
    ax.set_yticks(np.arange(0, 501, 100))
    ax.set_ylabel("Number of CNA bins", fontsize=13, fontweight="bold")
    ax.set_xlabel("")

    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        ["Gain bins", "Loss bins"],
        fontsize=13,
        fontweight="bold",
    )

    legend_elements = [
        mpatches.Patch(facecolor=PALETTE["BPE"], edgecolor="none", label="BPE"),
        mpatches.Patch(facecolor=PALETTE["MPE"], edgecolor="none", label="MPE"),
    ]

    leg = ax.legend(
        handles=legend_elements,
        title="",
        loc="upper right",
        frameon=False,
        fontsize=11,
    )

    for text in leg.get_texts():
        text.set_fontweight("bold")

    style_axis(ax)
    sns.despine()

    fig.subplots_adjust(left=0.22, right=0.95, bottom=0.15, top=0.92)
    fig.savefig(OUT_PDF, dpi=300)



if __name__ == "__main__":
    main()
