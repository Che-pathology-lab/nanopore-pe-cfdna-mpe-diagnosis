"""Purpose: Frag genome for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy.stats import mannwhitneyu


LEN_FILE = "results/12_genome_landscape/FRAG/length_220_270.csv"
MOTIF_FILE = "results/12_genome_landscape/FRAG/motif_1mer.csv"
META_FILE = "metadata/discovery_groupfile.csv"

OUT_DIR = "results/12_genome_landscape/FRAG"
OUT_PDF_LEN = os.path.join(OUT_DIR, "Train_Length_220_270.pdf")
OUT_PDF_C = os.path.join(OUT_DIR, "Train_Motif_C.pdf")

PALETTE = {
    "BPE": "#0D75BB",
    "MPE": "#d62728",
}

GROUP_MAP = {
    "bpe": "BPE",
    "mpe": "MPE",
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


def prepare_plot_df(feature_df, meta_df, metric):
    df = pd.merge(feature_df, meta_df, on="Sample", how="inner")

    df["split"] = df["split"].astype(str).str.strip().str.lower()
    df["group"] = df["group"].astype(str).str.strip().str.lower()

    df = df[
        (df["split"].isin(["train", "panelofnormal"]))
        & (df["group"].isin(["mpe", "bpe"]))
    ].copy()

    df["group_label"] = df["group"].map(GROUP_MAP)

    return df.dropna(subset=[metric, "group_label"]).copy()


def plot_violin_box(df_plot, metric, out_pdf, ylabel, ylim):
    bpe_data = df_plot[df_plot["group_label"] == "BPE"][metric].dropna()
    mpe_data = df_plot[df_plot["group_label"] == "MPE"][metric].dropna()

    _, p_value = mannwhitneyu(mpe_data, bpe_data, alternative="two-sided")
    sig_label = get_significance_label(p_value)

    sns.set_theme(style="ticks")
    fig, ax = plt.subplots(figsize=(5, 6))

    sns.violinplot(
        data=df_plot,
        x="group_label",
        y=metric,
        order=["BPE", "MPE"],
        palette=PALETTE,
        linewidth=1.5,
        ax=ax,
        inner=None,
        saturation=0.9,
        cut=2,
    )

    sns.boxplot(
        data=df_plot,
        x="group_label",
        y=metric,
        order=["BPE", "MPE"],
        width=0.12,
        ax=ax,
        showfliers=False,
        boxprops={
            "facecolor": "white",
            "edgecolor": "black",
            "alpha": 0.75,
            "linewidth": 1,
        },
        medianprops={"color": "black", "linewidth": 1},
        whiskerprops={"color": "black", "linewidth": 1},
        capprops={"color": "black", "linewidth": 1},
    )

    ax.set_ylim(ylim)

    if sig_label not in ["ns", ""]:
        y_span = ylim[1] - ylim[0]
        y = ylim[1] - 0.08 * y_span
        drop = 0.05 * y_span
        text_offset = 0.02 * y_span

        ax.plot(
            [0, 0, 1, 1],
            [y - drop, y, y, y - drop],
            lw=1.5,
            color="black",
        )

        ax.text(
            0.5,
            y + text_offset,
            sig_label,
            ha="center",
            va="bottom",
            fontsize=15,
            fontweight="bold",
            color="black",
        )

    ax.set_ylabel(ylabel, fontsize=13, fontweight="bold")
    ax.set_xlabel("")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["BPE", "MPE"], fontsize=13, fontweight="bold")

    legend = ax.legend(
        handles=[
            mpatches.Patch(facecolor=PALETTE["BPE"], edgecolor="none", label="BPE"),
            mpatches.Patch(facecolor=PALETTE["MPE"], edgecolor="none", label="MPE"),
        ],
        title="",
        loc="upper right",
        frameon=False,
        fontsize=11,
    )

    for text in legend.get_texts():
        text.set_fontweight("bold")

    style_axis(ax)
    sns.despine()

    fig.subplots_adjust(left=0.22, right=0.95, bottom=0.15, top=0.92)
    fig.savefig(out_pdf, dpi=300)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    plt.rcParams.update({
        "font.sans-serif": ["Arial"],
        "font.family": "Arial",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 2.0,
    })

    df_meta = pd.read_csv(META_FILE)
    df_meta["Sample"] = df_meta["Sample"].astype(str).str.strip()

    df_len = pd.read_csv(LEN_FILE)
    df_len["Sample"] = df_len["Sample"].astype(str).str.strip()

    df_motif = pd.read_csv(MOTIF_FILE)
    df_motif["Sample"] = df_motif["Sample"].astype(str).str.strip()

    if "Freq_220_270" in df_len.columns:
        df_plot_len = prepare_plot_df(df_len, df_meta, "Freq_220_270")
        plot_violin_box(
            df_plot=df_plot_len,
            metric="Freq_220_270",
            out_pdf=OUT_PDF_LEN,
            ylabel="Frequency of 220-270 bp",
            ylim=(0, 0.25),
        )

    if "C" in df_motif.columns:
        df_plot_c = prepare_plot_df(df_motif, df_meta, "C")
        plot_violin_box(
            df_plot=df_plot_c,
            metric="C",
            out_pdf=OUT_PDF_C,
            ylabel="Frequency of C motif",
            ylim=(0, 0.6),
        )



if __name__ == "__main__":
    main()
