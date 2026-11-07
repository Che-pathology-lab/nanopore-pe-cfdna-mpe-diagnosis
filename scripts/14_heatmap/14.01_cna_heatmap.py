"""Purpose: Cna heatmap for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial"]
mpl.rcParams["font.weight"] = "bold"
mpl.rcParams["axes.labelweight"] = "bold"
mpl.rcParams["axes.titleweight"] = "bold"

CNA_FILE = "results/03_nanopore_cna/discovery_Corrected_Copy_Number.csv"
GROUP_FILE = "metadata/discovery_groupfile.csv"
OUT_DIR = "results/14_heatmap"
OUT_PDF = os.path.join(OUT_DIR, "Train_CNA_subgroup_heatmap.pdf")

GROUP_ORDER = ["mpe", "bpe"]
SUBGROUP_ORDER = ["luad", "meso", "sclc", "scc", "tb", "ppe"]

group_label_map = {
    "mpe": "MPE",
    "bpe": "BPE",
}

group_colors = {
    "MPE": "#d62728",
    "BPE": "#0D75BB",
}

subgroup_colors = {
    "luad": "#D97C7C",
    "meso": "#B79AD6",
    "sclc": "#E6B65C",
    "scc": "#C07BA4",
    "tb": "#79B7D8",
    "ppe": "#7DBE9F",
}


def main():
    df_cna = pd.read_csv(CNA_FILE)
    df_group = pd.read_csv(GROUP_FILE)

    df_group["Sample"] = df_group["Sample"].astype(str).str.strip()
    df_group["split"] = df_group["split"].astype(str).str.strip().str.lower()
    df_group["group"] = df_group["group"].astype(str).str.strip().str.lower()
    df_group["subgroup"] = df_group["subgroup"].astype(str).str.strip().str.lower()

    df_train = df_group[
        (df_group["split"] == "train") &
        (df_group["group"].isin(GROUP_ORDER)) &
        (df_group["subgroup"].isin(SUBGROUP_ORDER))
    ].copy()

    df_train["Group"] = df_train["group"].map(group_label_map)

    df_train["group_rank"] = pd.Categorical(
        df_train["group"],
        categories=GROUP_ORDER,
        ordered=True
    )

    df_train["subgroup_rank"] = pd.Categorical(
        df_train["subgroup"],
        categories=SUBGROUP_ORDER,
        ordered=True
    )

    df_train = df_train.sort_values(
        ["group_rank", "subgroup_rank", "Sample"]
    ).reset_index(drop=True)

    pos_cols = df_cna.columns[:3].tolist()
    sample_raw_cols = df_cna.columns[3:].tolist()

    rename_map = {
        col: col.replace(".Corrected_Copy_Number", "")
        for col in sample_raw_cols
    }
    df_cna = df_cna.rename(columns=rename_map)

    sample_cols = [
        s for s in df_train["Sample"].tolist()
        if s in df_cna.columns
    ]
    sample_info = df_train[df_train["Sample"].isin(sample_cols)].copy()
    sample_info = sample_info.set_index("Sample").loc[sample_cols].reset_index()

    region_id = (
        df_cna[pos_cols[0]].astype(str) + "_" +
        df_cna[pos_cols[1]].astype(str) + "_" +
        df_cna[pos_cols[2]].astype(str)
    )

    mat = df_cna[sample_cols].copy()
    mat.index = region_id
    mat = mat.apply(pd.to_numeric, errors="coerce")

    mat_filled = mat.copy()
    row_median = mat_filled.median(axis=1)
    mat_filled = mat_filled.T.fillna(row_median).T

    col_median = mat_filled.median(axis=0)
    mat_filled = mat_filled.fillna(col_median)
    mat_filled = mat_filled.fillna(2.0)

    row_linkage = linkage(
        pdist(mat_filled.values, metric="euclidean"),
        method="average"
    )

    cmap = LinearSegmentedColormap.from_list(
        "paper_blue_white_red",
        ["#86C5DA", "#F0EEE9", "#E64B35"],
        N=256
    )
    norm = TwoSlopeNorm(vmin=1.0, vcenter=2.0, vmax=5.0)

    col_colors = pd.DataFrame({
        "Group": sample_info["Group"].map(group_colors).values,
        "Subgroup": sample_info["subgroup"].map(subgroup_colors).values,
    }, index=sample_info["Sample"])

    sns.set_theme(style="white", font_scale=0.75)

    g = sns.clustermap(
        mat_filled,
        row_cluster=True,
        col_cluster=False,
        row_linkage=row_linkage,
        col_colors=col_colors,
        cmap=cmap,
        norm=norm,
        xticklabels=False,
        yticklabels=False,
        figsize=(12, 10),
        dendrogram_ratio=(0.001, 0.001),
        cbar_kws={
            "label": "Corrected copy number",
            "ticks": [1, 2, 5],
        },
        cbar_pos=(0.98, 0.42, 0.018, 0.14),
        colors_ratio=(0.01, 0.014),
    )

    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)

    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("ichorCNA bins")
    g.ax_col_colors.set_ylabel("")

    g.ax_col_colors.axhline(
        1,
        color="white",
        linewidth=1.2,
        zorder=20
    )

    for spine in g.ax_col_colors.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.4)
        spine.set_color("#333333")

    g.ax_cbar.tick_params(labelsize=6, width=0.4, length=2)
    g.ax_cbar.set_ylabel(
        "Corrected copy number",
        fontsize=9,
        fontweight="bold",
        labelpad=4
    )

    for spine in g.ax_cbar.spines.values():
        spine.set_linewidth(0.4)

    group_handles = [
        Patch(facecolor=group_colors["MPE"], edgecolor="none", label="MPE"),
        Patch(facecolor=group_colors["BPE"], edgecolor="none", label="BPE"),
    ]

    subgroup_handles = [
        Patch(facecolor=subgroup_colors[x], edgecolor="none", label=x.upper())
        for x in SUBGROUP_ORDER
    ]

    leg1 = g.ax_heatmap.legend(
        handles=group_handles,
        title="Group",
        loc="upper left",
        bbox_to_anchor=(1.04, 1.00),
        frameon=False,
        borderaxespad=0,
        fontsize=8,
        title_fontsize=8,
    )

    g.ax_heatmap.add_artist(leg1)

    g.ax_heatmap.legend(
        handles=subgroup_handles,
        title="Subgroup",
        loc="upper left",
        bbox_to_anchor=(1.04, 0.95),
        frameon=False,
        borderaxespad=0,
        fontsize=8,
        title_fontsize=8,
    )

    plt.savefig(OUT_PDF, dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
