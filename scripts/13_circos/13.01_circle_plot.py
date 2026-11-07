"""Purpose: Circle plot for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
from tqdm import tqdm

mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6.0,
    "axes.linewidth": 0.6,
    "savefig.dpi": 600,
    "savefig.pad_inches": 0,
})

AXIS_LINEWIDTH = 0.55
METH_TRACE_LINEWIDTH = 0.4
METH_RED = "#E05A34"
METH_BLUE = "#2E86AB"
CNA_RED = "#C0392B"
CNA_BLUE = "#2C4A8C"
CNA_NEUTRAL = "#111111"
FRAG_TEAL = "#7ECFC3"
FRAG_MINT = "#4CB8A8"
DEPTH_GRAY = "#7F7F7F"
BG_COLOR = "#FFFFFF"
GRID_COLOR = "#D4D4D4"

meth_file = "results/13_circos/ONT_train/Train_methylation_5kb_mean.csv"
depth_file = "results/13_circos/ONT_train/Train_depth_5mb_mean.csv"
cna_file = "results/13_circos/ONT_train/Train_apparent_cn_mean.csv"
motif_file = "results/13_circos/ONT_train/Train_C_motif_5mb_mean.csv"
length_file = "results/13_circos/ONT_train/Train_length_220_270_5mb_mean.csv"
output_pdf = "results/13_circos/Train_Genome_wide_Multi_Circos.pdf"


def main():
    chrom_order = [f"chr{i}" for i in range(1, 23)]

    LEN_CLIP = 0.02
    MOTIF_CLIP = 0.08
    CNA_DIFF = 0.4
    METH_DIFF = 0.20
    MAX_DEPTH = 2
    CNA_MIN_BAR = 0.1
    BIN_EDGE_PAD_FRAC = 0.015

    ZERO_RS = [10.25, 13.65, 18.30, 24.20, 30.10, 33.50]
    TRACK_BANDS = [
        (9.0, 11.5),
        (12.4, 14.9),
        (15.8, 20.8),
        (21.7, 26.7),
        (27.6, 32.6),
        (33.5, 35.0),
    ]
    TRACK_LABELS = [
        "Frag. size\nΔ",
        "C-EM Δ",
        "Methylation\nΔβ",
        "CN (BPE)",
        "CN (MPE)",
        "Coverage",
    ]

    df_meth = pd.read_csv(meth_file)
    df_depth = pd.read_csv(depth_file)
    df_cna = pd.read_csv(cna_file)
    df_motif = pd.read_csv(motif_file)
    df_len = pd.read_csv(length_file)

    for data in [df_meth, df_depth, df_cna, df_motif, df_len]:
        data[["chr", "start", "end"]] = data["chr_start_end"].str.split(
            "_", expand=True
        )

    def prepare(data):
        data = data[data["chr"].isin(chrom_order)].copy()
        data["start"] = data["start"].astype(int)
        data["end"] = data["end"].astype(int)
        data["chr"] = pd.Categorical(
            data["chr"], categories=chrom_order, ordered=True
        )
        return data.sort_values(["chr", "start"]).reset_index(drop=True)

    df_meth, df_depth, df_cna, df_motif, df_len = map(
        prepare, [df_meth, df_depth, df_cna, df_motif, df_len]
    )

    df_meth["mpe_smooth"] = (
        df_meth.groupby("chr", observed=True)["mpe"]
        .transform(lambda values: values.rolling(40, center=True, min_periods=1).mean())
        / 100
    )
    df_meth["bpe_smooth"] = (
        df_meth.groupby("chr", observed=True)["bpe"]
        .transform(lambda values: values.rolling(40, center=True, min_periods=1).mean())
        / 100
    )
    df_meth["diff_smooth"] = df_meth["mpe_smooth"] - df_meth["bpe_smooth"]
    df_meth.loc[
        df_meth.groupby("chr", observed=True)["start"].diff().fillna(0) > 10000,
        ["mpe_smooth", "bpe_smooth", "diff_smooth"],
    ] = np.nan

    chrom_lengths = df_meth.groupby("chr", observed=True)["start"].max()
    gap_bp = chrom_lengths.sum() * 0.006
    chrom_offsets = {}
    current_offset = 0
    for chrom in chrom_order:
        chrom_offsets[chrom] = current_offset
        current_offset += chrom_lengths.get(chrom, 0) + gap_bp
    total_circle_len = current_offset
    use_radians = 1.75 * np.pi

    def to_theta(position, chrom):
        return (
            (position + chrom_offsets[str(chrom)])
            / total_circle_len
            * use_radians
        )

    df_meth["theta"] = [
        to_theta(start, chrom)
        for start, chrom in zip(df_meth["start"], df_meth["chr"].astype(str))
    ]

    def add_bin_theta(data):
        chrom_length = data["chr"].astype(str).map(chrom_lengths).astype(float)
        start = data["start"].clip(lower=0)
        end = np.minimum(data["end"].clip(lower=0), chrom_length)
        theta_start = np.array([
            to_theta(position, chrom)
            for position, chrom in zip(start, data["chr"].astype(str))
        ])
        theta_end = np.array([
            to_theta(position, chrom)
            for position, chrom in zip(end, data["chr"].astype(str))
        ])
        padding = (theta_end - theta_start) * BIN_EDGE_PAD_FRAC
        data["theta_s"] = theta_start + padding
        data["theta_e"] = theta_end - padding

    for data in [df_cna, df_motif, df_len]:
        add_bin_theta(data)

    df_meth["r_diff"] = ZERO_RS[2] + np.clip(
        df_meth["diff_smooth"] / METH_DIFF, -1, 1
    ) * 2.5

    axis_theta = 0.0
    figure_size = 98 / 25.4
    fig = plt.figure(figsize=(figure_size, figure_size), facecolor=BG_COLOR)
    ax = fig.add_axes(
        [0.025, 0.015, 0.95, 0.86], projection="polar", facecolor=BG_COLOR
    )
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 38)
    ax.axis("off")

    full_theta = np.linspace(0, use_radians, 800)
    for index, (inner_radius, outer_radius) in enumerate(TRACK_BANDS):
        ax.fill_between(
            full_theta,
            inner_radius,
            outer_radius,
            color="#EFEFEF",
            alpha=0.16 if index % 2 == 0 else 0.07,
            zorder=0,
            linewidth=0,
        )

    for chrom in tqdm(chrom_order, desc="Methylation and framework"):
        chrom_data = df_meth[df_meth["chr"] == chrom]
        if chrom_data.empty:
            continue

        theta_arc = np.linspace(
            to_theta(0, chrom), to_theta(chrom_lengths[chrom], chrom), 300
        )
        for zero_radius in ZERO_RS:
            ax.plot(
                theta_arc,
                [zero_radius] * 300,
                color=GRID_COLOR,
                lw=0.35,
                dashes=[3, 4],
                alpha=0.75,
                zorder=1,
            )

        for theta_edge in [
            to_theta(0, chrom),
            to_theta(chrom_lengths[chrom], chrom),
        ]:
            ax.plot(
                [theta_edge, theta_edge],
                [8.5, 35.0],
                color="#B8B8B8",
                lw=0.35,
                alpha=0.85,
                zorder=1,
            )

        positive = chrom_data["diff_smooth"].fillna(0) >= 0
        ax.fill_between(
            chrom_data["theta"], ZERO_RS[2], chrom_data["r_diff"],
            where=positive, color=METH_RED, alpha=1, zorder=2, linewidth=0,
        )
        ax.fill_between(
            chrom_data["theta"], ZERO_RS[2], chrom_data["r_diff"],
            where=~positive, color=METH_BLUE, alpha=1, zorder=2, linewidth=0,
        )
        ax.plot(
            chrom_data["theta"],
            np.where(positive, chrom_data["r_diff"], np.nan),
            color=METH_RED,
            lw=METH_TRACE_LINEWIDTH,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=3,
        )
        ax.plot(
            chrom_data["theta"],
            np.where(~positive, chrom_data["r_diff"], np.nan),
            color=METH_BLUE,
            lw=METH_TRACE_LINEWIDTH,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=3,
        )

        chrom_number = int(chrom.replace("chr", ""))
        ax.text(
            (theta_arc[0] + theta_arc[-1]) / 2,
            36.7,
            chrom.replace("chr", ""),
            ha="center",
            va="center",
            fontsize=6.2,
            fontweight="bold",
            color="#1A1A1A" if chrom_number % 2 == 1 else "#444444",
            clip_on=False,
        )

    def plot_diff_blocks(data, column, base_radius, maximum, color):
        for _, row in data.iterrows():
            value = row[column]
            if pd.isna(value) or value == 0:
                continue
            radial_change = np.clip(abs(value) / maximum, 0, 1) * 2.5
            target_radius = (
                base_radius + radial_change if value > 0 else base_radius - radial_change
            )
            ax.fill_between(
                [row["theta_s"], row["theta_e"]],
                base_radius,
                target_radius,
                color=color,
                alpha=0.88,
                zorder=2,
                linewidth=0,
            )

    plot_diff_blocks(df_len, "diff", ZERO_RS[0], LEN_CLIP, FRAG_TEAL)
    plot_diff_blocks(df_motif, "diff", ZERO_RS[1], MOTIF_CLIP, FRAG_MINT)

    def cna_delta_to_radius(copy_number):
        radial_change = np.clip((copy_number - 2) / CNA_DIFF, -1, 1) * 2.5
        if radial_change != 0 and abs(radial_change) < CNA_MIN_BAR:
            radial_change = np.sign(radial_change) * CNA_MIN_BAR
        return radial_change

    for _, row in df_cna.iterrows():
        theta_start = row["theta_s"]
        theta_end = row["theta_e"]

        for column, zero_radius, alpha in [
            ("bpe", ZERO_RS[3], 0.90),
            ("mpe", ZERO_RS[4], 0.85),
        ]:
            value = row.get(column, np.nan)
            if pd.isna(value):
                continue
            radial_change = cna_delta_to_radius(value)
            if radial_change == 0:
                radial_change = CNA_MIN_BAR
                color = CNA_NEUTRAL
            else:
                color = CNA_RED if radial_change > 0 else CNA_BLUE
            ax.fill_between(
                [theta_start, theta_end],
                zero_radius,
                zero_radius + radial_change,
                color=color,
                alpha=alpha,
                zorder=3 if column == "bpe" else 2,
                linewidth=0,
            )

    for _, row in tqdm(df_depth.iterrows(), total=len(df_depth), desc="Depth track"):
        chrom = str(row["chr"])
        if chrom not in chrom_offsets:
            continue
        bin_start = row["start"]
        genome_start = bin_start + chrom_offsets[chrom]
        genome_end = min(
            bin_start + 5_000_000 + chrom_offsets[chrom],
            chrom_offsets[chrom] + chrom_lengths.get(chrom, bin_start),
        )
        theta_start = genome_start / total_circle_len * use_radians
        theta_end = genome_end / total_circle_len * use_radians
        depth_radius = ZERO_RS[5] + (
            np.clip(row["mean_depth"], 0, MAX_DEPTH) / MAX_DEPTH
        ) * 1.5
        ax.fill_between(
            [theta_start, theta_end], ZERO_RS[5], depth_radius,
            color=DEPTH_GRAY, alpha=0.60, zorder=2, linewidth=0,
        )

    for inner_radius, outer_radius in TRACK_BANDS:
        ax.plot(
            [axis_theta, axis_theta],
            [inner_radius, outer_radius],
            color="#222222",
            lw=AXIS_LINEWIDTH,
            zorder=5,
            solid_capstyle="round",
        )

    ticks = [
        (9.0, f"−{LEN_CLIP}"), (10.25, "0"), (11.5, f"+{LEN_CLIP}"),
        (12.4, f"−{MOTIF_CLIP}"), (13.65, "0"), (14.9, f"+{MOTIF_CLIP}"),
        (15.8, f"−{METH_DIFF}"), (18.3, "0"), (20.8, f"+{METH_DIFF}"),
        (21.7, "1.6"), (24.2, "2.0"), (26.7, "2.4"),
        (27.6, "1.6"), (30.1, "2.0"), (32.6, "2.4"),
        (33.5, "0"), (35.0, str(MAX_DEPTH)),
    ]
    for radius, label in ticks:
        ax.text(
            axis_theta - 0.045,
            radius,
            label,
            ha="right",
            va="center",
            fontsize=4.8,
            color="#333333",
            clip_on=False,
        )

    open_sector = 2 * np.pi - use_radians
    label_theta = axis_theta - open_sector / 2
    label_rotation = np.degrees(open_sector / 2)
    for label, (inner_radius, outer_radius) in zip(TRACK_LABELS, TRACK_BANDS):
        ax.text(
            label_theta,
            (inner_radius + outer_radius) / 2,
            label,
            ha="center",
            va="center",
            rotation=label_rotation,
            rotation_mode="anchor",
            fontsize=5.5,
            fontweight="bold",
            color="#1A1A1A",
            linespacing=1.05,
            clip_on=False,
        )

    legend_elements = [
        mpatches.Patch(facecolor=CNA_RED, label="CN gain"),
        mpatches.Patch(facecolor=CNA_BLUE, label="CN loss"),
        mpatches.Patch(facecolor=METH_RED, label="Δβ > 0"),
        mpatches.Patch(facecolor=METH_BLUE, label="Δβ < 0"),
        mpatches.Patch(facecolor=FRAG_MINT, label="ΔC-EM"),
        mpatches.Patch(facecolor=FRAG_TEAL, label="Δ220–270-bp fraction"),
        mpatches.Patch(facecolor=DEPTH_GRAY, label="Coverage"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        ncol=4,
        frameon=False,
        fontsize=5.3,
        handlelength=0.9,
        handleheight=0.75,
        handletextpad=0.35,
        columnspacing=0.85,
        labelspacing=0.25,
        borderaxespad=0,
    )

    fig.savefig(
        output_pdf,
        format="pdf",
        facecolor=BG_COLOR,
        transparent=False,
        bbox_inches=None,
        pad_inches=0,
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
