import os
from itertools import product

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

CNA_FILE = "results/03_nanopore_cna/discovery_Corrected_Copy_Number.csv"
DMR_FILE = "results/02_nanopore_methylation/GLM_filtered_cpg_mean_diff_15_q_0.1_raw_data.tsv"
HYPER_BED = "results/02_nanopore_methylation/train_defined_hyper.bed"
HYPO_BED = "results/02_nanopore_methylation/train_defined_hypo.bed"
LENGTH_MATRIX_FILE = "results/12_genome_landscape/FRAG/cfdna_length_220_270_5mb_matrix.csv"
C_MOTIF_MATRIX_FILE = "results/12_genome_landscape/FRAG/cfdna_endmotif_5mb_autosomes_matrix.csv"
GROUP_FILE = "metadata/discovery_groupfile.csv"
OUTPUT_PDF = "results/12_genome_landscape/Figure_regional_abnormalities.pdf"

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]
CHROMOSOME_LENGTHS = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
}
CHROMOSOME_GAP = 5_000_000
EXPECTED_C_4MERS = {"C" + "".join(suffix) for suffix in product("ACGT", repeat=3)}

CENTROMERES = {
    "chr1": (120000000, 150000001), "chr2": (85000000, 100000001),
    "chr3": (85000000, 95000001), "chr4": (45000000, 55000001),
    "chr5": (45000000, 55000001), "chr6": (55000000, 65000001),
    "chr7": (55000000, 80000001), "chr8": (40000000, 50000001),
    "chr9": (35000000, 70000001), "chr10": (35000000, 50000001),
    "chr11": (45000000, 60000001), "chr12": (30000000, 40000001),
    "chr13": (0, 20000001), "chr14": (0, 20000001),
    "chr15": (0, 35000001), "chr16": (15000000, 50000001),
    "chr17": (15000000, 30000001), "chr18": (10000000, 25000001),
    "chr19": (20000000, 30000001), "chr20": (20000000, 35000001),
    "chr21": (0, 15000001), "chr22": (0, 25000001),
}

COLORS = {
    "gain": "#C0392B",
    "loss": "#2C4A8C",
    "hyper": "#E05A34",
    "hypo": "#2E86AB",
    "length": "#7ECFC3",
    "cem": "#4CB8A8",
    "boundary": "#B7B7B7",
    "alternate": "#F7F7F7",
}


def load_training_samples():
    group_df = pd.read_csv(GROUP_FILE)
    group_df["Sample"] = group_df["Sample"].astype(str).str.strip()
    group_df["split"] = group_df["split"].astype(str).str.strip().str.lower()
    group_df["group"] = group_df["group"].astype(str).str.strip().str.lower()
    group_df = group_df[
        (group_df["split"] == "train") & group_df["group"].isin(["mpe", "bpe"])
    ]
    mpe_samples = group_df.loc[group_df["group"] == "mpe", "Sample"].tolist()
    bpe_samples = group_df.loc[group_df["group"] == "bpe", "Sample"].tolist()
    return mpe_samples, bpe_samples


def chromosome_arm(chromosome, start, end):
    chromosome = str(chromosome)
    chromosome = chromosome if chromosome.startswith("chr") else "chr" + chromosome
    midpoint = (start + end) / 2
    return chromosome + ("p" if midpoint < CENTROMERES[chromosome][0] else "q")


def calculate_cna(mpe_samples):
    cna = pd.read_csv(CNA_FILE)
    sample_columns = [
        f"{sample}.Corrected_Copy_Number"
        for sample in mpe_samples
        if f"{sample}.Corrected_Copy_Number" in cna.columns
    ]
    matrix = cna[sample_columns].apply(pd.to_numeric, errors="coerce")

    cna_bins = cna[["chr", "start", "end"]].copy()
    cna_bins["Gain_Freq"] = (matrix > 2.01).mean(axis=1)
    cna_bins["Loss_Freq"] = (matrix < 1.99).mean(axis=1)

    arm_source = cna[["chr", "start", "end"]].copy()
    arm_source["chr"] = arm_source["chr"].astype(str).str.strip()
    arm_source["chr"] = arm_source["chr"].where(
        arm_source["chr"].str.startswith("chr"), "chr" + arm_source["chr"]
    )
    arm_source = arm_source[arm_source["chr"].isin(CHROMOSOMES)].copy()
    arm_source["Arm"] = arm_source.apply(
        lambda row: chromosome_arm(row["chr"], row["start"], row["end"]), axis=1
    )

    is_gain = matrix.loc[arm_source.index] > 2.1
    is_loss = matrix.loc[arm_source.index] < 1.9
    arm_statistics = {}

    for arm, arm_df in arm_source.groupby("Arm"):
        gain_samples = is_gain.loc[arm_df.index].mean(axis=0) > 0.8
        loss_samples = is_loss.loc[arm_df.index].mean(axis=0) > 0.8
        arm_statistics[arm] = {
            "Gain_Freq(%)": gain_samples.mean() * 100,
            "Loss_Freq(%)": loss_samples.mean() * 100,
        }

    return cna_bins, pd.DataFrame.from_dict(arm_statistics, orient="index")


def read_bed(filename):
    bed = pd.read_csv(filename, sep="\t", header=None, comment="#", usecols=[0, 1, 2])
    bed.columns = ["chr", "start", "end"]
    bed["chr"] = bed["chr"].astype(str)
    bed["start"] = pd.to_numeric(bed["start"], errors="coerce")
    bed["end"] = pd.to_numeric(bed["end"], errors="coerce")
    bed = bed.dropna(subset=["start", "end"])
    bed[["start", "end"]] = bed[["start", "end"]].astype(int)
    return bed[bed["chr"].isin(CHROMOSOMES)].copy()


def mark_overlap(regions, bed):
    overlap = np.zeros(len(regions), dtype=bool)

    for chromosome in CHROMOSOMES:
        region_indices = regions.index[regions["chr"] == chromosome].to_numpy()
        chromosome_bed = bed[bed["chr"] == chromosome].sort_values("start")
        bed_starts = chromosome_bed["start"].to_numpy()
        bed_ends = chromosome_bed["end"].to_numpy()

        for index in region_indices:
            candidate_count = np.searchsorted(bed_starts, regions.at[index, "end"], side="left")
            overlap[index] = np.any(bed_ends[:candidate_count] > regions.at[index, "start"])

    return overlap


def calculate_dmr(mpe_samples, bpe_samples):
    dmr = pd.read_csv(DMR_FILE, sep="\t")
    valid_mpe = [sample for sample in mpe_samples if sample in dmr.columns]
    valid_bpe = [sample for sample in bpe_samples if sample in dmr.columns]

    positions = dmr.iloc[:, 0].astype(str).str.split("_", expand=True).iloc[:, :3]
    positions.columns = ["chr", "start", "end"]
    positions["start"] = pd.to_numeric(positions["start"], errors="coerce")
    positions["end"] = pd.to_numeric(positions["end"], errors="coerce")
    keep = (
        positions["chr"].isin(CHROMOSOMES)
        & positions["start"].notna()
        & positions["end"].notna()
    )
    positions = positions.loc[keep].reset_index(drop=True)
    positions[["start", "end"]] = positions[["start", "end"]].astype(int)
    dmr = dmr.loc[keep].reset_index(drop=True)

    positions["is_hyper_region"] = mark_overlap(positions, read_bed(HYPER_BED))
    positions["is_hypo_region"] = mark_overlap(positions, read_bed(HYPO_BED))

    methylation = dmr[valid_mpe + valid_bpe].apply(pd.to_numeric, errors="coerce")
    bpe_matrix = methylation[valid_bpe]
    mpe_matrix = methylation[valid_mpe].to_numpy(dtype=float)
    mpe_valid = ~np.isnan(mpe_matrix)
    mpe_valid_count = mpe_valid.sum(axis=1)
    bpe_q05 = bpe_matrix.quantile(0.05, axis=1, interpolation="linear").to_numpy()
    bpe_q95 = bpe_matrix.quantile(0.95, axis=1, interpolation="linear").to_numpy()

    hyper_change = (
        (mpe_matrix >= bpe_q95[:, None])
        & mpe_valid
        & positions["is_hyper_region"].to_numpy()[:, None]
    )
    hypo_change = (
        (mpe_matrix <= bpe_q05[:, None])
        & mpe_valid
        & positions["is_hypo_region"].to_numpy()[:, None]
    )
    hyper_frequency = np.divide(
        hyper_change.sum(axis=1),
        mpe_valid_count,
        out=np.full(len(positions), np.nan),
        where=mpe_valid_count > 0,
    )
    hypo_frequency = np.divide(
        hypo_change.sum(axis=1),
        mpe_valid_count,
        out=np.full(len(positions), np.nan),
        where=mpe_valid_count > 0,
    )

    positions["hyper_freq"] = np.where(
        positions["is_hyper_region"], hyper_frequency, np.nan
    )
    positions["hypo_freq"] = np.where(
        positions["is_hypo_region"], hypo_frequency, np.nan
    )
    return positions[
        (positions["is_hyper_region"] | positions["is_hypo_region"])
        & pd.Series(bpe_q05).notna()
        & pd.Series(bpe_q95).notna()
        & (mpe_valid_count > 0)
    ].copy()


def parse_region(region):
    chromosome, start, end = str(region).split("_")[:3]
    return chromosome, int(start), int(end)


def reduction_frequency(matrix, mpe_samples, bpe_samples):
    valid_mpe = [sample for sample in mpe_samples if sample in matrix.index]
    valid_bpe = [sample for sample in bpe_samples if sample in matrix.index]
    mpe_matrix = matrix.loc[valid_mpe]
    bpe_matrix = matrix.loc[valid_bpe]
    threshold = bpe_matrix.quantile(0.05, axis=0, interpolation="linear")
    valid = mpe_matrix.notna()
    abnormal = mpe_matrix.le(threshold, axis=1) & valid
    return abnormal.sum(axis=0).divide(valid.sum(axis=0).replace(0, np.nan))


def calculate_length(mpe_samples, bpe_samples):
    matrix = pd.read_csv(LENGTH_MATRIX_FILE, index_col=0)
    matrix.index = matrix.index.astype(str).str.strip()
    feature_columns = [
        column
        for column in matrix.columns
        if len(str(column).split("_")) >= 3
        and str(column).split("_")[0] in CHROMOSOMES
    ]
    matrix = matrix[feature_columns].apply(pd.to_numeric, errors="coerce")
    frequency = reduction_frequency(matrix, mpe_samples, bpe_samples)

    rows = []
    for feature in feature_columns:
        chromosome, start, end = parse_region(feature)
        rows.append({
            "chr": chromosome,
            "start": start,
            "end": end,
            "MPE_abnormal_freq": frequency[feature],
        })
    return pd.DataFrame(rows)


def calculate_c_motif(mpe_samples, bpe_samples):
    matrix = pd.read_csv(C_MOTIF_MATRIX_FILE, index_col=0)
    matrix.index = matrix.index.astype(str).str.strip()
    valid_samples = [
        sample for sample in mpe_samples + bpe_samples if sample in matrix.index
    ]
    matrix = matrix.loc[valid_samples].copy()
    feature_columns = [str(column) for column in matrix.columns]
    regions = sorted({
        column.rsplit("_", 1)[0]
        for column in feature_columns
        if "_" in column
        and len(column.rsplit("_", 1)[-1]) == 4
        and column.split("_")[0] in CHROMOSOMES
    })

    c_matrix = pd.DataFrame(index=matrix.index)
    for region in regions:
        motif_columns = {
            column.rsplit("_", 1)[-1].upper(): column
            for column in feature_columns
            if column.startswith(f"{region}_")
            and column.rsplit("_", 1)[-1].upper() in EXPECTED_C_4MERS
        }
        if set(motif_columns) == EXPECTED_C_4MERS:
            selected = [motif_columns[motif] for motif in sorted(EXPECTED_C_4MERS)]
            c_matrix[region] = matrix[selected].apply(
                pd.to_numeric, errors="coerce"
            ).sum(axis=1, min_count=64)

    frequency = reduction_frequency(c_matrix, mpe_samples, bpe_samples)
    rows = []
    for region in c_matrix.columns:
        chromosome, start, end = parse_region(region)
        rows.append({
            "chr": chromosome,
            "start": start,
            "end": end,
            "MPE_abnormal_freq": frequency[region],
        })
    return pd.DataFrame(rows)


def normalise_chromosome(values):
    values = values.astype(str).str.strip()
    return values.where(values.str.startswith("chr"), "chr" + values)


def load_genome_layout():
    offsets = {}
    centers = []
    spans = []
    boundaries = []
    current_position = 0.0

    for chromosome in CHROMOSOMES:
        chromosome_end = current_position + CHROMOSOME_LENGTHS[chromosome]
        offsets[chromosome] = current_position
        centers.append((current_position + chromosome_end) / 2)
        spans.append((current_position, chromosome_end))
        boundaries.append(chromosome_end)
        current_position = chromosome_end + CHROMOSOME_GAP

    return CHROMOSOME_LENGTHS, offsets, centers, spans, boundaries, boundaries[-1]


def add_plot_positions(table, offsets):
    table = table.copy()
    table["chr"] = normalise_chromosome(table["chr"])
    table["start"] = pd.to_numeric(table["start"], errors="coerce")
    table["end"] = pd.to_numeric(table["end"], errors="coerce")
    table = table[
        table["chr"].isin(CHROMOSOMES)
        & table["start"].notna()
        & table["end"].notna()
    ].copy()
    table["x_pos"] = [
        offsets[chromosome] + (start + end) / 2
        for chromosome, start, end in zip(table["chr"], table["start"], table["end"])
    ]
    table["plot_width"] = np.maximum(table["end"] - table["start"], 1)
    return table


def draw_binned_track(ax, table, values, color, baseline=0.0):
    plot_data = table[["chr", "start", "end", "x_pos", "plot_width"]].copy()
    plot_data["value"] = pd.to_numeric(values, errors="coerce").to_numpy()

    for chromosome in CHROMOSOMES:
        chromosome_data = plot_data[plot_data["chr"] == chromosome].sort_values(
            ["start", "end"]
        )
        edges = []
        run_values = []
        previous_right = None

        for row in chromosome_data.itertuples(index=False):
            left = float(row.x_pos - row.plot_width / 2)
            right = float(row.x_pos + row.plot_width / 2)
            contiguous = previous_right is not None and abs(left - previous_right) <= 10

            if run_values and (not np.isfinite(row.value) or not contiguous):
                ax.stairs(
                    run_values, edges, baseline=baseline, fill=True,
                    facecolor=color, edgecolor="none", linewidth=0, antialiased=False,
                )
                edges, run_values, previous_right = [], [], None

            if np.isfinite(row.value):
                if not run_values:
                    edges = [left]
                else:
                    left = previous_right
                run_values.append(float(row.value))
                edges.append(right)
                previous_right = right

        if run_values:
            ax.stairs(
                run_values, edges, baseline=baseline, fill=True,
                facecolor=color, edgecolor="none", linewidth=0, antialiased=False,
            )


def add_genome_background(ax, spans, boundaries, x_end):
    for index, (start, end) in enumerate(spans):
        if index % 2 == 1:
            ax.axvspan(start, end, color=COLORS["alternate"], linewidth=0, zorder=0)
    for boundary in boundaries[:-1]:
        ax.axvline(boundary, color=COLORS["boundary"], linewidth=0.35, zorder=1)
    ax.set_xlim(0, x_end)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", width=0.60, length=2.2, labelsize=5.8)


def format_diverging_axis(ax):
    ax.set_ylim(-1, 1)
    ax.set_yticks([-1, -0.5, 0, 0.5, 1])
    ax.set_yticklabels(["100", "50", "0", "50", "100"])
    ax.axhline(0, color="black", linewidth=0.65, zorder=5)


def arm_span(arm, chromosome_lengths, offsets):
    chromosome = arm[:-1]
    if arm.endswith("p"):
        start, end = 0, CENTROMERES[chromosome][0]
    else:
        start, end = CENTROMERES[chromosome][1], chromosome_lengths[chromosome]
    return offsets[chromosome] + start, offsets[chromosome] + end


def plot_cna(ax, cna, arm_frequency, chromosome_lengths, offsets):
    draw_binned_track(ax, cna, cna["Gain_Freq"], COLORS["gain"])
    draw_binned_track(ax, cna, -cna["Loss_Freq"], COLORS["loss"])
    format_diverging_axis(ax)

    for arm, row in arm_frequency.iterrows():
        start, end = arm_span(arm, chromosome_lengths, offsets)
        gain = row["Gain_Freq(%)"] / 100
        loss = row["Loss_Freq(%)"] / 100
        label = arm.replace("chr", "")

        if gain >= 0.30:
            ax.plot([start, end], [0.94, 0.94], color=COLORS["gain"], lw=0.8)
            ax.text((start + end) / 2, 0.89, label, ha="center", va="top", fontsize=5)
        if loss >= 0.30:
            ax.plot([start, end], [-0.94, -0.94], color=COLORS["loss"], lw=0.8)
            ax.text((start + end) / 2, -0.89, label, ha="center", va="bottom", fontsize=5)

    ax.set_ylabel("MPE samples with CNA\ngain or loss per bin (%)", fontsize=6.2, labelpad=7)
    ax.legend(
        handles=[
            Patch(facecolor=COLORS["gain"], label="Bin-level gain"),
            Patch(facecolor=COLORS["loss"], label="Bin-level loss"),
            Line2D([0], [0], color=COLORS["gain"], lw=1, label="Recurrent arm-level gain (>=30%)"),
            Line2D([0], [0], color=COLORS["loss"], lw=1, label="Recurrent arm-level loss (>=30%)"),
        ],
        loc="upper right", ncol=2, fontsize=5, frameon=False,
        handlelength=1.3, columnspacing=0.8, labelspacing=0.25,
    )


def plot_dmr(ax, dmr):
    hyper = dmr[dmr["is_hyper_region"] & dmr["hyper_freq"].notna()]
    hypo = dmr[dmr["is_hypo_region"] & dmr["hypo_freq"].notna()]
    ax.scatter(hyper["x_pos"], hyper["hyper_freq"], s=2.5, color=COLORS["hyper"], edgecolors="none")
    ax.scatter(hypo["x_pos"], -hypo["hypo_freq"], s=2.5, color=COLORS["hypo"], edgecolors="none")
    format_diverging_axis(ax)
    ax.set_ylabel("MPE samples with abnormal\nmethylation per DMR (%)", fontsize=6.2, labelpad=7)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", linestyle="none", markersize=3.2,
                   markerfacecolor=COLORS["hyper"], markeredgewidth=0,
                   label="Hyper-DMR (>= BPE 95th percentile)"),
            Line2D([0], [0], marker="o", linestyle="none", markersize=3.2,
                   markerfacecolor=COLORS["hypo"], markeredgewidth=0,
                   label="Hypo-DMR (<= BPE 5th percentile)"),
        ],
        loc="upper right", ncol=2, fontsize=5, frameon=False,
        handletextpad=0.35, columnspacing=0.8,
    )


def plot_fragment_track(ax, table, color, ylabel, label, y_max):
    draw_binned_track(ax, table, table["MPE_abnormal_freq"], color)
    ax.set_ylim(0, y_max)
    step = 0.10 if y_max <= 0.50 else 0.20
    ticks = np.arange(0, y_max + step / 2, step)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{value * 100:.0f}" for value in ticks])
    ax.axhline(0, color="black", linewidth=0.65)
    ax.set_ylabel(ylabel, fontsize=6.2, labelpad=7)
    ax.legend(
        handles=[Patch(facecolor=color, edgecolor="none", label=label)],
        loc="upper right", fontsize=5, frameon=False,
    )


def main():
    os.makedirs(os.path.dirname(OUTPUT_PDF), exist_ok=True)
    mpe_samples, bpe_samples = load_training_samples()

    cna, arm_frequency = calculate_cna(mpe_samples)
    dmr = calculate_dmr(mpe_samples, bpe_samples)
    length = calculate_length(mpe_samples, bpe_samples)
    c_motif = calculate_c_motif(mpe_samples, bpe_samples)

    chromosome_lengths, offsets, centers, spans, boundaries, x_end = load_genome_layout()
    cna = add_plot_positions(cna, offsets)
    dmr = add_plot_positions(dmr, offsets)
    length = add_plot_positions(length, offsets)
    c_motif = add_plot_positions(c_motif, offsets)

    fragment_values = pd.concat(
        [length["MPE_abnormal_freq"], c_motif["MPE_abnormal_freq"]]
    ).dropna()
    fragment_y_max = min(1, max(0.30, np.ceil(fragment_values.max() * 1.08 / 0.05) * 0.05))

    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 5.8,
        "axes.linewidth": 0.65,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "pdf.compression": 6,
    })

    fig, axes = plt.subplots(
        4, 1, sharex=True, figsize=(176 / 25.4, 126 / 25.4),
        gridspec_kw={"height_ratios": [1.30, 1.30, 1.00, 1.00]},
        facecolor="white",
    )

    for ax in axes:
        ax.set_facecolor("white")
        add_genome_background(ax, spans, boundaries, x_end)

    plot_cna(axes[0], cna, arm_frequency, chromosome_lengths, offsets)
    plot_dmr(axes[1], dmr)
    plot_fragment_track(
        axes[2], c_motif, COLORS["cem"],
        "MPE samples with reduced\nC-EM frequency per bin (%)",
        "Reduced C-EM frequency", fragment_y_max,
    )
    plot_fragment_track(
        axes[3], length, COLORS["length"],
        "MPE samples with reduced\n220-270-bp fraction per bin (%)",
        "Reduced 220-270-bp fragment fraction", fragment_y_max,
    )

    for ax in axes[:-1]:
        ax.tick_params(axis="x", which="both", bottom=False, labelbottom=False)

    axes[-1].set_xticks(centers)
    axes[-1].set_xticklabels([chromosome.replace("chr", "") for chromosome in CHROMOSOMES], fontsize=5.6)
    axes[-1].set_xlabel("Chromosome", fontsize=6.2, labelpad=2.5)

    fig.align_ylabels(axes)
    fig.subplots_adjust(left=0.145, right=0.995, bottom=0.072, top=0.985, hspace=0.20)
    fig.savefig(
        OUTPUT_PDF,
        format="pdf",
        facecolor="white",
        transparent=False,
        bbox_inches=None,
        pad_inches=0,
    )
    plt.close(fig)


if __name__ == "__main__":
    main()
