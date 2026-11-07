"""Purpose: Downsample fastq for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import os


COV_FILE = "results/01_nanopore_preprocessing/1.5_coverage/coverage.csv"  # 1.5 mosdepth output; requires Sample and MeanCov.
GROUP_FILE = "metadata/discovery_groupfile.csv"
INPUT_FASTQ_DIR = "results/01_nanopore_preprocessing/1.3_fastq"
OUTPUT_FASTQ_DIR = "results/06_downsampling/1.3_fastq"

TARGET_DEPTHS = (0.01, 0.1)
NUM_REPLICATES = 10
BASE_SEED = 42


def load_test_coverage():
    df_cov = pd.read_csv(COV_FILE)
    df_group = pd.read_csv(GROUP_FILE)

    df_cov["Sample"] = df_cov["Sample"].astype(str).str.strip()
    df_group["Sample"] = df_group["Sample"].astype(str).str.strip()
    df_group["split"] = df_group["split"].astype(str).str.strip().str.lower()

    test_samples = set(
        df_group.loc[df_group["split"] == "test", "Sample"]
    )
    if not test_samples:
        raise ValueError(
            f"No samples with split='test' were found in {GROUP_FILE}"
        )

    df_cov = df_cov[df_cov["Sample"].isin(test_samples)].copy()
    missing_samples = sorted(test_samples - set(df_cov["Sample"]))
    if missing_samples:
        raise ValueError(
            f"Test samples missing from {COV_FILE}: {missing_samples}"
        )

    if df_cov["Sample"].duplicated().any():
        duplicates = sorted(
            df_cov.loc[df_cov["Sample"].duplicated(), "Sample"].unique()
        )
        raise ValueError(f"Duplicate coverage rows found: {duplicates}")

    mean_cov = pd.to_numeric(df_cov["MeanCov"], errors="coerce")
    if mean_cov.isna().any() or (mean_cov <= 0).any():
        raise ValueError("MeanCov must be a positive number for every test sample")
    df_cov["MeanCov"] = mean_cov

    return df_cov.sort_values("Sample").reset_index(drop=True)


def generate_repeated_downsample_script(df, target_depth):
    output_script = f"run_{target_depth}x_downsample.sh"
    
    shell_commands = [
        "#!/bin/bash",
        f"mkdir -p {OUTPUT_FASTQ_DIR}"
    ]
    
    for _, row in df.iterrows():
        sample = str(row['Sample']).strip()
        current_cov = float(row['MeanCov'])
        fraction = target_depth / current_cov

        input_fq = os.path.join(INPUT_FASTQ_DIR, f"{sample}.fastq.gz")
        
        shell_commands.append(f"\n# --- Sample: {sample} (Cov: {current_cov:.2f}X) ---")
        
        for rep in range(1, NUM_REPLICATES + 1):
            current_seed = BASE_SEED + rep
            output_fq_name = f"{sample}_{target_depth}x_rep{rep}.fastq.gz"
            output_fq = os.path.join(OUTPUT_FASTQ_DIR, output_fq_name)

            if fraction >= 1.0:
                cmd = f"ln -sf {input_fq} {output_fq}"
            else:
                cmd = f"seqtk sample -s {current_seed} {input_fq} {fraction:.6f} | gzip > {output_fq}"
                
            shell_commands.append(cmd)
            
    with open(output_script, "w") as f:
        f.write("\n".join(shell_commands))
        
    os.chmod(output_script, 0o755)
    print(
        f"Wrote {output_script} for {len(df)} test samples "
        f"and {NUM_REPLICATES} replicates"
    )

def main() -> None:
    df = load_test_coverage()
    for target_depth in TARGET_DEPTHS:
        generate_repeated_downsample_script(df, target_depth)


if __name__ == "__main__":
    main()
