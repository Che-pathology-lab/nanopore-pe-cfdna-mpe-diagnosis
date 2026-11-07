"""Purpose: Merge cna segments for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import glob
import os

INPUT_DIR = "results/03_nanopore_cna/3.3_ichorCNA"
OUTPUT_DIR = "results/03_nanopore_cna"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, 'discovery_Corrected_Copy_Number.csv')

def main():
    search_pattern = os.path.join(INPUT_DIR, "*.cna.seg")
    files = glob.glob(search_pattern)
    dfs = []
    for file in files:
        df = pd.read_csv(file, sep='\t')
        target_cols = [col for col in df.columns if 'Corrected_Copy_Number' in col]
        target_col = target_cols[0]
        df_sub = df[['chr', 'start', 'end', target_col]].copy()
        df_sub.set_index(['chr', 'start', 'end'], inplace=True)
        dfs.append(df_sub)
    merged_df = pd.concat(dfs, axis=1).reset_index()
    chrom_order = [f"chr{i}" for i in range(1, 23)]
    merged_df['chr'] = pd.Categorical(merged_df['chr'], categories=chrom_order, ordered=True)
    merged_df = merged_df.sort_values(by=['chr', 'start', 'end'])
    merged_df = merged_df.dropna(subset=['chr'])
    merged_df = merged_df.reset_index(drop=True)
    merged_df.to_csv(OUTPUT_FILE, index=False)


if __name__ == "__main__":
    main()
