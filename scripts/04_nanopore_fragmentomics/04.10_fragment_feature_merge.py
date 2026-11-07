"""Purpose: Fragment feature merge for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import os
from functools import reduce


INPUT_FILES = {
    '1mer': "results/04_nanopore_fragmentomics/motif_1mer_significant.csv",
    '2mer': "results/04_nanopore_fragmentomics/motif_2mer_significant.csv",
    '3mer': "results/04_nanopore_fragmentomics/motif_3mer_significant.csv",
    '4mer': "results/04_nanopore_fragmentomics/motif_4mer_significant.csv",
    'entropy': "results/04_nanopore_fragmentomics/motif_entropy_matrix.csv",
    'size': "results/04_nanopore_fragmentomics/discovery_size_feature_significant.csv"
}
GROUP_FILE = "metadata/discovery_groupfile.csv"
OUTPUT_DIR = "results/04_nanopore_fragmentomics"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "discovery_FRAG_features_significant.csv")

def main():


    df_group = pd.read_csv(GROUP_FILE)
    target_samples = df_group['Sample'].astype(str).str.strip().unique()

    df_list = []
    for filepath in INPUT_FILES.values():
        df = pd.read_csv(filepath)          
        df['Sample'] = df['Sample'].astype(str).str.strip()
        df_list.append(df)
    merged_df = reduce(lambda left, right: pd.merge(left, right, on='Sample', how='outer'), df_list)
    merged_df = merged_df[merged_df['Sample'].isin(target_samples)].copy()
    merged_df.sort_values(by='Sample', inplace=True)
    merged_df.reset_index(drop=True, inplace=True)

    merged_df.to_csv(OUTPUT_FILE, index=False)


if __name__ == "__main__":
    main()
