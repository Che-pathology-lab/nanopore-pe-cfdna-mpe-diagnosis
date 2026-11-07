"""Purpose: Ngs multimodal feature merge for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd

CNA_FILE = "results/07_ngs_validation/CNA/NGS_CNA_features_wide_table.csv"
METH_FILE = "results/07_ngs_validation/METH/NGS_validation_global_methylation_matrix.csv"
FRAG_FILE = "results/07_ngs_validation/FRAG/NGS_FRAG_features_wide_table.csv"
OUTPUT_FILE = "results/07_ngs_validation/Multi/NGS_Multi_features_wide_table.csv"
def main():
    
    df_cnv = pd.read_csv(CNA_FILE)
    df_meth = pd.read_csv(METH_FILE)
    df_frag = pd.read_csv(FRAG_FILE)

    for df in [df_cnv, df_meth, df_frag]:
        df['Sample'] = df['Sample'].astype(str).str.strip()

    df_merge = pd.merge(df_cnv, df_frag, on='Sample', how='inner')
    df_final = pd.merge(df_merge, df_meth, on='Sample', how='inner')

    df_final.sort_values(by='Sample', inplace=True)
    df_final.reset_index(drop=True, inplace=True)
    df_final.to_csv(OUTPUT_FILE, index=False)


if __name__ == "__main__":
    main()
