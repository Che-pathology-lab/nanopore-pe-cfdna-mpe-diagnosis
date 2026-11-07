"""Purpose: Nanopore for chip multimodal feature merge for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd

# input: 2组学 vif筛选后的特征矩阵
# output: 2组学合并的筛选后的特征矩阵

CNA_FILE = "results/08_chip_validation/ONT_for_CHIP_CNA/discovery_ONT_for_CHIP_CNA_features_vif_filtered.csv"
METH_FILE = "results/02_nanopore_methylation/discovery_METH_features_vif_filtered.csv"
GROUP_INFO_FILE = "metadata/discovery_groupfile.csv"
CLEAN_FEATURE_TABLE = "results/08_chip_validation/ONT_for_CHIP_Multi/ONT_for_CHIP_Multi_features_wide_table.csv"

def main():
    df_cnv = pd.read_csv(CNA_FILE)
    df_meth = pd.read_csv(METH_FILE)
    group_info_df = pd.read_csv(GROUP_INFO_FILE)

    for df in [df_cnv, df_meth, group_info_df]:
        df['Sample'] = df['Sample'].astype(str).str.strip()
    df_merge = pd.merge(df_cnv, df_meth, on='Sample', how='inner')

    cols_to_use = ['Sample', 'group', 'split']
        
    df_full = pd.merge(df_merge, group_info_df[cols_to_use], on='Sample', how='inner')
    df_full = df_full[df_full['group'].isin(['mpe', 'bpe'])].copy()
    
    df_full.sort_values(by='Sample', inplace=True)
    df_full.reset_index(drop=True, inplace=True)
    
    exclude_cols = ['Sample', 'group', 'split', 'subgroup']
    FEATURE_COLS = [col for col in df_full.columns if col not in exclude_cols]

    columns_to_save = ['Sample'] + FEATURE_COLS
    df_final_clean = df_full[columns_to_save]

    df_final_clean.to_csv(CLEAN_FEATURE_TABLE, index=False)


if __name__ == "__main__":
    main()
