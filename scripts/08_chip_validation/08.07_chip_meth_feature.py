"""Purpose: Chip meth feature for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import numpy as np

NGS_METH_FILE = "results/08_chip_validation/CHIP_METH/CHIP_DMR.tsv" 
NGS_GROUP_FILE = "metadata/chip_group.csv"
HYPER_BED = "results/02_nanopore_methylation/hyper_DMR.bed"
HYPO_BED = "results/02_nanopore_methylation/hypo_DMR.bed"
OUTPUT_MATRIX = "results/08_chip_validation/CHIP_METH/CHIP_validation_global_methylation_matrix.csv"


def load_bed_as_index(bed_file):
    bed_df = pd.read_csv(bed_file, sep='\t', header=None, comment='#')
    index_list = bed_df[0].astype(str) + '_' + bed_df[1].astype(str) + '_' + bed_df[2].astype(str)
    return index_list.tolist()

def load_and_format_ngs_data():
    meth_df = pd.read_csv(NGS_METH_FILE, sep='\t')

    meth_df['chr_start_end'] = meth_df.iloc[:, 0].astype(str) + '_' + \
                               meth_df.iloc[:, 1].astype(str) + '_' + \
                               meth_df.iloc[:, 2].astype(str)
    
    meth_df = meth_df.set_index('chr_start_end')
    meth_df = meth_df.drop(meth_df.columns[0:3], axis=1)
    
    available_samples = meth_df.columns.tolist()
    group_df = pd.read_csv(NGS_GROUP_FILE)
    target_samples = [s for s in group_df['Sample'].tolist() if s in available_samples]
    
    return meth_df, target_samples

def calculate_global_matrix(df, hyper_ids, hypo_ids, target_samples):
    result_matrix = pd.DataFrame(index=target_samples)
    result_matrix.index.name = 'Sample'
    valid_hyper = [i for i in hyper_ids if i in df.index]
    valid_hypo = [i for i in hypo_ids if i in df.index]
    print(f" valid Hyper DMRs: {len(valid_hyper)} / {len(hyper_ids)}")
    print(f" valid Hypo DMRs: {len(valid_hypo)} / {len(hypo_ids)}")

    if valid_hyper:
        result_matrix['Global_Hyper_Mean'] = df.loc[valid_hyper, target_samples].mean(axis=0, skipna=True)
    else:
        result_matrix['Global_Hyper_Mean'] = np.nan
        
    if valid_hypo:
        result_matrix['Global_Hypo_Mean'] = df.loc[valid_hypo, target_samples].mean(axis=0, skipna=True)
    else:
        result_matrix['Global_Hypo_Mean'] = np.nan
    return result_matrix

def main() -> None:
    meth_df, valid_samples = load_and_format_ngs_data()
    hyper_list = load_bed_as_index(HYPER_BED)
    hypo_list = load_bed_as_index(HYPO_BED)
    final_matrix = calculate_global_matrix(meth_df, hyper_list, hypo_list, valid_samples)
    final_matrix.to_csv(OUTPUT_MATRIX)


if __name__ == "__main__":
    main()
