"""Purpose: Methylation feature for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd


METH_FILE = "results/02_nanopore_methylation/GLM_filtered_METH_mean_diff_15_q_0.1_raw_data.tsv"
GROUP_FILE = "metadata/discovery_groupfile.csv"
OUTPUT_HYPER_BED = "results/02_nanopore_methylation/hyper_DMR.bed"
OUTPUT_HYPO_BED = "results/02_nanopore_methylation/hypo_DMR.bed"
OUTPUT_MATRIX = "results/02_nanopore_methylation/nanopore_discovery_global_methylation_matrix.csv"

def load_and_map_samples():
    group_df = pd.read_csv(GROUP_FILE)
    meth_df = pd.read_csv(METH_FILE, sep='\t')
    feature_col = meth_df.columns[0]
    meth_df = meth_df.set_index(feature_col)
    available_samples = meth_df.columns.tolist()
    train_df = group_df[group_df['split'] == 'train']
    train_ca_samples = train_df[train_df['group'] == 'mpe']['Sample'].tolist()
    train_be_samples = train_df[train_df['group'] == 'bpe']['Sample'].tolist()
    valid_train_ca = [s for s in train_ca_samples if s in available_samples]
    valid_train_be = [s for s in train_be_samples if s in available_samples]
    all_target_samples = [s for s in group_df['Sample'].tolist() if s in available_samples]
    
    return meth_df, valid_train_ca, valid_train_be, all_target_samples

def filter_dmrs(df, ca_cols, be_cols):
    mean_ca = df[ca_cols].mean(axis=1)
    mean_be = df[be_cols].mean(axis=1)
    diff = mean_ca - mean_be
    mask_hyper = (diff > 0)
    hyper_dmrs = df.index[mask_hyper].tolist()
    mask_hypo = (diff < 0)
    hypo_dmrs = df.index[mask_hypo].tolist()
    return hyper_dmrs, hypo_dmrs

def save_bed(dmr_list, filename):
    coords = pd.Series(dmr_list).str.split('_', expand=True)
    coords = coords.iloc[:, :3]
    coords.to_csv(filename, sep='\t', header=False, index=False)

def calculate_global_matrix(df, hyper_ids, hypo_ids, target_samples):
    
    result_matrix = pd.DataFrame(index=target_samples)
    if hyper_ids:
        result_matrix['Global_Hyper_Mean'] = df.loc[hyper_ids, target_samples].mean(axis=0, skipna=True)
    else:
        result_matrix['Global_Hyper_Mean'] = 0
        
    if hypo_ids:
        result_matrix['Global_Hypo_Mean'] = df.loc[hypo_ids, target_samples].mean(axis=0, skipna=True)
    else:
        result_matrix['Global_Hypo_Mean'] = 0
        
    return result_matrix
def main() -> None:
    meth_df, train_ca, train_be, all_samples = load_and_map_samples()
    hyper_list, hypo_list = filter_dmrs(meth_df, train_ca, train_be)
    save_bed(hyper_list, OUTPUT_HYPER_BED)
    save_bed(hypo_list, OUTPUT_HYPO_BED)
    final_matrix = calculate_global_matrix(meth_df, hyper_list, hypo_list, all_samples)
    final_matrix.to_csv(OUTPUT_MATRIX, index_label='Sample')


if __name__ == "__main__":
    main()
