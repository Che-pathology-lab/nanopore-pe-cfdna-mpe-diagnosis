"""Purpose: Glm dmr sigtest for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings("ignore")

FILE_PATH = "results/02_nanopore_methylation/threshold_10_discovery_cohort_methylationfile.tsv.gz"
GROUP_INFO_FILE = "metadata/discovery_groupfile.csv"
NA_PERCENT = 0.1
AGE_THRESHOLD = 60
Q_VALUE = 0.1
MEAN_DIFF = 15
OUTPUT_RAW_DATA = f"results/02_nanopore_methylation/GLM_filtered_METH_mean_diff_{MEAN_DIFF}_q_{Q_VALUE}_raw_data.tsv"
OUTPUT_STATS = "results/02_nanopore_methylation/GLM_statistics_results.tsv"


def main() -> None:
    df = pd.read_csv(FILE_PATH, sep='\t')
    valid_chroms = [f'chr{i}' for i in range(1, 23)]
    chrom_col_name = df.columns[0] 
    df = df[df[chrom_col_name].isin(valid_chroms)]

    original_columns = df.columns.tolist()
    df['chr_start_end'] = df.iloc[:, :3].astype(str).apply('_'.join, axis=1)
    cols = ['chr_start_end'] + [col for col in df.columns if col not in original_columns[:3] and col != 'chr_start_end']
    df = df[cols]

    data_columns = [col for col in df.columns if col != 'chr_start_end']
    row_na_counts_data = df[data_columns].isna().sum(axis=1)
    df_cleaned = df[row_na_counts_data / len(data_columns) <= NA_PERCENT]


    group_info_df = pd.read_csv(GROUP_INFO_FILE)

    group_info_df['Sample'] = group_info_df['Sample'].astype(str).str.strip()
    group_info_df['split'] = group_info_df['split'].astype(str).str.strip().str.lower()
    group_info_df['split'] = group_info_df['split'].replace({
        'panelofnormal': 'train',
    })
    group_info_df['group'] = group_info_df['group'].astype(str).str.strip().str.lower()

    train_samples_df = group_info_df[
        group_info_df['split'] == 'train'
    ].copy()

    train_samples_df = train_samples_df[
        train_samples_df['group'].isin(['mpe', 'bpe'])
    ].copy()
    train_samples_df['group_bin'] = train_samples_df['group'].map({'mpe': 1, 'bpe': 0})
    train_samples_df['gender_bin'] = pd.to_numeric(train_samples_df['gender'], errors='coerce')
    train_samples_df['age'] = pd.to_numeric(train_samples_df['age'], errors='coerce')
    train_samples_df = train_samples_df.dropna(subset=['Sample', 'group_bin', 'gender_bin', 'age'])
    train_samples_df['age_group'] = (train_samples_df['age'] >= AGE_THRESHOLD).astype(int)

    age_dist = train_samples_df['age_group'].value_counts()
    available_columns_in_df = set(df_cleaned.columns)
    train_samples_df = train_samples_df[train_samples_df['Sample'].isin(available_columns_in_df)]
    sample_list = train_samples_df['Sample'].tolist()
    y_group = train_samples_df['group_bin'].values
    covariates_df = train_samples_df[['age_group', 'gender_bin']].copy()

    test_results = []
    processed_count = 0
    total_rows = len(df_cleaned)
    for index, row in df_cleaned.iterrows():
        meth_values = row[sample_list].astype(float).values

        temp_df = pd.DataFrame({
            'group': y_group,
            'meth': meth_values,
            'age': covariates_df['age_group'].values,
            'gender': covariates_df['gender_bin'].values
        })
    
        temp_df = temp_df.dropna()

        group_counts = temp_df['group'].value_counts()
        if len(group_counts) < 2 or group_counts.min() < 5:
            continue

        ca_values = temp_df[temp_df['group'] == 1]['meth'].values
        be_values = temp_df[temp_df['group'] == 0]['meth'].values
        mean_ca = np.mean(ca_values)
        mean_be = np.mean(be_values)
        mean_diff = abs(mean_ca - mean_be)

        X = temp_df[['meth', 'age', 'gender']]
        X = sm.add_constant(X)
        y = temp_df['group']
    
        model = sm.GLM(y, X, family=sm.families.Binomial()).fit()

        pvalue = model.pvalues['meth']
        coef = model.params['meth']
        bse = model.bse['meth']
        
        test_results.append({
            'CpG_Index': index, 
            'chr_start_end': row['chr_start_end'],
            'P-value': pvalue,
            'GLM_Coefficient': coef,
            'Mean_CA': mean_ca,
            'Mean_BE': mean_be,
            'Mean_Difference': mean_diff
        })
        processed_count += 1
    results_df = pd.DataFrame(test_results)
    results_df = results_df.dropna(subset=['P-value'])
    raw_pvalues = results_df['P-value'].values
    rejected, qvalues, _, _ = multipletests(raw_pvalues, alpha=0.1, method='fdr_bh')
    results_df['q-value'] = qvalues

    significantBH_cpg_df = results_df[results_df['q-value'] < Q_VALUE]
    filtered_diff_df_mean = significantBH_cpg_df[significantBH_cpg_df['Mean_Difference'] >= MEAN_DIFF]
    mean_cpg_indices = filtered_diff_df_mean['CpG_Index'].tolist()
    filtered_cpg_raw_data = df_cleaned.loc[mean_cpg_indices].copy()

    other_columns = [col for col in filtered_cpg_raw_data.columns if col != 'chr_start_end']
    final_output = filtered_cpg_raw_data[['chr_start_end'] + other_columns]
    final_output.to_csv(OUTPUT_RAW_DATA, sep='\t', index=False, na_rep='NA')
    stat_output = results_df.loc[results_df['CpG_Index'].isin(mean_cpg_indices)].copy()
    stat_output.to_csv(OUTPUT_STATS, sep='\t', index=False)



if __name__ == "__main__":
    main()
