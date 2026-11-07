"""Purpose: Roc plot for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib as mpl
from sklearn.metrics import roc_curve, roc_auc_score
mpl.rcParams['font.family'] = 'Arial'
GROUP_FILE = "metadata/ULP_groupfile.csv"
FRAG_CSV = "results/09_ulp_validation/ULP_FRAG_prediction_results_with_scores.csv"
CNA_CSV = "results/09_ulp_validation/ULP_CNA_prediction_results_with_scores.csv"
METH_CSV = "results/09_ulp_validation/ULP_METH_prediction_results_with_scores.csv"
MULTI_CSV = "results/09_ulp_validation/ULP_Multi_prediction_results_with_scores.csv"
OUTPUT_ROC = "results/11_roc_plot/ULP_Models_ROC.pdf"

def main():


    df_group = pd.read_csv(GROUP_FILE)
    df_group['Sample'] = df_group['Sample'].astype(str).str.strip()
    df_group = df_group[df_group['group'].isin(['mpe', 'bpe'])].copy()
    df_group['y_true'] = df_group['group'].map({'mpe': 1, 'bpe': 0})
    df_merged = df_group[['Sample', 'split', 'y_true']].copy()
    models_info = [
        ('FRAG', FRAG_CSV, 'FRAG_Score_Prob'),
        ('CNA', CNA_CSV, 'CNA_Score_Prob'),
        ('METH', METH_CSV, 'METH_Score_Prob'),
        ('Multi', MULTI_CSV, 'Multi_Score_Prob')
    ]

    for model_name, csv_path, prob_col in models_info:    
            
        df_temp = pd.read_csv(csv_path)
        df_temp['Sample'] = df_temp['Sample'].astype(str).str.strip()
        
        df_temp = df_temp[['Sample', prob_col]]
        df_merged = pd.merge(df_merged, df_temp, on='Sample', how='left')
    df_ulp = df_merged[df_merged['split'] == 'ulp'].copy()
    def plot_combined_roc(df, output_path):
        plt.figure(figsize=(8, 8))
        ax = plt.gca()
        
        axis_width = 2.0
        ax.spines['bottom'].set_linewidth(axis_width)
        ax.spines['left'].set_linewidth(axis_width)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.tick_params(width=axis_width, length=6, labelsize=11)
        colors = ['#0072B2', '#D55E00', '#009E73', "#CC0F0C"] 
        
        for i, (model_name, _, prob_col) in enumerate(models_info):  
            valid_df = df.dropna(subset=[prob_col, 'y_true'])
                
            y_true = valid_df['y_true']
            y_prob = valid_df[prob_col]
            
            fpr, tpr, _ = roc_curve(y_true, y_prob)
            auc_score = roc_auc_score(y_true, y_prob)

            sensitivity = tpr
            specificity = 1 - fpr
            line_width = 2.5
            
            plt.plot(specificity, sensitivity, color=colors[i], lw=line_width, 
                     label=f'{model_name} (AUC={auc_score:.3f})')

        plt.plot([1, 0], [0, 1], color='gray', lw=1.5, linestyle='--', alpha=0.7)
        target_ticks = [0.0, 0.25, 0.50, 0.75, 1.00]
        ax.set_xticks(target_ticks)
        ax.set_yticks(target_ticks)
        formatter = ticker.FormatStrFormatter('%.2f')
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)
        plt.xlim([1.05, -0.05])
        plt.ylim([-0.05, 1.05])
        
        plt.xlabel('Specificity', fontsize=14)
        plt.ylabel('Sensitivity', fontsize=14)
        legend = plt.legend(loc="lower right", fontsize=12, frameon=False)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300)
    plot_combined_roc(df_ulp, OUTPUT_ROC)


if __name__ == "__main__":
    main()
