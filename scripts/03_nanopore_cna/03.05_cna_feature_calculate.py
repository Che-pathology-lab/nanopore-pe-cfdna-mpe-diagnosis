"""Purpose: Cna feature calculate for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import numpy as np
import os
import warnings
from scipy.stats import entropy
warnings.filterwarnings('ignore')


INPUT_FILE = "results/03_nanopore_cna/discovery_Corrected_Copy_Number.csv"
OUTPUT_DIR = "results/03_nanopore_cna"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, 'discovery_CNA_features_wide_table.csv')

ICHORCNA_REGIONS = {
    'chr1': (5000001, 250000000), 'chr2': (1, 245000000), 'chr3': (1, 195000000),
    'chr4': (1, 190000000), 'chr5': (1,  175000000), 'chr6': (1, 170000000),
    'chr7': (1, 160000000), 'chr8': (15000001, 145000000), 'chr9': (1, 140000000),
    'chr10': (1, 135000000), 'chr11': (5000001, 135000000), 'chr12': (1, 135000000),
    'chr13': (20000001, 115000000), 'chr14': (20000001, 105000000), 'chr15': (35000001, 100000000),
    'chr16': (1, 90000000), 'chr17': (1, 85000000), 'chr18': (1, 80000000),
    'chr19': (1, 60000000), 'chr20': (1, 65000000), 'chr21': (15000001, 50000000), 
    'chr22': (25000001, 55000000)
}

CENTROMERES = {
    'chr1': (120000000, 150000001), 'chr2': (85000000, 100000001), 'chr3': (85000000, 95000001),
    'chr4': (45000000, 55000001), 'chr5': (45000000, 55000001), 'chr6': (55000000, 65000001),
    'chr7': (55000000, 80000001), 'chr8': (40000000, 50000001),'chr9': (35000000, 70000001),
    'chr10': (35000000, 50000001), 'chr11': (45000000, 60000001), 'chr12': (30000000, 40000001),
    'chr13': (0, 20000001), 'chr14': (0, 20000001), 'chr15': (0, 35000001), 'chr16': (15000000, 50000001),
    'chr17': (15000000, 30000001), 'chr18': (10000000, 25000001), 'chr19': (20000000, 30000001),
    'chr20': (20000000, 35000001), 'chr21': (0, 15000001), 'chr22': (0, 25000001)
}

def get_arm(chrom, start, end):
    mid = (start + end) / 2
    if chrom not in CENTROMERES: return None
    cen_start, cen_end = CENTROMERES[chrom]
    if mid < cen_start: return f"{chrom}p"
    elif mid > cen_end: return f"{chrom}q"
    else: return None

def get_strict_oscillating_lengths(seg_cns):
    if len(seg_cns) < 3: return []
    lengths = []
    i = 0
    while i < len(seg_cns) - 2:
        if seg_cns[i] == seg_cns[i+2] and seg_cns[i] != seg_cns[i+1]:
            state_a, state_b = seg_cns[i], seg_cns[i+1]
            jump_count = 1
            j = i + 3
            while j < len(seg_cns):
                expected = state_a if (j - i) % 2 == 0 else state_b
                if seg_cns[j] == expected:
                    jump_count += 1
                    j += 1
                else:
                    break
            lengths.append(jump_count)
            i = j - 1 
        else:
            i += 1
    return lengths

def extract_wide_features(df):
    
    df['arm'] = df.apply(lambda row: get_arm(row['chr'], row['start'], row['end']), axis=1)
    df = df.dropna(subset=['arm'])
    
    sample_cols = [col for col in df.columns if col not in ['chr', 'start', 'end', 'arm']]
    wide_data = []
    
    for col in sample_cols:
        sample_name = col.replace('.Corrected_Copy_Number', '')
        sample_features = {'Sample': sample_name}
        
        df[col] = pd.to_numeric(df[col], errors='coerce')
        s_df = df[['chr', 'start', 'end', 'arm', col]].dropna(subset=[col]).sort_values(by=['chr', 'start'])
        total_genome_bins = len(s_df)
        altered_bins = 0 
        lost_bins = 0
        gain_bins = 0

        aberrant_sizes_mb = []
        all_cns = []
        changepoints = []
        bp_per_arm = []
        bp_10mb = []
        osc_lengths = []
        
        for arm in s_df['arm'].unique():
            arm_df = s_df[s_df['arm'] == arm]
            if arm_df.empty: continue
            
            cns = arm_df[col].values
            ends = arm_df['end'].values
            
            all_cns.extend(cns)
            altered_bins += np.sum(cns != 2)
            lost_bins += np.sum(cns < 2)
            gain_bins += np.sum(cns > 2)
            
            segments = []       
            arm_bp_coords = []  
            current_cn = cns[0]
            current_bins = 1
            
            for i in range(1, len(cns)):
                cn = cns[i]
                if cn == current_cn:
                    current_bins += 1
                else:
                    segments.append((current_cn, current_bins))
                    arm_bp_coords.append(ends[i-1]) 
                    current_cn = cn
                    current_bins = 1
            segments.append((current_cn, current_bins))
            
            for cn, bins_count in segments:
                if cn != 2:
                    aberrant_sizes_mb.append(bins_count * 5)
            
            arm_changepoints = [abs(segments[i][0] - segments[i-1][0]) for i in range(1, len(segments))]
            changepoints.extend(arm_changepoints)
            
            bp_per_arm.append(len(arm_changepoints))
            
            chrom = arm_df['chr'].iloc[0]
            if chrom in ICHORCNA_REGIONS and arm_bp_coords:
                chr_start, chr_len = ICHORCNA_REGIONS[chrom]
                bins_10mb = np.arange(chr_start, chr_len + 10000000, 10000000)
                counts, _ = np.histogram(arm_bp_coords, bins=bins_10mb)
                bp_10mb.extend(counts)
            seg_cns = [cn for cn, bins_count in segments]
            osc_lengths.extend(get_strict_oscillating_lengths(seg_cns))
            
        sample_features['Max_Aberrant_Seg_Size'] = np.max(aberrant_sizes_mb) if aberrant_sizes_mb else 0
        sample_features['Mean_Aberrant_Seg_Size'] = np.mean(aberrant_sizes_mb) if aberrant_sizes_mb else 0
        
        sample_features['Mean_Copy_Number'] = np.mean(all_cns) if all_cns else 2.0
        sample_features['Max_Copy_Number'] = np.max(all_cns) if all_cns else 2.0
        
        sample_features['Total_Changepoints'] = len(changepoints)
        sample_features['Mean_Changepoint_Jump'] = np.mean(changepoints) if changepoints else 0
        sample_features['Max_Changepoint_Jump'] = np.max(changepoints) if changepoints else 0
        
        sample_features['Mean_BP_per_Arm'] = np.mean(bp_per_arm) if bp_per_arm else 0 
        sample_features['Max_BP_per_Arm'] = np.max(bp_per_arm) if bp_per_arm else 0
        
        sample_features['Mean_BP_per_10MB'] = np.mean(bp_10mb) if bp_10mb else 0
        sample_features['Max_BP_per_10MB'] = np.max(bp_10mb) if bp_10mb else 0
        
        sample_features['Total_Oscillating_Chains'] = len(osc_lengths)
        sample_features['Max_Oscillating_Length'] = np.max(osc_lengths) if osc_lengths else 0

        if total_genome_bins > 0:
            sample_features['Fraction_Genome_Altered'] = altered_bins / total_genome_bins
            sample_features['Fraction_Genome_Lost'] = lost_bins / total_genome_bins
            sample_features['Fraction_Genome_Gain'] = gain_bins / total_genome_bins
            
            value_counts = pd.Series(all_cns).value_counts()
            probabilities = value_counts / total_genome_bins
            sample_features['CN_Shannon_Entropy'] = entropy(probabilities, base=2)
        else:
            sample_features['Fraction_Genome_Altered'] = 0
            sample_features['Fraction_Genome_Gain'] = 0
            sample_features['Fraction_Genome_Lost'] = 0
            sample_features['CN_Shannon_Entropy'] = 0
            
        wide_data.append(sample_features) 
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    wide_df = pd.DataFrame(wide_data)
    wide_df.to_csv(OUTPUT_FILE, index=False)

def main() -> None:
    df_merged = pd.read_csv(INPUT_FILE)
    extract_wide_features(df_merged)


if __name__ == "__main__":
    main()
