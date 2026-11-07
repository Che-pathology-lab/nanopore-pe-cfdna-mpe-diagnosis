"""Purpose: Nanopore for chip cna feature for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import numpy as np
import warnings
from scipy.stats import entropy
warnings.filterwarnings('ignore')


INPUT_FILE = "results/03_nanopore_cna/discovery_Corrected_Copy_Number.csv"
OUTPUT_FILE = "results/08_chip_validation/ONT_for_CHIP_CNA/discovery_ONT_for_CHIP_CNA_features_wide_table.csv"
    
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

def get_strict_oscillating_lengths(seg_states):
    if len(seg_states) < 3: return []
    lengths = []
    i = 0
    while i < len(seg_states) - 2:
        if seg_states[i] == seg_states[i+2] and seg_states[i] != seg_states[i+1]:
            state_a, state_b = seg_states[i], seg_states[i+1]
            jump_count = 1
            j = i + 3
            while j < len(seg_states):
                expected = state_a if (j - i) % 2 == 0 else state_b
                if seg_states[j] == expected:
                    jump_count += 1
                    j += 1
                else:
                    break
            lengths.append(jump_count)
            i = j - 1 
        else:
            i += 1
    return lengths

def map_cn_to_state(cn):
    if pd.isna(cn):
        return np.nan
    cn = round(float(cn))
    if cn <= 0: return 0     # Homozygous deletion
    elif cn == 1: return 1   # Hemizygous deletion
    elif cn == 2: return 2   # Neutral
    elif cn in [3,4]: return 3 # Gain
    elif cn >= 5: return 4   # High-level amplification
    return np.nan

def extract_wide_features(df):
    
    df['arm'] = df.apply(lambda row: get_arm(row['chr'], row['start'], row['end']), axis=1)
    df = df.dropna(subset=['arm'])
    
    sample_cols = [col for col in df.columns if col not in ['chr', 'start', 'end', 'arm']]
    wide_data = []
    
    for col in sample_cols:
        sample_name = col.replace('.Corrected_Copy_Number', '')
        sample_features = {'Sample': sample_name}
        df[col] = pd.to_numeric(df[col], errors='coerce')
        s_df = df[['chr', 'start', 'end', 'arm', col]].dropna(subset=[col]).copy()
        s_df[col] = s_df[col].apply(map_cn_to_state)
        s_df = s_df.dropna(subset=[col]).sort_values(by=['chr', 'start'])
        
        total_genome_bins = len(s_df)
        altered_bins = 0 # State != 2
        lost_bins = 0    # State < 2
        gain_bins = 0    # State > 2

        aberrant_sizes_mb = []
        all_states = []
        total_changepoints_count = 0
        bp_per_arm = []
        bp_10mb = []
        osc_lengths = []
        
        for arm in s_df['arm'].unique():
            arm_df = s_df[s_df['arm'] == arm]
            if arm_df.empty: continue
            
            states = arm_df[col].values
            ends = arm_df['end'].values
            
            all_states.extend(states)
            altered_bins += np.sum(states != 2)
            lost_bins += np.sum(states < 2)
            gain_bins += np.sum(states > 2)
            
            segments = []       
            arm_bp_coords = []  
            current_state = states[0]
            current_bins = 1
            for i in range(1, len(states)):
                state = states[i]
                if state == current_state:
                    current_bins += 1
                else:
                    segments.append((current_state, current_bins))
                    arm_bp_coords.append(ends[i-1]) 
                    current_state = state
                    current_bins = 1
            segments.append((current_state, current_bins))
            
            for state, bins_count in segments:
                if state != 2:
                    aberrant_sizes_mb.append(bins_count * 5)
            
            num_changepoints = len(segments) - 1
            total_changepoints_count += num_changepoints
            bp_per_arm.append(num_changepoints)
            
            chrom = arm_df['chr'].iloc[0]
            if chrom in ICHORCNA_REGIONS and arm_bp_coords:
                chr_start, chr_len = ICHORCNA_REGIONS[chrom]
                bins_10mb = np.arange(chr_start, chr_len + 10000000, 10000000)
                counts, _ = np.histogram(arm_bp_coords, bins=bins_10mb)
                bp_10mb.extend(counts)
            

            seg_states = [state for state, bins_count in segments]
            osc_lengths.extend(get_strict_oscillating_lengths(seg_states))

        sample_features['Max_Aberrant_Seg_Size'] = np.max(aberrant_sizes_mb) if aberrant_sizes_mb else 0
        sample_features['Mean_Aberrant_Seg_Size'] = np.mean(aberrant_sizes_mb) if aberrant_sizes_mb else 0

        sample_features['Total_Changepoints'] = total_changepoints_count
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
            
            value_counts = pd.Series(all_states).value_counts()
            probabilities = value_counts / total_genome_bins
            sample_features['CN_Shannon_Entropy'] = entropy(probabilities, base=2)
        else:
            sample_features['Fraction_Genome_Altered'] = 0
            sample_features['Fraction_Genome_Gain'] = 0
            sample_features['Fraction_Genome_Lost'] = 0
            sample_features['CN_Shannon_Entropy'] = 0
            
        wide_data.append(sample_features)
        
    wide_df = pd.DataFrame(wide_data)
    wide_df.to_csv(OUTPUT_FILE, index=False)

def main() -> None:
    df_merged = pd.read_csv(INPUT_FILE)
    extract_wide_features(df_merged)


if __name__ == "__main__":
    main()
