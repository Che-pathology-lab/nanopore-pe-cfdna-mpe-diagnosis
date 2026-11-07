"""Purpose: Chip cna feature for the manuscript analysis workflow.

Edit the paths and parameters below, then run this script directly.
"""

import pandas as pd
import numpy as np
from scipy.stats import entropy


INPUT_FILE = "results/08_chip_validation/CHIP_CNA/CHIP_shift_ratio.csv"
OUTPUT_FILE = "results/08_chip_validation/CHIP_CNA/CHIP_CNA_features_wide_table.csv"

def map_to_cn_state(shift_val):
    if pd.isna(shift_val): return np.nan
    if shift_val <= -0.96: return 0
    elif shift_val < -0.24: return 1
    elif shift_val <= 0.12: return 2
    elif shift_val <= 0.72: return 3
    else: return 4


ICHORCNA_REGIONS = {
    'chr1': (585989, 248946422), 'chr2': (10001, 242183529), 'chr3': (10001, 198235559),
    'chr4': (10001, 190000000), 'chr5': (10001,  181478259), 'chr6': (60001, 170745979),
    'chr7': (240243, 159335973), 'chr8': (60001, 145000000), 'chr9': (10001, 138334717),
    'chr10': (10001, 133690466), 'chr11': (60001, 130000000), 'chr12': (10001, 133265309),
    'chr13': (20000001, 114354328), 'chr14': (20000001, 106883718), 'chr15': (25000001, 101981189),
    'chr16': (10001, 90000000), 'chr17': (491112, 83247441), 'chr18': (10001, 80000000),
    'chr19': (60001, 58607616), 'chr20': (66336, 64334167), 'chr21': (15000001, 46699983), 
    'chr22': (16314011, 50808468)
}

CENTROMERES = {
    'chr1': (120000000, 145000001), 'chr2': (89330679, 95000001), 'chr3': (75000000, 95000001),
    'chr4': (49336924, 51793952), 'chr5': (45000000, 50109808), 'chr6': (58453888, 70000001),
    'chr7': (58119653, 65000001), 'chr8': (43983744, 45927266),'chr9': (35000000, 70000001),
    'chr10': (38529907, 42066266), 'chr11': (50000000, 60000001), 'chr12': (34719407, 40000001),
    'chr13': (0, 20000001), 'chr14': (0, 20000001), 'chr15': (0, 25000001), 'chr16': (33214595, 46380683),
    'chr17': (21795850, 26935981), 'chr18': (10000000, 20911207), 'chr19': (24448980, 27240875),
    'chr20': (26348365, 31161626), 'chr21': (0, 15000001), 'chr22': (0, 16314011)
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
    exclude_cols = ['chr', 'start', 'end', 'arm', 'probes']
    sample_cols = [col for col in df.columns if col not in exclude_cols]
    
    wide_data = []
    for col in sample_cols:
        sample_name = col
        sample_features = {'Sample': sample_name}
        
        df[col] = pd.to_numeric(df[col], errors='coerce')
        df['mapped_cn'] = df[col].apply(map_to_cn_state)
        s_df = df[['chr', 'start', 'end', 'arm', 'mapped_cn']].dropna(subset=['mapped_cn']).sort_values(by=['chr', 'start'])
        total_genome_bins = len(s_df)
        altered_bins = 0 
        lost_bins = 0    
        gain_bins = 0

        aberrant_sizes_mb = []
        all_cns = []
        total_changepoints_count = 0
        bp_per_arm = []
        bp_10mb = []
        osc_lengths = []
        
        for arm in s_df['arm'].unique():
            arm_df = s_df[s_df['arm'] == arm]
            if arm_df.empty: continue
            
            cns = arm_df['mapped_cn'].values
            starts = arm_df['start'].values
            ends = arm_df['end'].values
            all_cns.extend(cns)
            
            altered_bins += np.sum(cns != 2)
            lost_bins += np.sum(cns < 2)
            gain_bins += np.sum(cns > 2)
            
            segments = []       
            arm_bp_coords = []  
            current_cn = cns[0]
            seg_start = starts[0]
            
            for i in range(1, len(cns)):
                cn = cns[i]
                if cn != current_cn:
                    seg_end = ends[i-1]
                    segments.append((current_cn, seg_start, seg_end))
                    arm_bp_coords.append(seg_end) 
                    current_cn = cn
                    seg_start = starts[i]
                    
            segments.append((current_cn, seg_start, ends[-1]))

            for cn, s_start, s_end in segments:
                if cn != 2:
                    true_length_mb = (s_end - s_start) / 1000000.0
                    aberrant_sizes_mb.append(true_length_mb)
            
            num_changepoints = len(segments) - 1
            total_changepoints_count += num_changepoints
            bp_per_arm.append(num_changepoints)
            
            chrom = arm_df['chr'].iloc[0]
            if chrom in ICHORCNA_REGIONS and arm_bp_coords:
                chr_start, chr_len = ICHORCNA_REGIONS[chrom]
                bins_10mb = np.arange(chr_start, chr_len + 10000000, 10000000)
                counts, _ = np.histogram(arm_bp_coords, bins=bins_10mb)
                bp_10mb.extend(counts)
            
            seg_cns = [seg[0] for seg in segments]
            osc_lengths.extend(get_strict_oscillating_lengths(seg_cns))

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
            value_counts = pd.Series(all_cns).value_counts()
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
