#!/usr/bin/env python3

import sys
import glob
import argparse
import warnings
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')

def main():
    # Argument Parsing
    parser = argparse.ArgumentParser(description = 'Script to merge MixALiME output and multiple BED files into one final TSV table.')
    parser.add_argument('--mixalime', required = True, help = 'TSV file with MixALiME output.')
    parser.add_argument('--bed', required = True, help = 'Glob pattern to find non-archived BED files.')
    parser.add_argument('--output', required = True, help = 'Name of the final TSV table.')
    args = parser.parse_args()

    # Reading MixALiME TSV File
    try:
        df_mix = pd.read_csv(args.mixalime, sep = '\t')
    except Exception as e:
        sys.exit('Error reading MixALiME file: ' + str(e))

    # Creating final DataFrame with selected columns from MixALiME file
    try:
        df_final = pd.DataFrame({
            'chr': df_mix['#chr'],
            'start': df_mix['start'],
            'end': df_mix['end'],
            'ID': df_mix['id'],
            'ref': df_mix['ref'],
            'alt': df_mix['alt'],
            'mean_BAD': df_mix['mean_bad'],
            'n_aggregated': df_mix['n_reps'],
            'es_mean_ref': df_mix['ref_comb_es'],
            'es_mean_alt': df_mix['alt_comb_es'],
            'fdrp_bh_ref': df_mix['ref_fdr_comb_pval'],
            'fdrp_bh_alt': df_mix['alt_fdr_comb_pval']
        })
    except KeyError as err:
        sys.exit('Required column not found in MixALiME file: ' + str(err))

    # Adding empty columns for placeholders
    empty_cols = ['repeat_type', 'motif_log_pref', 'motif_log_palt', 'motif_fc', 'motif_pos', 'motif_orient', 'motif_conc', 'motif_index']
    for col in empty_cols:
        df_final[col] = ''

    # Creating temporary key column 'id_ref_alt' in final DataFrame
    df_final['id_ref_alt'] = df_final['ID'].astype(str) + '_' + df_final['ref'].astype(str) + '_' + df_final['alt'].astype(str)

    # Reading and combining BED files
    bed_files = glob.glob(args.bed)
    if not bed_files:
        sys.exit('No BED files found with pattern: ' + args.bed)

    bed_dfs = []
    for bf in bed_files:
        try:
            df_bed = pd.read_csv(bf, delim_whitespace = True)
            bed_dfs.append(df_bed)
        except Exception as e:
            sys.exit('Error reading BED file ' + bf + ': ' + str(e))
    df_bed_all = pd.concat(bed_dfs, ignore_index = True)

    # Sorting BED data and creating temporary key column 'id_ref_alt'
    if '#chr' not in df_bed_all.columns:
        sys.exit('BED files missing column "#chr".')
    df_bed_all.sort_values(by = ['#chr', 'start'], inplace = True)
    for col in ['id', 'ref', 'alt']:
        if col not in df_bed_all.columns:
            sys.exit('BED files missing column ' + col + '.')
    df_bed_all['id_ref_alt'] = df_bed_all['id'].astype(str) + '_' + df_bed_all['ref'].astype(str) + '_' + df_bed_all['alt'].astype(str)

    # Aggregating BED data: sum of total_cover and mean of SNP_per_segment
    for col in ['total_cover', 'SNP_per_segment']:
        if col not in df_bed_all.columns:
            sys.exit('BED files missing column ' + col + '.')
    bed_agg = df_bed_all.groupby('id_ref_alt').agg({'total_cover': 'sum', 'SNP_per_segment': 'mean'}).reset_index()

    # Merging aggregated BED data with final DataFrame
    df_final = df_final.merge(bed_agg, on = 'id_ref_alt', how = 'left')

    # Checking for missing aggregated data
    missing = df_final['total_cover'].isnull() | df_final['SNP_per_segment'].isnull()
    if missing.any():
        missing_ids = df_final.loc[missing, 'id_ref_alt'].unique()
        sys.exit('No matching BED entries found for id_ref_alt: ' + ', '.join(missing_ids))

    # Renaming aggregated column and cleaning up temporary columns
    df_final.rename(columns = {'SNP_per_segment': 'mean_SNP_per_segment'}, inplace = True)
    df_final.drop(columns = ['id_ref_alt'], inplace = True)

    # Reordering final DataFrame columns
    final_order = ['chr', 'start', 'end', 'ID', 'ref', 'alt', 'repeat_type', 'mean_BAD', 'mean_SNP_per_segment', 'n_aggregated', 'total_cover', 'es_mean_ref', 'es_mean_alt', 'fdrp_bh_ref', 'fdrp_bh_alt', 'motif_log_pref', 'motif_log_palt', 'motif_fc', 'motif_pos', 'motif_orient', 'motif_conc', 'motif_index']
    df_final = df_final[final_order]

    # Rounding columns to desired precision
    df_final['mean_SNP_per_segment'] = df_final['mean_SNP_per_segment'].round(1)
    df_final['es_mean_ref'] = df_final['es_mean_ref'].round(3)
    df_final['es_mean_alt'] = df_final['es_mean_alt'].round(3)
    df_final['mean_BAD'] = df_final['mean_BAD'].round(2)

    # Saving final TSV file
    try:
        df_final.to_csv(args.output, sep = '\t', index = False)
    except Exception as e:
        sys.exit('Error saving final file: ' + str(e))

if __name__=='__main__':
    main()
