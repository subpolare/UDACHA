#!/usr/bin/env python3

import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Transfers raw p-values into ADASTRA table from MixALiME output table.")
    parser.add_argument('--adastra', required=True, help='Path to ADASTRA file')
    parser.add_argument('--mixalime', required=True, help='Path to MixALiME file')
    args = parser.parse_args()

    adastra = pd.read_csv(args.adastra, sep='\t')
    mixalime = pd.read_csv(args.mixalime, sep='\t')
    
    adastra['key'] = adastra['ID'].astype(str) + "_" + adastra['alt'].astype(str)
    mixalime['key'] = mixalime['id'].astype(str) + "_" + mixalime['alt'].astype(str)
    
    merged = pd.merge(adastra, mixalime[['key', 'ref_comb_pval', 'alt_comb_pval']], on='key', how='left')
    
    merged.rename(columns={'ref_comb_pval': 'pval_mean_ref',
                           'alt_comb_pval': 'pval_mean_alt'}, inplace=True)
    
    ordered_columns = [
        'chr', 'start', 'end', 'ID', 'ref', 'alt', 'repeat_type', 
        'mean_BAD', 'mean_SNP_per_segment', 'n_aggregated', 'total_cover', 
        'es_mean_ref', 'es_mean_alt', 'pval_mean_ref', 'pval_mean_alt', 
        'fdrp_bh_ref', 'fdrp_bh_alt', 'motif_log_pref', 'motif_log_palt', 
        'motif_fc', 'motif_pos', 'motif_orient', 'motif_conc', 'motif_index'
    ]
    
    result = merged[ordered_columns]
    
    result.to_csv(args.adastra, sep='\t', index=False)

if __name__ == '__main__':
    main()
