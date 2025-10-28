from tqdm import tqdm
import pandas as pd
import numpy as np
import warnings
import shutil
import sys
import os
from concurrent.futures import ProcessPoolExecutor
warnings.simplefilter(action = 'ignore', category = Warning)

def process_file(file):
    tf = pd.read_csv(f'/home/subpolare/adastra-v7/new-version/TF/{file.split(".")[0]}_HUMAN.tsv', sep = '\t')
    tf['UniqID_1'] = tf['ID'] + tf['ref'] + tf['alt']

    my = pd.read_csv(f'/home/subpolare/adastra-v7/SNPScan/merged_results/{file}', sep = '\t')
    my['UniqID_2'] = my['SNP name'] + my['allele 1/allele 2'].str.split('/').str[0] + my['allele 1/allele 2'].str.split('/').str[1]

    merged = pd.merge(tf, my, left_on = 'UniqID_1', right_on = 'UniqID_2', how = 'inner')
    merged['Fold change'] = np.log2(merged['Fold change'])

    condition_1 = (merged[['fdrp_bh_ref', 'fdrp_bh_alt']].min(axis = 1) < 0.1)
    condition_2 = (merged[['P-value 1', 'P-value 2']].min(axis = 1) < 0.001)
    condition_3 = (merged['fdrp_bh_ref'] - merged['fdrp_bh_alt']) * (merged['P-value 1'] - merged['P-value 2']) < 0
    condition_4 = (merged['fdrp_bh_ref'] - merged['fdrp_bh_alt']) * (merged['P-value 1'] - merged['P-value 2']) > 0

    merged['motif_conc'] = np.select(
        [
            condition_1 & condition_2 & condition_3 & (abs(merged['Fold change']) < 2),
            condition_1 & condition_2 & condition_3 & (abs(merged['Fold change']) >= 2),
            condition_1 & condition_2 & condition_4 & (abs(merged['Fold change']) < 2),
            condition_1 & condition_2 & condition_4 & (abs(merged['Fold change']) >= 2),
            condition_1
        ],
        ['Weak Discordant', 'Discordant', 'Weak Concordant', 'Concordant', 'No Hit'],
        default = ''
    )

    mask = merged['P-value 1'] < merged['P-value 2']
    merged['motif_pos'] = mask * merged['position 1'] + ~mask * merged['position 2']
    merged['motif_orient'] = mask * merged['orientation 1'] + ~mask * merged['orientation 2']
    revcomp = merged['motif_orient'] == 'revcomp'
    merged['motif_orient'] = merged['motif_orient'].replace({'direct': '+', 'revcomp': '-'})
    merged['motif_pos'] = revcomp * (merged['word 1'].str.len() + merged['motif_pos']) + ~revcomp * (-merged['motif_pos'] + 1)
    merged['motif_pos'] = merged['motif_pos'] - 1

    merged['motif_log_pref'] = -np.log10(merged['P-value 1'])
    merged['motif_log_palt'] = -np.log10(merged['P-value 2'])
    merged['motif_fc'] = merged['Fold change']
    merged['motif_index'] = merged['index'].astype(int)

    merged.drop(columns=['SNP name', 'motif', 'position 1', 'orientation 1', 'word 1', 
                         'position 2', 'orientation 2', 'word 2', 'allele 1/allele 2', 
                         'P-value 1', 'P-value 2', 'Fold change', 'index', 'UniqID_1', 
                         'UniqID_2'], inplace = True)
    merged = merged.sort_values(by = ['chr', 'start'])
    merged.to_csv(f'/home/subpolare/adastra-v7/new-version/TF/{file.split(".")[0]}_HUMAN.tsv', sep = '\t', index = False)

files = [file for file in os.listdir('/home/subpolare/adastra-v7/SNPScan/merged_results/') if not file.startswith('.')]
with ProcessPoolExecutor(max_workers = 50) as executor:
    executor.map(process_file, files)
