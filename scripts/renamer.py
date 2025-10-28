import os
import pandas as pd
import re

meta_file = '/home/subpolare/gtrd/meta/meta_6_may.tsv'
vcf_dir = '/home/subpolare/adastra-v7/VCFs'

meta_df = pd.read_csv(meta_file, sep = '\t')
algn_to_gsm = dict(zip(meta_df['algn_id'], meta_df['geo_gsm']))
vcf_pattern = re.compile(r'(?P<tf>[^_]+)_(?P<cell>[^_]+)_(?P<algn>ALIGNS\d+)\.vcf\.gz')

for filename in os.listdir(vcf_dir):
    match = vcf_pattern.match(filename)
    if match:
        tf = match.group('tf')
        cell = match.group('cell')
        algn = match.group('algn')
        
        geo_gsm = algn_to_gsm.get(algn)
        if geo_gsm:
            new_filename = f'{tf}_{cell}_{algn}_{geo_gsm}.vcf.gz'
            old_path = os.path.join(vcf_dir, filename)
            new_path = os.path.join(vcf_dir, new_filename)
            
            os.rename(old_path, new_path)
        else:
            print(f'Пропущено (нет geo_gsm): {filename}')
