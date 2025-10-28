import os
import subprocess
from collections import defaultdict

VCF_DIR = '/home/subpolare/adastra-v7/VCFs'
OUTPUT_DIR = '/home/subpolare/GEO_GSE'

os.makedirs(OUTPUT_DIR, exist_ok = True)

geo_groups = defaultdict(list)
nan_files = []

for filename in os.listdir(VCF_DIR):
    if filename.endswith('.vcf.gz'):
        parts = filename.split('_')
        geo_id = parts[-1]
        
        if geo_id == 'nan.vcf.gz':
            unique_id = parts[-2]
            new_name = f'NaN_{unique_id}.vcf.gz'
            nan_files.append((os.path.join(VCF_DIR, filename), os.path.join(OUTPUT_DIR, new_name)))
        else:
            geo_groups[geo_id].append(os.path.join(VCF_DIR, filename))

for geo_id, file_list in geo_groups.items():
    if len(file_list) <= 1:
        continue

    output_file = os.path.join(OUTPUT_DIR, geo_id)
    temp_file = f'/tmp/{geo_id}.txt'
    
    with open(temp_file, 'w') as f:
        f.write('\n'.join(file_list) + '\n')

    subprocess.run(['bcftools', 'merge', '-l', temp_file, '-Oz', '-o', output_file])
    os.remove(temp_file)

for input_file, output_file in nan_files:
    subprocess.run(['bcftools', 'view', '-Oz', '-o', output_file, input_file])
