from tqdm import tqdm
import pandas as pd
import warnings
import os

warnings.simplefilter(action = 'ignore', category = Warning)
folders = ['/home/subpolare/adastra-v7/SNPScan/pwm_results_0',
           '/home/subpolare/adastra-v7/SNPScan/pwm_results_1',
           '/home/subpolare/adastra-v7/SNPScan/pwm_results_2',
           '/home/subpolare/adastra-v7/SNPScan/pwm_results_3']

file_paths = dict()

for folder in folders:
    files = os.listdir(folder)
    for file in files:
        file_path = os.path.join(folder, file)
        if file in file_paths:
            file_paths[file].append(file_path)
        else:
            file_paths[file] = [file_path]

for file, paths in tqdm(file_paths.items(), desc = 'Обработка файлов', colour = 'green'):
    if len(paths) == 1:
        data = pd.read_csv(paths[0], sep = '\t')
        data['index'] = paths[0].split('/')[-2][-1]
        data.to_csv('/home/subpolare/adastra-v7/SNPScan/merged_results/' + paths[0].split('/')[-1], sep = '\t', index = False)

    else:
        df = pd.DataFrame()
        file_name = paths[0].split('/')[-1]
        for i, path in enumerate(paths):
            data = pd.read_csv(path, sep = '\t')
            data['index'] = paths[i].split('/')[-2][-1]
            df = pd.concat([df, data])
        df['Abs fold change'] = df['Fold change'].abs()
        df['min_P_value'] = df[['P-value 1', 'P-value 2']].min(axis = 1)
        df['round_P_value'] = df['min_P_value'].apply(lambda x: 0.001 if x < 0.001 else x)
        df = df.sort_values(by = ['SNP name', 'round_P_value', 'Abs fold change'], ascending = [True, True, False])
        df['UniqID'] = df['SNP name'] + df['allele 1/allele 2'].str.split('/').str[0] + df['allele 1/allele 2'].str.split('/').str[1]
        df.drop_duplicates(subset = ['UniqID'], keep = 'first', inplace = True)
        df.drop(columns = ['UniqID', 'Abs fold change', 'min_P_value', 'round_P_value'], inplace = True)
        df.to_csv('/home/subpolare/adastra-v7/SNPScan/merged_results/' + file_name, sep = '\t', index = False)
