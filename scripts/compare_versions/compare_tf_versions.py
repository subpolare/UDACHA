import argparse, sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def parse_size(size_str):
    try:
        width, height = size_str.split(':')
        return (float(width), float(height))
    except Exception as e:
        print(f'Error parsing size: {e}')
        sys.exit(1)


def main(old_path, new_path, style, size):
    plt.style.use(style)
    
    df_old = pd.read_csv(old_path, sep='\t')
    df_new = pd.read_csv(new_path, sep='\t')
    
    df_old.columns = [col.strip() for col in df_old.columns]
    df_new.columns = [col.strip() for col in df_new.columns]
    
    ids_old = set(df_old['ID'])
    ids_new = set(df_new['ID'])
    common_ids = ids_old.intersection(ids_new)
    
    pct_old = len(common_ids) / len(ids_old) * 100
    pct_new = len(common_ids) / len(ids_new) * 100
    print(f'Common SNPs: {len(common_ids)}')
    print(f'Percentage of common SNPs in first table (--old): {pct_old:.2f}%')
    print(f'Percentage of common SNPs in second table (--new): {pct_new:.2f}%')
    
    if len(common_ids) == 0:
        print('No common SNPs for comparison. Exiting.')
        sys.exit(0)
    
    df_old_common = df_old[df_old['ID'].isin(common_ids)].copy()
    df_new_common = df_new[df_new['ID'].isin(common_ids)].copy()
    df_old_common.set_index('ID', inplace=True)
    df_new_common.set_index('ID', inplace=True)
    
    for df in [df_old_common, df_new_common]:
        df['motif_fc'] = pd.to_numeric(df['motif_fc'], errors='coerce')
        df['motif_conc'] = pd.to_numeric(df['motif_conc'], errors='coerce')
    
    df_old_common = df_old_common[['motif_fc', 'motif_conc']].rename(
        columns={'motif_fc': 'motif_fc_old', 'motif_conc': 'motif_conc_old'})
    df_new_common = df_new_common[['motif_fc', 'motif_conc']].rename(
        columns={'motif_fc': 'motif_fc_new', 'motif_conc': 'motif_conc_new'})
    
    df_merged = df_old_common.join(df_new_common, how='inner')
    df_merged = df_merged[(df_merged['motif_conc_old'].notnull()) | (df_merged['motif_conc_new'].notnull())].copy()
    if df_merged.empty:
        print('No SNP pairs meeting the condition (at least one motif_conc value present).')
        sys.exit(0)
    
    df_merged['fc_diff'] = df_merged['motif_fc_old'] - df_merged['motif_fc_new']
    df_merged['abs_fc_diff'] = df_merged['fc_diff'].abs()
    
    def determine_sign(row):
        if pd.notnull(row['motif_conc_old']) and pd.notnull(row['motif_conc_new']):
            return 1 if row['motif_conc_old'] == row['motif_conc_new'] else -1
        else:
            return -1
            
    df_merged['sign'] = df_merged.apply(determine_sign, axis=1)
    df_merged['y_value'] = df_merged['sign'] * df_merged['abs_fc_diff']
    
    plot_data = pd.DataFrame({'Comparison': 'SNPs', 'Similarity': df_merged['y_value']})
    fig, ax = plt.subplots(figsize=size)
    sns.violinplot(x='Comparison', y='Similarity', data=plot_data, inner=None, color='0.8', ax=ax)
    sns.stripplot(x='Comparison', y='Similarity', data=plot_data, color='black', size=4, jitter=True, ax=ax)
    
    ax.set_title('SNP Comparison (motif_fc difference with sign)', fontsize=14)
    ax.set_ylabel('SNP similarity (signed motif_fc difference)', fontsize=12)
    ax.set_xlabel('')
    
    sns.despine(left=True, bottom=True)
    plt.tight_layout()
    
    plt.savefig('comparison.png', transparent=True, dpi=300)
    print('Plot saved as comparison.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script for comparing TSV SNP tables and plotting violin plot.')
    parser.add_argument('--old', required=True, help='Path to the first (old) table')
    parser.add_argument('--new', required=True, help='Path to the second (new) table')
    parser.add_argument('--style', default='ggplot', help='Plot style (default: ggplot)')
    parser.add_argument('--size', default='16:9', help='Plot size in WIDTH:HEIGHT format (default: 16:9)')
    args = parser.parse_args()
    
    fig_size = parse_size(args.size)
    main(args.old, args.new, args.style, fig_size)
