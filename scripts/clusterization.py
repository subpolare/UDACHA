import argparse
import numpy as np
import pandas as pd

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster

def main():
    ap = argparse.ArgumentParser(description = 'Genotype-based clustering from PLINK --export A (complete-linkage, correlation metric)')
    ap.add_argument('--geno',    required = True, help = 'PLINK --export A .raw file')
    ap.add_argument('--output',  required = True, help = 'Output TSV: sample_id\tgeno_cluster_id')
    ap.add_argument('--cut',     type = float, default = 0.1, help = 'Dendrogram cut height in correlation distance')
    args = ap.parse_args()

    X = pd.read_csv(args.geno, sep = r'\s+', engine = 'python')
    print(f'✅ [1/6] Load PLINK raw with shape {X.shape}')

    if 'IID' not in X.columns:
        raise ValueError('IID column not found in PLINK raw file')
    sample_id = X['IID'].astype(str).values
    G = X.iloc[:, 6:].replace(-9, np.nan).astype('float32')
    print(f'✅ [2/6] Extract genotype matrix with shape {G.shape}')

    nan_frac = G.isna().mean(axis = 1).values
    mask = nan_frac < 1.0
    dropped = (~mask).sum()
    if dropped > 0:
        print(f'⚠️  Drop {dropped} samples with all NaN genotypes')
    G = G.loc[mask].copy()
    sample_id_clust = sample_id[mask]
    print(f'✅ [3/6] Keep {G.shape[0]} samples for clustering')

    col_mean = G.mean(axis = 0)
    G = G.fillna(col_mean)
    std = G.std(axis = 0).replace(0, 1.0)
    G = (G - col_mean) / std
    print('✅ [4/6] Standardize genotypes')

    D = squareform(pdist(G.values, metric = 'correlation'))
    Z = linkage(squareform(D, checks = False), method = 'complete')
    labels_clust = fcluster(Z, t = args.cut, criterion = 'distance')
    print(f'✅ [5/6] Perform hierarchical clustering with cut {args.cut}')

    labels = np.zeros(sample_id.shape[0], dtype = int)
    labels[mask] = labels_clust
    if dropped > 0:
        max_lab = labels_clust.max()
        idx = np.where(~mask)[0]
        labels[idx] = np.arange(max_lab + 1, max_lab + 1 + dropped)
        print(f'⚠️  Assign unique cluster IDs to {dropped} samples with all-missing genotypes')

    out = pd.DataFrame({'sample_id' : sample_id, 'geno_cluster_id' : labels.astype(int)})
    out.to_csv(args.output, sep = '\t', index = False)
    print(f'✅ [6/6] Save clusters to {args.output}')

if __name__ == '__main__':
    main()
