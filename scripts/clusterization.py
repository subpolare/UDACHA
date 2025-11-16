import argparse
import numpy as np
import pandas as pd

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster

def main():
    ap = argparse.ArgumentParser(description = 'Genotype-based clustering from KING matrix (complete-linkage, correlation metric)')
    ap.add_argument('--king',    required = True, help = 'PLINK2 --make-king square output (.king)')
    ap.add_argument('--king-id', required = True, help = 'Sample IDs for KING matrix (.king.id)')
    ap.add_argument('--output',  required = True, help = 'Output TSV: sample_id\tgeno_cluster_id')
    ap.add_argument('--dist0',   type = float, default = 0.4, help = 'Threshold on inter-sample KING distances to be zeroed')
    ap.add_argument('--cut',     type = float, default = 0.1, help = 'Dendrogram cut height in correlation distance')
    args = ap.parse_args()

    K = np.loadtxt(args.king, dtype = np.float32)
    n_k = K.shape[0]
    print(f'✔️ [1/5] Load KING matrix with shape {K.shape}')

    ids = pd.read_csv(args.king_id, sep = r'\s+', engine = 'python')
    if ids.shape[0] == n_k + 1:
        ids = ids.iloc[1:, :]
    if ids.shape[0] != n_k:
        raise ValueError(f'KING ID count {ids.shape[0]} != {n_k} rows in {args.king}')
    if ids.shape[1] == 1:
        sample_id = ids.iloc[:, 0].astype(str).values
    else:
        sample_id = ids.iloc[:, 1].astype(str).values
    print(f'✔️ [2/5] Load {len(sample_id)} sample IDs')

    np.fill_diagonal(K, 0.0)
    K[K < args.dist0] = 0.0
    print(f'✔️ [3/5] Zero KING distances < {args.dist0}')

    corr_dist = pdist(K, metric = 'correlation')
    Z = linkage(corr_dist, method = 'complete')
    labels = fcluster(Z, t = args.cut, criterion = 'distance')
    print(f'✔️ [4/5] Perform hierarchical clustering with cut {args.cut}')

    out = pd.DataFrame({'sample_id' : sample_id, 'geno_cluster_id' : labels.astype(int)})
    out.to_csv(args.output, sep = '\t', index = False)
    print(f'✔️ [5/5] Save clusters to {args.output}')

if __name__ == '__main__':
    main()
