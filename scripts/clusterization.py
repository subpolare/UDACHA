import numpy as np
import pandas as pd
import os, sys, argaprs
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster

def main():
    ap = argparse.ArgumentParser(description = 'Genotype-based clustering (complete-linkage, correlation)')
    ap.add_argument("--geno",   required = True, help = 'PLINK --export A file (geno.raw)')
    ap.add_argument("--output", required = True, help = 'Output TSV: sample_id\tgeno_cluster_id')
    args = ap.parse_args()

    X = pd.read_csv(args.geno, sep = r'\s+', engine = 'python')
    sid = X['IID']
    G = X.iloc[:,6:].replace(-9, np.nan).astype('float32')
    col_mean = G.mean(axis = 0); G = G.fillna(col_mean)
    std = G.std(axis = 0).replace(0, 1.0)
    G = (G - col_mean) / std

    D = squareform(pdist(G.values, metric = 'correlation'))
    D[D < 0.4] = 0.0
    Z = linkage(squareform(D, checks=False), method = 'complete')
    labels = fcluster(Z, t=0.1, criterion = 'distance')

    out = pd.DataFrame({'sample_id' : sample_id, 'geno_cluster_id' : labels})
    out.to_csv(args.output, sep = '\t', index = False)

if __name__ == "__main__":
    main()
