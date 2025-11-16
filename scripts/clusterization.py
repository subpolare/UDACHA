import argparse
import numpy as np
import pandas as pd

from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster

def load_sample_ids(king_id_path: str) -> np.ndarray:
    ids = pd.read_csv(king_id_path, sep = r"\s+", header = None, engine = "python")
    if ids.shape[1] == 1:
        sid = ids.iloc[:, 0].astype(str).values
    else:
        sid = ids.iloc[:, 1].astype(str).values
    return sid

def main():
    ap = argparse.ArgumentParser(description = "Genotype-based clustering from KING matrix (complete-linkage, correlation metric)")
    ap.add_argument("--king",    required = True, help = "PLINK2 --make-king square output (.king)")
    ap.add_argument("--king-id", required = True, help = "Sample IDs for KING matrix (.king.id)")
    ap.add_argument("--output",  required = True, help = "Output TSV: sample_id\\tgeno_cluster_id")
    ap.add_argument("--dist0",   type = float, default = 0.4, help = "Threshold on inter-sample KING distances to be zeroed")
    ap.add_argument("--cut",     type = float, default = 0.1, help = "Dendrogram cut height in correlation distance")
    args = ap.parse_args()

    sample_id = load_sample_ids(args.king_id)
    n = len(sample_id)

    K = np.loadtxt(args.king, dtype = np.float32)
    if K.shape != (n, n):
        raise ValueError(f"KING matrix shape {K.shape} != ({n}, {n}) from {args.king_id}")

    np.fill_diagonal(K, 0.0)
    K[K < args.dist0] = 0.0

    corr_dist = pdist(K, metric = "correlation")
    Z = linkage(corr_dist, method = "complete")
    labels = fcluster(Z, t = args.cut, criterion = "distance")

    out = pd.DataFrame({"sample_id" : sample_id, "geno_cluster_id" : labels.astype(int)})
    out.to_csv(args.output, sep = "\t", index = False)

if __name__ == "__main__":
    main()
