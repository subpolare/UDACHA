import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist


def compute_clusters_for_cutoff(args):
    mat, zero_cutoff, cluster_cutoff = args
    x = mat.copy()
    x[x < zero_cutoff] = 0.0
    x = np.nan_to_num(x, nan = 0.0, posinf = 0.0, neginf = 0.0)
    d = pdist(x, metric = 'correlation')
    d = np.nan_to_num(d, nan = 1.0, posinf = 1.0, neginf = 1.0)
    z = hierarchy.linkage(d, method = 'complete')
    cl = hierarchy.fcluster(z, cluster_cutoff, criterion = 'distance')
    return float(zero_cutoff), int(np.unique(cl).size)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', type = str, default = 'king_matrix.king')
    parser.add_argument('--threads', type = int, default = 1)
    parser.add_argument('--outpath', type = str, default = './clustering')
    parser.add_argument('--cutoff-min', type = float, default = 0.0)
    parser.add_argument('--cutoff-max', type = float, default = 0.8)
    parser.add_argument('--cutoff-step', type = float, default = 0.05)
    parser.add_argument('--cluster-cutoff', type = float, default = 0.1)
    args = parser.parse_args()

    matrix_path = Path(args.matrix)
    outpath = Path(args.outpath)
    outpath.mkdir(parents = True, exist_ok = True)

    mat = np.loadtxt(matrix_path, dtype = np.float32)
    mat = np.nan_to_num(mat, nan = 0.0, posinf = 0.0, neginf = 0.0)

    cutoffs = np.arange(args.cutoff_min, args.cutoff_max + args.cutoff_step / 2, args.cutoff_step, dtype = np.float32)

    tasks = [(mat, c, args.cluster_cutoff) for c in cutoffs]

    if args.threads <= 1:
        results = [compute_clusters_for_cutoff(t) for t in tasks]
    else:
        with ThreadPoolExecutor(max_workers = args.threads) as ex:
            results = list(ex.map(compute_clusters_for_cutoff, tasks))

    results = sorted(results, key = lambda x: x[0])
    xs = np.array([r[0] for r in results], dtype = float)
    ys = np.array([r[1] for r in results], dtype = int)

    tsv_path = outpath / 'cutoff_vs_clusters.tsv'
    with open(tsv_path, 'w') as f:
        f.write('zero_cutoff\tclusters\n')
        for x, y in results:
            f.write(f'{x}\t{y}\n')

    plt.style.use('ggplot')
    plt.figure(figsize = (18, 9))
    plt.plot(xs, ys)
    plt.xlabel('Zeroing cutoff for KING matrix (mat[mat < cutoff] = 0)')
    plt.ylabel(f'Number of clusters at dendrogram cutoff {args.cluster_cutoff}')
    plt.tight_layout()
    plt.savefig(outpath / 'cutoff_vs_clusters.png', dpi = 200)
    plt.show()


if __name__ == '__main__':
    main()
