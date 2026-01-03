import os, argparse
import numpy as np
import pandas as pd
import seaborn as sns
from curses import meta
import matplotlib.pyplot as plt
from scipy.cluster import hierarchy
from sklearn.metrics import silhouette_samples
plt.style.use('ggplot')


def visualize_clustering(mat, linkage, out_path):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 10))
    dendro = hierarchy.dendrogram(linkage, no_plot=False, ax=ax1)
    g = sns.heatmap(mat.iloc[dendro['leaves'],dendro['leaves']], cmap='Blues',
     square=True, xticklabels=False, yticklabels=False, cbar=False, ax=ax2)
    for l in ['left','right','top','bottom']:
        g.spines[l].set_visible(True)
        g.spines[l].set_color('k')
    
    plt.savefig(f"{out_path}.png")
    plt.close(fig)


def visualize_silhouette(sim_mat, labels, out_dir):
    labels = np.asarray(labels)
    uniq = np.unique(labels)

    if uniq.size < 2 or uniq.size >= labels.size:
        return

    dist = sim_mat.to_numpy(copy=False).astype(np.float32, copy=False)
    max_sim = float(np.nanmax(dist))
    dist[:] = max_sim - dist
    dist[dist < 0] = 0.0
    np.fill_diagonal(dist, 0.0)

    s = silhouette_samples(dist, labels, metric='precomputed')
    order = np.argsort(-s)
    s_sorted = s[order]

    fig = plt.figure(figsize=(18, 6))
    plt.plot(np.arange(s_sorted.shape[0]), s_sorted)
    plt.xlabel('Samples (sorted by silhouette score, desc)')
    plt.ylabel('Silhouette score')
    plt.tight_layout()

    out_file = os.path.join(out_dir, 'silhoette-score.png')
    plt.savefig(out_file, dpi=200)
    plt.close(fig)


def main(input_matrix, input_matrix_ids, meta_path, outpath):
    new_meta_path = os.path.join(outpath, "metadata.clustered.tsv") 
    indivs = np.loadtxt(input_matrix_ids, skiprows=0, dtype=str)
    rel_mat = np.loadtxt(input_matrix)
    mat = pd.DataFrame(rel_mat, index=indivs, columns=indivs)
    mat[mat < 0.4] = 0
    mat[np.isnan(mat)] = 0
    linkage = hierarchy.linkage(mat, method='complete', metric='correlation')
    cl = hierarchy.fcluster(linkage, 0.1, criterion='distance')
	visualize_silhouette(mat, cl, outpath)
    clusters = pd.DataFrame({'indiv_id': mat.index, 'genotype_cluster': cl}).sort_values(
        by='genotype_cluster')
    
    metadata = pd.read_table(meta_path, header=0, dtype={'indiv_id': str})
    # metadata = metadata.rename(columns={'indiv_id': 'ds_number'})
    metadata = metadata.merge(clusters, on='indiv_id').sort_values(by='genotype_cluster')
    metadata.rename(columns={'indiv_id': 'old_indiv_id'}, inplace=True)
    metadata.rename(columns={'genotype_cluster': 'indiv_id'}, inplace=True)
    metadata['indiv_id'] = 'INDIV_' + metadata['indiv_id'].astype(str).str.zfill(4)
    metadata.to_csv(new_meta_path, header=True, index=False, sep='\t')
    # visualizations_path = os.path.join(outpath, 'clustering')
    # visualize_clustering(mat, linkage, out_path=visualizations_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Count tags by allele")
    parser.add_argument("--matrix", type=str, help="Result of plink2 --king with suffix .king")

    parser.add_argument("--matrix-ids", type=str,
						help="Result of plink2 --king with suffix .king.id")

    parser.add_argument("--meta-file", type=str, 
						help="Path to meta file")

    parser.add_argument("--outpath", type=str, 
						help="Path to directory to save updated metafile and visualizations")
    

    args = parser.parse_args()

    main(args.matrix, args.matrix_ids, args.meta_file, args.outpath)
