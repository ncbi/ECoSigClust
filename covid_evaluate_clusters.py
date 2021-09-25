#!/usr/bin/env python
# coding: utf-8


# Import packages
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from sklearn import metrics
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


# Input: type(n_clusters)     = range
#        type(corr_consensus) = DataFrame
#        type(cluster_ids)    = DataFrame
#        type(s_type)         = string, possibilities are silhouette, calinski, davies
# Output: list
# Goal: To evaluate the optimal number of clusters to use by giving each cluster number a score
def eval_cluster_scores(n_clusters, corr_consensus, cluster_ids, s_type='silhouette'):
    scores = []
    for cluster_num in n_clusters:
        if s_type == "silhouette":
            scores.append(
                metrics.silhouette_score(corr_consensus, cluster_ids[f"k={cluster_num}"].values, metric='euclidean'))
        elif s_type == "calinski":
            scores.append(metrics.calinski_harabasz_score(corr_consensus, cluster_ids[f"k={cluster_num}"].values))
        elif s_type == "davies":
            scores.append(metrics.davies_bouldin_score(corr_consensus, cluster_ids[f"k={cluster_num}"].values))
    return scores


# Input: type(n_clusters) = range, the list of cluster numbers that were evaluated with the given score type
#        type(scores)     = list of floats, the score for each cluster number in the n_clusters range
#        type(s_type)     = string, score_type used in evaluating clusters
#        type(filename)   = string, name under which to save graphs
# Output: None
# Goal: Create a lineplot with the scores for each cluster number
def plot_scores(n_clusters, scores, s_type):
    fig = plt.figure()
    ax = sns.lineplot(x=n_clusters, y=scores, markers="o")
    ax.set_xlabel('Cluster Label')
    ax.set_ylabel(f"{s_type.capitalize()} Scores")
    ax.set_title(f"Consensus of KMeans {s_type.capitalize()} Scores")
    plt.xticks(n_clusters)
    pdf.savefig(fig)


if __name__ == "__main__":
    ## Initialization
    # Grab file paths and parameters
    corr_consensus_file = snakemake.input[0]
    cluster_id_file = snakemake.input[1]
    cluster_eval_scores_file = snakemake.output[0]
    # score_type = snakemake.params.score_type
    # Set number of clusters to evaluate on
    num_clusters = range(5, 21)

    # Read in files
    corr_consensus_df = pd.read_csv(corr_consensus_file, sep="\t", index_col=0)
    cluster_ids_df = pd.read_csv(cluster_id_file, sep="\t", index_col=0)

    score_type = ["silhouette", "calinski", "davies"]

    pdf = PdfPages(cluster_eval_scores_file)
    for s_type in score_type:
        clustering_scores = eval_cluster_scores(num_clusters, corr_consensus_df, cluster_ids_df, s_type)
        plot_scores(num_clusters, clustering_scores, s_type)
    pdf.close()