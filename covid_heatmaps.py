#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import scipy.stats as ss
import statsmodels.stats.multitest as ssm

import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


# input files
corr_count_file = snakemake.input[0]
cluster_id_file = snakemake.input[1]

home_dir = "/home/kimy3/Projects/Covid/"
heatmap_file = snakemake.output[0]

# read correlation file
corr_count_df = pd.read_csv(corr_count_file, sep="\t", index_col=0)
gene_names = corr_count_df.iloc[:,0]
gene_name_dic = dict(zip(corr_count_df.index, gene_names)) # ensembl -> gene mapping
corr_count_df  = corr_count_df.iloc[:,1:]

cluster_id_df = pd.read_csv(cluster_id_file, sep="\t", index_col=0)


matrix_to_cluster = corr_count_df

start = snakemake.params.start
end = snakemake.params.end
nrange = range(start, end)

pdf = PdfPages(heatmap_file)
for n_clusters in nrange:
    groups  = cluster_id_df.groupby("k="+str(n_clusters))
    rows, lens = {}, {}
    for i, group in groups:
        rows[i] = matrix_to_cluster.loc[group.index].mean()
        lens[i] = group.shape[0]
    sns.set(font_scale=1)
    df = pd.DataFrame([rows[i] for i  in range(1, n_clusters+1)])
    df.index=["CL"+str(i+1)+"("+str(lens[i+1])+")" for i in range(n_clusters)]
    plt.clf()
    #print(df)
    plt.tight_layout()
    cm = sns.clustermap(df, cmap=plt.cm.bwr, center=0, annot=True, square=False, annot_kws={"size": 10}, vmin=-0.4, vmax=0.5)
    cm.fig.suptitle("k="+str(n_clusters))
    pdf.savefig(cm.fig)
pdf.close()


# print gene files
for n_clusters in nrange:
    gene_file = snakemake.output[1+n_clusters-nrange[0]]
    print(n_clusters, gene_file)
    cluster_ids = cluster_id_df["k="+str(n_clusters)]
    gene_df = pd.DataFrame(index=matrix_to_cluster.index)
    gene_df["gene"] = [gene_name_dic[id] for id in matrix_to_cluster.index]
    gene_df["label"] = cluster_ids[cluster_ids > 0]
    gene_groups = {}
    for (l, group_df) in gene_df.groupby("label"):
        gene_groups[l] = group_df.gene.values

    gene_rows = []
    for i in range(n_clusters):
        print(i)
        cluster_id = i+1
        my_genes = list(gene_groups[i+1])
        my_genes.sort()
        gene_rows.append(("CL"+str(i+1), ";".join(my_genes)))
    pd.DataFrame(gene_rows, columns=["CL", "gene"]).set_index("CL").to_csv(gene_file)
