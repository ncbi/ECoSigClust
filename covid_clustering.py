#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy 
import scipy.stats as ss
import scipy.cluster.hierarchy as sch

from sklearn.cluster import KMeans, MiniBatchKMeans
import sklearn.cluster as cluster
from sklearn import metrics
import statsmodels.stats.multitest as ssm


# input files
corr_count_file = snakemake.input[0]
pv_count_file = snakemake.input[1]

# output files
corr_consensus_file = snakemake.output[0]
cluster_id_file = snakemake.output[1]
selected_file = snakemake.output[2]

# parameters
pv_th = 0.05 

# read correlation file
pv_count_df = pd.read_csv(pv_count_file, sep="\t", index_col=0)
corr_count_df = pd.read_csv(corr_count_file, sep="\t", index_col=0)

gene_names = pv_count_df.iloc[:, 0]
pv_count_df = pv_count_df.iloc[:,1:]

corr_count_df = corr_count_df.iloc[:,1:]

# remove Sig 45 (artifact signature)
corr_count_df.drop("SBS45", axis=1, inplace=True)
pv_count_df.drop("SBS45", axis=1, inplace=True)


# SELECT GENES to cluster

selected_corr_df = corr_count_df[(pv_count_df < pv_th).sum(axis=1)>0]
selected_names = gene_names.loc[selected_corr_df.index]
selected_corr_df.insert(0, "Name", selected_names)


# write the correlation for the selected genes
selected_corr_df.to_csv(selected_file, sep="\t")
# END SELECT GENES

# Options for clustering

matrix_to_cluster = selected_corr_df.iloc[:,1:]
subset_genes =  list(selected_corr_df.index)
denom=10
step=5



# Clustering genes:
# Step 1: k-means clustering
labels_all = []
for ri in range(100):
    if ri%10 == 0:
        print(ri)
    clusterer_ri = MiniBatchKMeans(n_clusters=(int(ri/denom)+1)*step, random_state=ri)
    labels = clusterer_ri.fit_predict(matrix_to_cluster)
    labels_all.append(labels)

labels_df = pd.DataFrame(labels_all, columns=matrix_to_cluster.index)





# Step 2: find consensus score for each 
def count_consensus(df, label):
    return (df == label).astype(int).sum(1)

consensus_df = pd.DataFrame()
labels_T = labels_df.T
for i in range(labels_df.shape[1]):
    if (i%1000==0): # print progress
        print(i)
    consensus_df[i] = count_consensus(labels_T, labels_T.iloc[i,:])

consensus_df.columns=consensus_df.index
consensus_df.to_csv(corr_consensus_file, sep="\t", compression="gzip")


# Step 3: hierarchical clustering
consensus_mat = consensus_df.values
link_matrix = sch.linkage(consensus_mat, method="complete", metric="cosine") 
cluster_id_df = pd.DataFrame(index=matrix_to_cluster.index)
nrange = range(5, 21)
for n_clusters in nrange:
    print(n_clusters)
    th = [x[2] for x in link_matrix][-n_clusters]
    labels = sch.fcluster(link_matrix, th, criterion='distance')
    cluster_id_df["k="+str(n_clusters)] = labels


cluster_id_df.to_csv(cluster_id_file, sep="\t")

