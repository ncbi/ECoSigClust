#!/usr/bin/env python
# coding: utf-8


import scipy.stats as ss
import pandas as pd


n_clusters = int(snakemake.params.nc_param)

enrich_pv_th = 0.01

# input files
corr_count_file = snakemake.input[0]
selected_corr_file = snakemake.input[1]
go_file = snakemake.input[2]
cluster_id_file = snakemake.input[3]

print("reading files...")
def read_corr_file(corr_file):
    corr_df = pd.read_csv(corr_file, sep="\t", index_col=0)
    gene_names = corr_df.iloc[:, 0]
    gene_name_dic = dict(zip(corr_df.index, gene_names))  # ensembl -> gene mapping
    corr_df = corr_df.iloc[:, 1:]

    return corr_df, gene_name_dic


selected_corr_df, _ = read_corr_file(selected_corr_file)
corr_count_df, gene_name_dic = read_corr_file(corr_count_file)

cluster_id_df = pd.read_csv(cluster_id_file, sep="\t", index_col=0)


# GO genes
lines = [l.split() for l in open(go_file).readlines()]
go_dic = {}
for l in lines:
    go_dic[l[0]] = l[2:]


def my_gsea(gene_set, go_dic, background, enrich_pv_th):
    """
    given a background geneset, discover significance (pval) of a gene_set for each GO
    geneset is a subset of background
    for each GO term, compute
    1) overlap with the GO genes and background --> n
    2) overlap with the GO genes and gene set --> k
    3) size of gene set --> N
    4) size of background --> M
    """

    M = len(background)
    enriched = []
    for go in go_dic:
        go_genes = go_dic[go]
        N = len(background.intersection(gene_set))
        if N < 3:
            continue
        n = len(background.intersection(go_genes))
        k = len(set(go_genes).intersection(gene_set))

        p = ss.hypergeom.sf(k - 1, M, n, N)  # fix precision error
        common_genes = list(set(go_genes).intersection(gene_set))
        if p < enrich_pv_th:
            common_genes.sort()
            enriched.append((go, p, k, n, N,";".join(common_genes)))
            # print(enriched)
    if len(enriched) > 0:
        enriched_df = pd.DataFrame(list(zip(*enriched)), ["go", "pv", "k", "n", "N", "genes"]).T
        enriched_df.sort_values(by=["pv"], inplace=True)
        return enriched_df
    return pd.DataFrame(columns=["go", "pv", "k", "n", "N", "genes"])



cluster_ids = cluster_id_df["k="+str(n_clusters)]

# create gene_group dict: cluster id -> list of genes in the cluster
gene_df = pd.DataFrame(index=selected_corr_df.index)
gene_df["gene"] = [gene_name_dic[id] for id in selected_corr_df.index]
gene_df["label"] = cluster_ids[cluster_ids > 0]
gene_groups = {}
for (l, group_df) in gene_df.groupby("label"):
    gene_groups[l] = group_df.gene.values


gene_rows = []

# different background
all_genes = set([gene_name_dic[eid] for eid in list(corr_count_df.index)])
corr_genes = set([gene_name_dic[eid] for eid in  list(selected_corr_df.index)])
bg_genes = corr_genes


# compute enrichment for each cluster
for i in range(n_clusters):
    print(i)
    cluster_id = i+1
    my_genes = list(gene_groups[i+1])
    my_genes.sort()
    gene_rows.append(("CL"+str(i+1), ";".join(my_genes)))
    enriched_df  = my_gsea(my_genes, go_dic, bg_genes, enrich_pv_th)
    filename = snakemake.output[i]
    sys.stderr.write("%s\n" % filename)
    enriched_df.sort_values(by=["pv"], inplace=True)
    enriched_df.to_csv(filename, sep=',')





