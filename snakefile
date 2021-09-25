# parameters for clustering
start, end = 5, 21 
n_clusters_range = range(start, end) 
n_clusters_list = [str(i) for i in n_clusters_range]

# cluster IDs for enrichment files (K=9)
cluster_ids = range(1, 10) 

rule all:
    input:
        "results/LUAD_norm_covid_corr_consensus.tsv.gz",
        "results/LUAD_norm_covid_cluster_id.tsv", 
        "figs/LUAD_norm_covid_corr_heatmap_clusters.pdf", 
        "figs/LUAD_norm_covid_evaluated_clusters.pdf",
        expand("results/LUAD_cluster_genes/LUAD_norm_covid_K{n_clusters}_genes.txt",  
               n_clusters=n_clusters_list),
        expand("results/LUAD_norm_K9/LUAD_norm_covid_GO_K9_CL{cluster_id}.csv",
               cluster_id=cluster_ids)

rule clustering:
    input:
        "data/LUAD_norm_corr_expr_count.tsv",
        "data/LUAD_norm_corr_pv_expr_count.tsv"
    output:
        "results/LUAD_norm_covid_corr_consensus.tsv.gz",
        "results/LUAD_norm_covid_cluster_id.tsv",
        "results/LUAD_norm_selected_corr.tsv"
    script:
        "covid_clustering.py"


rule vis_clusters:
    input:
        "results/LUAD_norm_selected_corr.tsv",
        "results/LUAD_norm_covid_cluster_id.tsv"
    params:
        start = start,
        end = end,
    output:
        "figs/LUAD_norm_covid_corr_heatmap_clusters.pdf",
        expand("results/LUAD_cluster_genes/LUAD_norm_covid_K{n_clusters}_genes.txt",
               n_clusters=n_clusters_range)
    script:
        "covid_heatmaps.py"

rule eval_clusters:
    input:
        "results/LUAD_norm_covid_corr_consensus.tsv.gz",
        "results/LUAD_norm_covid_cluster_id.tsv"
    output:
        "figs/LUAD_norm_covid_evaluated_clusters.pdf"
    script:
        "covid_evaluate_clusters.py"


rule covid_enrich:
    input:
        "data/LUAD_norm_corr_expr_count.tsv",
        "results/LUAD_norm_selected_corr.tsv",
        "data/GO.v7.1.symbols.gmt",
        "results/LUAD_norm_covid_cluster_id.tsv"
    params:
        nc_param = 9,
    output:
        expand("results/LUAD_norm_K9/LUAD_norm_covid_GO_K9_CL{cluster_id}.csv",
               cluster_id=cluster_ids),
    script:
        "covid_enrich.py"


