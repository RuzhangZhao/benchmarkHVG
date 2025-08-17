# source("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/rarecode/cellsius/myCellSIUS.R")
myCellSIUS = function (mat.norm, group_id, min_n_cells = 10, min_fc = 2, corr_cutoff = NULL, 
                       iter = 0, max_perc_cells = 50, fc_between_cutoff = 1, 
                       mcl_path = "/users/jlu2/.conda/envs/CellSIUS_env/bin/mcl") 
{
  options(warn = -1)
  # if (class(mat.norm) != "matrix") 
  #   stop("mat.norm is not of 'matrix' class")
  if (length(unique(group_id)) < 2) 
    stop("The number of clusters in 'group_id' must be > 1")
  if (!file.exists(mcl_path)) 
    stop("I don't find mcl executable. External UNIX installation of MCL is required: https://micans.org/mcl/")
  if (is.null(rownames(mat.norm)) | is.null(colnames(mat.norm))) 
    stop("Column and row names of mat.norm matrix must be non-emptys")
  group_id = data.table::data.table(cell_idx = colnames(mat.norm), 
                                    main_cluster = as.character(group_id))
  expr_dt = data.table::data.table(gene_id = rownames(mat.norm), 
                                   mat.norm)
  expr_dt_melt = data.table::melt(expr_dt, id.vars = "gene_id", 
                                  value.name = "expr", variable.name = "cell_idx")
  expr_dt_melt = merge(expr_dt_melt, group_id, by = "cell_idx")
  expr_dt_melt[, `:=`(c("N_cells", "within_p", "pos0", "pos1", 
                        "Dpos"), CellSIUS:::cellsius_find_bimodal_genes(expr, min_n_cells = min_n_cells, 
                                                             max_perc_cells = max_perc_cells)), by = c("gene_id", 
                                                                                                       "main_cluster")]
  expr_dt_melt[, `:=`(sig, within_p < 100 & Dpos > min_fc)]
  expr_dt_melt[sig == T, `:=`(within_adj_p, p.adjust(within_p)), 
               by = c("cell_idx")]
  expr_dt_melt[, `:=`(sig, within_adj_p < 0.1)]
  expr_dt_melt = expr_dt_melt[gene_id %in% expr_dt_melt[!is.na(sig) & 
                                                          sig == T]$gene_id]
  if (dim(expr_dt_melt)[1] == 0) {
    print("No genes with bimodal distribution found, returning NA.")
    return(NA)
  }
  for (clust in unique(expr_dt_melt$main_cluster)) {
    expr_dt_melt = expr_dt_melt[, `:=`(paste0(clust, "_", 
                                              c("p_between", "fc")), CellSIUS:::cellsius_test_cluster_specificity(expr, 
                                                                                                       main_cluster, clust, fc_between_cutoff = fc_between_cutoff)), 
                                by = "gene_id"]
    expr_dt_melt[main_cluster == clust, `:=`(keep, (expr_dt_melt[main_cluster == 
                                                                   clust][[paste0(clust, "_p_between")]] < 0.1))]
  }
  expr_dt_melt = expr_dt_melt[keep == TRUE & !is.na(sig)]
  expr_dt_melt[, `:=`(n_clust_per_gene, length(unique(main_cluster))), 
               by = "gene_id"]
  expr_dt_melt = expr_dt_melt[n_clust_per_gene == 1]
  expr_dt_melt[, `:=`(n_clust_per_gene, NULL)]
  expr_dt_melt = expr_dt_melt[, `:=`(gene_cluster, 0)]
  expr_dt_melt = CellSIUS:::cellsius_find_gene_sets(expr_dt_melt, corr_cutoff = corr_cutoff, 
                                         mcl_path = mcl_path)
  expr_dt_melt = expr_dt_melt[gene_cluster != 0]
  if (dim(expr_dt_melt)[1] == 0) {
    print("No subclusters found, returning NA.")
    return(NA)
  }
  expr_dt_melt[, `:=`(sub_cluster, main_cluster)]
  expr_dt_melt[, `:=`(mean_expr, mean(expr)), by = c("main_cluster", 
                                                     "gene_cluster", "cell_idx")]
  expr_dt_melt[, `:=`(sub_cluster, CellSIUS:::cellsius_sub_cluster(mean_expr, 
                                                        sub_cluster, gene_cluster, iter = iter)), by = c("main_cluster", 
                                                                                                         "gene_cluster")]
  clust_list = expr_dt_melt[, list(sub = length(unique(cell_idx))), 
                            by = c("sub_cluster", "main_cluster")]
  clust_list[, `:=`(tot, sum(sub)/(length(sub_cluster)/2)), 
             by = "main_cluster"]
  clust_list = clust_list[grep("_1$", sub_cluster)]
  clust_list[, `:=`(perc, sub/tot * 100)]
  discard_sub_clust = clust_list[perc > max_perc_cells]$sub_cluster
  discard_sub_clust = append(discard_sub_clust, gsub("_1$", 
                                                     "_0", discard_sub_clust))
  expr_dt_melt = expr_dt_melt[!sub_cluster %in% discard_sub_clust]
  keep.columns = c("cell_idx", "gene_id", "expr", "main_cluster", 
                   "N_cells", "Dpos", "sub_cluster", "gene_cluster")
  expr_dt_melt = expr_dt_melt[, ..keep.columns]
  data.table::setnames(expr_dt_melt, "Dpos", "log2FC")
  return(expr_dt_melt)
}