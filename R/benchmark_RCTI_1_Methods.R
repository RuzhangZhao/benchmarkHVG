rcti_pca <- function(rna_mat,
                   nfeatures = 2000){
  
  outputs <-list()
  outputs_types <-list()
  
  ######## 1. CellSIUS
  
  library(CellSIUS)
  library(SingleCellExperiment)
  library(scuttle)
  library(scran)
  sce <- SingleCellExperiment(list(counts=rna_mat))
  clusters <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters=clusters)
  sce <- logNormCounts(sce)
  normalized_data_scran <- logcounts(sce)
  
  seurat.obj <- as.Seurat(sce, counts = "counts", data = "logcounts")
  # seurat.obj<-CreateSeuratObject(rna_mat)
  # seurat.obj<-NormalizeData(seurat.obj,verbose = F)
  seurat.obj<-FindVariableFeaturesMix(seurat.obj,method.names = "scran",nfeatures = nfeatures)
  seurat.obj<-ScaleData(seurat.obj,verbose = F)
  seurat.obj<-RunPCA(seurat.obj,npcs=30,verbose=F)
  snn_ <- FindNeighbors(
    object = seurat.obj@reductions$pca@cell.embeddings,
    nn.method = "rann",
    verbose = F
  )$snn
  cluster_label <- FindClusters(snn_,
                                verbose = F
  )[[1]]
  
  names(cluster_label) = colnames(rna_mat)
  
  
  
  source("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/rarecode/cellsius/myCellSIUS.R")
  CellSIUS.out <- try({myCellSIUS(mat.norm = as.matrix(normalized_data_scran),group_id = cluster_label)},
                      silent = TRUE)
  if (!"data.table" %in% class(CellSIUS.out)) {
    Final_Clusters = cluster_label
  } else {
    Final_Clusters = CellSIUS_final_cluster_assignment(CellSIUS.out=CellSIUS.out,
                                                       group_id=cluster_label)
  }
  
  outputs[[1]] <- Final_Clusters
  outputs_types[[1]] <- "label"
  
  
  ######## 2. GiniClust3
  
  library(reticulate)
  sc <- import("scanpy")
  gini <- import("giniclust3")
  sp <- import("scipy.sparse")
  # seurat.obj16<-CreateSeuratObject(rna_mat)
  # scanpy_rna=sc$AnnData(sp$csc_matrix(Matrix::t(rna_mat)))
  scanpy_rna <- sc$AnnData(as.matrix(Matrix::t(rna_mat)))
  scanpy_rna$obs_names <- as.character(1:ncol(rna_mat))
  scanpy_rna$var_names <- rownames(rna_mat)
  # sc$pp$filter_cells(scanpy_rna, min_genes=200)
  sc$pp$filter_genes(scanpy_rna, min_cells = 1)
  if (dataset_name %in% c("cbmc_pbmc","cbmc8k")) {
    sc$pp$normalize_total(scanpy_rna, target_sum = 1e4)
  } else {
    sc$pp$normalize_per_cell(scanpy_rna, counts_per_cell_after = 1e4)
  }
  
  ap <- 0.0001
  bp <- 1
  pp <- 0.5
  gini$gini$calGini(scanpy_rna, p_value = pp)
  count <- 0
  while ((abs(sum(scanpy_rna$var$gini) - nfeatures) > 0) & (count < 30)) {
    gini$gini$calGini(scanpy_rna, p_value = pp)
    print(sum(scanpy_rna$var$gini))
    if (sum(scanpy_rna$var$gini) > nfeatures) {
      bp <- pp
    } else {
      ap <- pp
    }
    pp <- (ap + bp) / 2
    count <- count + 1
  }
  seurat.obj <- CreateSeuratObject(rna_mat)
  seurat.obj <- NormalizeData(seurat.obj, verbose = F)
  VariableFeatures(seurat.obj) <- scanpy_rna$var_names$to_list()[scanpy_rna$var$gini]
  seurat.obj <- ScaleData(seurat.obj, verbose = F)
  seurat.obj <- RunPCA(seurat.obj, npcs = 30, verbose = F)
  
  scaleMatrix <- seurat.obj@reductions$pca@cell.embeddings
  
  outputs[[2]] <- scaleMatrix
  outputs_types[[2]] <- "embedding"
  
  
  ######## 3. SCA
  
  library(reticulate)
  sca = import("shannonca")
  reduction = sca$dimred$reduce(t(as.matrix(rna_mat)), n_comps=as.integer(30), iters = 5)
  
  embedding = reduction
  rownames(embedding) = 1:nrow(reduction)
  colnames(embedding) = 1:ncol(reduction)
  
  outputs[[3]] <- embedding
  outputs_types[[3]] <- "embedding"
  
  
  ######## 4. EDGE
  
  library(EDGE)
  source("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/rarecode/edge30/myEDGE.R")
  
  seurat.obj <- CreateSeuratObject(rna_mat)
  seurat.obj <- NormalizeData(seurat.obj,verbose = F)
  dat_norm <- GetAssayData(seurat.obj, slot = "data")
  mat_norm <- as.matrix(dat_norm)
  dat_log <- t(log2(mat_norm+1))
  
  custom_defs <- endr_defs
  custom_defs$n_comps <- 30
  res_endr = myendr(dat_log, custom_defs)
  
  embedding = as.matrix(res_endr)
  rownames(embedding) = 1:nrow(res_endr)
  colnames(embedding) = 1:ncol(res_endr)
  
  outputs[[4]] <- embedding
  outputs_types[[4]] <- "embedding"
  
  
  ######## 5. HIG
  
  library(reticulate)
  use_condaenv("scCAD_env", required = TRUE)
  scCAD <- import("scCAD")
  np = import("numpy")
  py_float <- py_eval("float")
  data <- np$vectorize(py_float)(t(as.matrix(rna_mat)))
  geneNames = np$array(rownames(rna_mat))
  cellNames = np$array(colnames(rna_mat))
  HIGres = scCAD$myscCAD$getHIG(data=data, 
                                dataName=dataset_name,
                                cellNames=cellNames, geneNames=geneNames) 
  
  selected_gene_names = HIGres[[1]]
  selected_gene_imp = HIGres[[2]]
  
  seurat.obj <- CreateSeuratObject(rna_mat)
  seurat.obj <- NormalizeData(seurat.obj, verbose = F)
  VariableFeatures(seurat.obj) <- selected_gene_names[order(selected_gene_imp, decreasing = TRUE)[1:nfeatures]]
  seurat.obj <- ScaleData(seurat.obj, verbose = F)
  seurat.obj <- RunPCA(seurat.obj, npcs = 30, verbose = F)
  scaleMatrix <- seurat.obj@reductions$pca@cell.embeddings
  
  outputs[[5]] <- scaleMatrix
  outputs_types[[5]] <- "embedding"
  
  
  list("outputs" = outputs,
       "outputs_types" = outputs_types)
}

