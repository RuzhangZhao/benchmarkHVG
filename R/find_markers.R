### Get marker genes for each dataset
# find_markers.R: find marker genes for each dataset using different min.pct and logfc.threshold values

library(Seurat)
library(parallel)

# Get marker genes given clusters -----------------------------------------

dataset_name_range = c("duo4",
                       "duo4un",
                       "duo8",
                       "GBM_sd",
                       "human",
                       "mouse",
                       "zheng",
                       "cbmc_pbmc", #1
                       "cbmc8k",#2
                       "fetalBM",#3
                       "CD34",#4
                       "bmcite",#5
                       "seurat_cite",#6
                       "Sucovid",
                       "pbmc3k",
                       "human_brain_3k",
                       "mouse_brain_fresh_5k",
                       "pbmc10k",
                       "lymphoma_14k")

get_marker <- function(expr, cell_label, min.pct = 0.1, logfc.threshold = 0.25) {
  # Normalize
  obj<-CreateSeuratObject(expr)
  obj$cell_label<-cell_label
  Idents(obj)<-obj$cell_label
  obj<-NormalizeData(obj,verbose=F)

  # Find markers
  FindAllMarkers(
    object = obj,
    only.pos = TRUE,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    test.use = "wilcox"
  )
}

get_expr <- function(dataset_name) {
  if (dataset_name == "duo4") {
    duo_path = "/dcs04/hongkai/data/rzhao/pbmc/duo/"
    expr<-readRDS(paste0(duo_path,"duo4_expr.rds"))
  } else if (dataset_name == "duo4un") {
    duo_path = "/dcs04/hongkai/data/rzhao/pbmc/duo/"
    expr<-readRDS(paste0(duo_path,"duo4un_expr.rds"))
  } else if (dataset_name == "duo8") {
    duo_path = "/dcs04/hongkai/data/rzhao/pbmc/duo/"
    expr<-readRDS(paste0(duo_path,"duo8_expr.rds"))
  } else if (dataset_name == "GBM_sd") {
    expr<-readRDS("/dcs04/hongkai/data/rzhao/pbmc/GBM_sd/GBM_rna.rds")
  } else if (dataset_name == "human") {
    human_path<-"/dcs04/hongkai/data/rzhao/pbmc/cell_sorting/human/"
    expr<-readRDS(paste0(human_path,"human_cell_ENAR_final_500.rds"))
  } else if (dataset_name == "mouse") {
    mouse_path<-"/dcs04/hongkai/data/rzhao/pbmc/cell_sorting/mouse/"
    expr<-readRDS(paste0(mouse_path,"mouse_cell_ENAR_final.rds"))
  } else if (dataset_name == "zheng") {
    zheng_path<-"/dcs04/hongkai/data/rzhao/pbmc/cell_sorting/zheng/"
    expr<-readRDS(paste0(zheng_path,"zheng_cells_new.rds"))
  } else if (dataset_name == "cbmc_pbmc") {
    CBMCpath<-"/dcs04/hongkai/data/rzhao/pbmc/cbmc/"
    expr<-readRDS(paste0(CBMCpath,dataset_name,"_rna_filter.rds"))
  } else if (dataset_name == "cbmc8k") {
    CBMCpath<-"/dcs04/hongkai/data/rzhao/pbmc/cbmc/"
    expr<-readRDS(paste0(CBMCpath,dataset_name,"_rna_filter.rds"))
  } else if (dataset_name == "fetalBM") {
    ADT_file<-"/dcs04/hongkai/data/rzhao/pbmc/cellatlas/GSE166895_postQC_ADT_raw_FBM-MNCs.csv.gz"
    pro<-t(read.table(file =ADT_file, header = T, row.names=1,sep=",", as.is=T))
    colnames(pro)<-paste0(colnames(pro),'-1')
    RNA_file<-"/dcs04/hongkai/data/rzhao/pbmc/cellatlas/GSE166895_postQC_mRNAraw_FBM-MNCs.csv.gz"
    expr<-t(read.table(file =RNA_file, header = T, row.names=1,sep=",", as.is=T))
    expr<-expr[,colnames(expr)%in%colnames(pro)]
  } else if (dataset_name == "CD34") {
    cellatlas_path<-"/dcs04/hongkai/data/rzhao/pbmc/cellatlas"
    expr<-readRDS(file.path(cellatlas_path,"RNA_FL-FBM-CB.rds"))
  } else if (dataset_name == "bmcite") {
    expr<-readRDS("/dcs04/hongkai/data/rzhao/pbmc/bmcite_expr.rds")
  } else if (dataset_name == "seurat_cite") {
    expr<-readRDS("/dcs04/hongkai/data/rzhao/pbmc/seurat_cite/pbmc_multimodal_count.rds")
  } else if (dataset_name == "Sucovid") {
    expr <- readRDS("/dcs04/hongkai/data/rzhao/Sucovid/count.rds")
  } else if (dataset_name == "pbmc3k") {
    expr<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_rna.mat.rds"))
  } else if (dataset_name == "human_brain_3k") {
    expr<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_rna.mat.rds"))
  } else if (dataset_name == "mouse_brain_fresh_5k") {
    expr<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_rna.mat.rds"))
  } else if (dataset_name == "pbmc10k") {
    expr<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_rna.mat.rds"))
  } else if (dataset_name == "lymphoma_14k") {
    expr<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_rna.mat.rds"))
  }
  expr
}

get_cell_label <- function(dataset_name, resol = NULL) {
  clu_dir <- "/dcs05/hongkai/data/yli6/rzhao/mixhvg/cluster"
  if (dataset_name == "duo4") {
    duo_path = "/dcs04/hongkai/data/rzhao/pbmc/duo/"
    cell_label<-readRDS(paste0(duo_path,"duo4_label.rds"))
  } else if (dataset_name == "duo4un") {
    duo_path = "/dcs04/hongkai/data/rzhao/pbmc/duo/"
    cell_label<-readRDS(paste0(duo_path,"duo4un_label.rds"))
  } else if (dataset_name == "duo8") {
    duo_path = "/dcs04/hongkai/data/rzhao/pbmc/duo/"
    cell_label<-readRDS(paste0(duo_path,"duo8_label.rds"))
  } else if (dataset_name == "GBM_sd") {
    cell_label<-readRDS("/dcs04/hongkai/data/rzhao/pbmc/GBM_sd/GBM_cell_label.rds")
  } else if (dataset_name == "human") {
    human_path<-"/dcs04/hongkai/data/rzhao/pbmc/cell_sorting/human/"
    cell_label<-readRDS(paste0(human_path,"human_label_ENAR_final_500_short.rds"
    ))
  } else if (dataset_name == "mouse") {
    mouse_path<-"/dcs04/hongkai/data/rzhao/pbmc/cell_sorting/mouse/"
    cell_label<-readRDS(paste0(mouse_path,"mouse_label_ENAR_final.rds"))
  } else if (dataset_name == "zheng") {
    zheng_path<-"/dcs04/hongkai/data/rzhao/pbmc/cell_sorting/zheng/"
    cell_label<-readRDS(paste0(zheng_path,"zheng_label_new.rds"))
  } else {
    cell_label <- readRDS(file.path(clu_dir, paste0(dataset_name, "_clusters.rds")))[, as.character(resol)]
  }
  return(cell_label)
}

min.pct <- c(0, 0.05, 0.1)
logfc.threshold <- c(0, 0.1, 0.25)
params <- expand.grid(min.pct, logfc.threshold)

major_resol <- readRDS("/dcs05/hongkai/data/yli6/rzhao/mixhvg/major_resol.rds")

mclapply(1:nrow(params), function(i) {
  print(params[i, 1])
  print(params[i, 2])
  marker_list <- lapply(dataset_name_range, function(dataset) {
    print(dataset)
    expr <- get_expr(dataset)
    if (dataset == "Sucovid") {
      # Subsample dataset
      subsamp_idx <- readRDS("/dcs05/hongkai/data/yli6/rzhao/mixhvg/sucovid_clu_subsample_idx.rds")
      expr <- expr[, subsamp_idx]
    }
    print("Looking for cell label")
    cell_label <- get_cell_label(dataset,
                                 resol = ifelse(dataset %in% names(major_resol),
                                                major_resol[dataset], NULL))
    colnames(expr)<-1:ncol(expr)
    names(cell_label)<-1:ncol(expr)
    print("Finding marker genes")
    get_marker(expr, cell_label, params[i, 1], params[i, 2])
  })
  names(marker_list) <- dataset_name_range
  print(str(marker_list))
  saveRDS(marker_list, paste0("/dcs05/hongkai/data/yli6/rzhao/mixhvg/marker/pct",
                              params[i, 1], "_logfc", params[i, 2], ".rds"))
})

