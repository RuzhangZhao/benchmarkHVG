args = commandArgs(trailingOnly = TRUE)
ratioid = as.numeric(args[1]) 
labelid = as.numeric(args[2]) 


library(Seurat)
library(SeuratObject)
library(mixhvg)
library(Seurat)
library(reticulate)
library(Rfast)
library(FNN)
library(pdist)
library(caret)
library(mclust)
library(cluster)
library(NMI)
library(lisi)

if(0){
    zheng_path<-"/dcs04/hongkai/data/rzhao/pbmc/cell_sorting/zheng/"
    expr<-readRDS(paste0(zheng_path,"zheng_cells_new.rds"))
    colnames(expr)<-make.names(colnames(expr), unique = TRUE)
    cell_label<-readRDS(paste0(zheng_path,"zheng_label_new.rds"))
    
    dataset_name = "cbmc8k"
    CBMCpath<-"/dcs04/hongkai/data/rzhao/pbmc/cbmc/"
    expr<-readRDS(paste0(CBMCpath,dataset_name,"_rna_filter.rds"))
    cell_label<-readRDS("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/jiuyaoplot/rare_cbmc8k_pro_cluster.rds")
    
    dataset_name = "pbmc10k"
    expr<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_rna.mat.rds"))
    cell_label<-readRDS("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/jiuyaoplot/rare_pbmc10k_pro_cluster.rds")
    
    table(cell_label)
    dim(expr)
    length(cell_label)
    
}
dataset_name = "pbmc10k"
expr<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_rna.mat.rds"))
cell_label<-readRDS("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/jiuyaoplot/rare_pbmc10k_pro_cluster.rds")

rownames(expr)<-1:nrow(expr)
colnames(expr)<-1:ncol(expr)
unique_cell_labels = c("B cells","Monocytes","CD34+","NK cells",
"Cytotoxic","CD4+ Memory","Naive Cytotoxic","CD4+ Naive T","Regulatory")
unique_cell_labels = sort(unique(cell_label))
print(unique_cell_labels)
set.seed(1)
ratio_list = c(1/2,1/4,1/8,1/16,1/32,1/64)
ratio = ratio_list[ratioid]
print(ratio)
print(labelid)
rm_index = sample(which(cell_label == unique_cell_labels[labelid]),
       round(sum(cell_label == unique_cell_labels[labelid])*(1-ratio)))
cell_label_new = cell_label[-rm_index]
expr_new = expr[,-rm_index]
cell_label_new_binary = rep("Others",length(cell_label_new))
cell_label_new_binary[cell_label_new == unique_cell_labels[labelid]] = unique_cell_labels[labelid]
saveRDS(list(cell_label_new_binary,cell_label_new),paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/","rare","/",dataset_name,"_true_cell_label_","base","_ratio",ratio,"_label",labelid,".rds"))

source("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/code/func24.R")
pcafile<-paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/","rare","/",dataset_name,"_pcalist_","base","_ratio",ratio,"_label",labelid,".rds")
if(file.exists(pcafile)){
    pcalist<-readRDS(pcafile)
}else{
    pcalist<-hvg_pca(expr_new)
    saveRDS(pcalist,pcafile)   
}
pcafile2<-paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/","rare","/",dataset_name,"_pcalist_","mix","_ratio",ratio,"_label",labelid,".rds")
if(file.exists(pcafile2)){
    pcalist2<-readRDS(pcafile2)    
}else{
    seurat.obj.pca<-list()
    var.seurat.obj<-list()
    
    seurat.obj1<-CreateSeuratObject(expr_new)
    seurat.obj1<-NormalizeData(seurat.obj1,verbose = F)
    seurat.obj1<-FindVariableFeaturesMix(seurat.obj1,nfeatures = 2000,verbose = F)
    seurat.obj1<-ScaleData(seurat.obj1,verbose = F)
    seurat.obj1<-RunPCA(seurat.obj1,npcs=30,verbose=F)
    
    seurat.obj.pca[[1]]<-seurat.obj1@reductions$pca@cell.embeddings
    var.seurat.obj[[1]]<-VariableFeatures(seurat.obj1)
    
    pcalist2<-list("seurat.obj.pca"=seurat.obj.pca,
                   "var.seurat.obj"=var.seurat.obj)
    
    #pcalist2<-#mixture_hvg_pca(expr_new,method_list=c("mv_lognc","logmv_lognc","scran_pos","seuratv1","mean_max_nc"))
    saveRDS(pcalist2,pcafile2)   
}

pcalist12=c(pcalist[[1]],pcalist2[[1]])
if(is.null(rownames(pcalist12[[1]])[1])){
    rownames(pcalist12[[1]]) = 1:nrow(pcalist12[[1]])
}
for(i in 1:length(pcalist12)){
    rownames(pcalist12[[i]]) = rownames(pcalist12[[1]])
}

filename = paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/","rare","/",dataset_name,"_eval_","base","_ratio",ratio,"_label",labelid,".rds")
if(file.exists(filename)){
    res<-readRDS(filename)   
}else{
    res<-list()
}
if(length(res) < length(pcalist12)){
    for(i in (length(res)+1):length(pcalist12)){
        res[[i]]<-ARI_NMI_func_Max(pcalist12[[i]],cell_label_new_binary)  
        saveRDS(res,filename)
    }
    
}
### section 2 giniclust

ARI_NMI_func_Max<-function(
        embedding,
        cell_label
){
    #cell_label = as.numeric(as.factor(cell_label))
    unique_names = unique(cell_label)
    unique_names = unique_names[unique_names!='Others']
    snn_<- FindNeighbors(object = embedding,
                         nn.method = "rann",
                         verbose = F)$snn
    Nresolution=20
    resolution_range=seq(0,2,0.1)
    res<-sapply(resolution_range,function(cur_resolution){
        cluster_label <- FindClusters(snn_,
                                      resolution = cur_resolution,
                                      verbose = F)[[1]]
        x<-sapply(unique(cluster_label), function(i){
            cluster_cell_label = rep("Others",length(cell_label))
            cluster_cell_label[cluster_label == i] = unique_names
            compute_f1(cell_label,cluster_cell_label,unique_names)
        })
        max(x)
    })
    res
}


compute_f1 <- function(true_labels, predictions, positive_label) {
    conf_matrix <- table(Predicted = predictions, Actual = true_labels)
    TP <- conf_matrix[positive_label, positive_label]
    FP <- sum(conf_matrix[positive_label, ]) - TP
    FN <- sum(conf_matrix[, positive_label]) - TP
    
    precision <- TP / (TP + FP)
    recall <- TP / (TP + FN)
    f1_score <- 2 * (precision * recall) / max(precision + recall,1e-10)
    
    return(f1_score)
}

sc = import("scanpy")
gini = import("giniclust3")
sp = import("scipy.sparse")
#seurat.obj16<-CreateSeuratObject(rna_mat)
#scanpy_rna=sc$AnnData(sp$csc_matrix(Matrix::t(expr_new)))
#expr_new = expr_new[rowSums(expr_new) != 0,]
scanpy_rna=sc$AnnData(as.matrix(Matrix::t(expr_new)))
scanpy_rna$obs_names = as.character(1:ncol(expr_new))
scanpy_rna$var_names = rownames(expr_new)
#sc$pp$filter_cells(scanpy_rna, min_genes=200)
sc$pp$filter_genes(scanpy_rna, min_cells=1)
sc$pp$normalize_total(scanpy_rna, target_sum=1e4)
gini$gini$calGini(scanpy_rna) 
print(sum(scanpy_rna$var$gini))
if(sum(scanpy_rna$var$gini) != 0 ){
    adataGini = scanpy_rna$X[,scanpy_rna$var$gini]
    scaleMatrix0=gini$gini$arctanTransform(adataGini)
    rownames(scaleMatrix0)<-colnames(expr_new)   
}else{
    scaleMatrix0 = NA
}

ap = 0.0001
bp = 1
pp = 0.5
gini$gini$calGini(scanpy_rna,p_value = pp) 
count = 0
while((abs(sum(scanpy_rna$var$gini) - 2000) > 100)&(count<10) ){
    gini$gini$calGini(scanpy_rna,p_value = pp) 
    print(sum(scanpy_rna$var$gini))
    if (sum(scanpy_rna$var$gini) > 2000){
        bp = pp
    }else{
        ap = pp
    }
    pp = (ap+bp)/2
    count = count+1
}
print(sum(scanpy_rna$var$gini))
seurat.obj15<-CreateSeuratObject(expr_new)    
seurat.obj15<-NormalizeData(seurat.obj15,verbose = F)
VariableFeatures(seurat.obj15)<-rownames(expr_new)[scanpy_rna$var$gini]
seurat.obj15<-ScaleData(seurat.obj15,verbose = F)
seurat.obj15<-RunPCA(seurat.obj15,npcs=30,verbose=F)

scaleMatrix<-seurat.obj15@reductions$pca@cell.embeddings

xx0<-ARI_NMI_func_Max(scaleMatrix0,cell_label_new_binary)
xx<-ARI_NMI_func_Max(scaleMatrix,cell_label_new_binary)

ginipcafile<-paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/","rare","/",dataset_name,"_giniclust_","base","_ratio",ratio,"_label",labelid,".rds")
ginievafile<-paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/","rare","/",dataset_name,"_giniclust_evaluate_","base","_ratio",ratio,"_label",labelid,".rds")
saveRDS(list(scaleMatrix0,scaleMatrix),ginipcafile)
saveRDS(list(xx0,xx),ginipcafile)



       