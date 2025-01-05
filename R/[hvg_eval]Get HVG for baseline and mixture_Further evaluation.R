
if(0){
    remove.packages("mixhvg")
    devtools::install_github("RuzhangZhao/mixhvg",force=T)
    remove.packages("benchmarkHVG")
    devtools::install_github("RuzhangZhao/benchmarkHVG",force=T)
}
args <- commandArgs(trailingOnly = TRUE)
print(args)
datatype = args[1]
type0 = 'base'
type1 = 'mix'
id = as.numeric(args[2])
print(id)
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

if(datatype == 'sort'){
    dataset_name_range = c("duo4",
                           "duo4un",
                           "duo8",
                           "GBM_sd",
                           "human",
                           "mouse",
                           "zheng")[id]
    set.seed(1)
    
    ## cell sorting
    for (dataset_name in dataset_name_range) {
        cat(dataset_name,"\n")
        if (dataset_name == "duo4") {
            duo_path = "/dcs04/hongkai/data/rzhao/pbmc/duo/"
            expr<-readRDS(paste0(duo_path,"duo4_expr.rds"))
            cell_label<-readRDS(paste0(duo_path,"duo4_label.rds"))
        } else if (dataset_name == "duo4un") {
            duo_path = "/dcs04/hongkai/data/rzhao/pbmc/duo/"
            expr<-readRDS(paste0(duo_path,"duo4un_expr.rds"))
            cell_label<-readRDS(paste0(duo_path,"duo4un_label.rds"))
        } else if (dataset_name == "duo8") {
            duo_path = "/dcs04/hongkai/data/rzhao/pbmc/duo/"
            expr<-readRDS(paste0(duo_path,"duo8_expr.rds"))
            cell_label<-readRDS(paste0(duo_path,"duo8_label.rds"))
        } else if (dataset_name == "GBM_sd") {
            expr<-readRDS("/dcs04/hongkai/data/rzhao/pbmc/GBM_sd/GBM_rna.rds")
            cell_label<-readRDS("/dcs04/hongkai/data/rzhao/pbmc/GBM_sd/GBM_cell_label.rds")
        } else if (dataset_name == "human") {
            human_path<-"/dcs04/hongkai/data/rzhao/pbmc/cell_sorting/human/"
            expr<-readRDS(paste0(human_path,"human_cell_ENAR_final_500.rds"))
            cell_label<-readRDS(paste0(human_path,"human_label_ENAR_final_500_short.rds"
            ))
        } else if (dataset_name == "mouse") {
            mouse_path<-"/dcs04/hongkai/data/rzhao/pbmc/cell_sorting/mouse/"
            expr<-readRDS(paste0(mouse_path,"mouse_cell_ENAR_final.rds"))
            cell_label<-readRDS(paste0(mouse_path,"mouse_label_ENAR_final.rds"))
        } else if (dataset_name == "zheng") {
            zheng_path<-"/dcs04/hongkai/data/rzhao/pbmc/cell_sorting/zheng/"
            expr<-readRDS(paste0(zheng_path,"zheng_cells_new.rds"))
            cell_label<-readRDS(paste0(zheng_path,"zheng_label_new.rds"))
        }
        #cat("expr: ",dim(as.matrix(expr)),"\n")
        pcafile<-paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/",dataset_name,"/final_",dataset_name,"_pcalist_",type0,".rds")
        pcafile2<-paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/",dataset_name,"/final_",dataset_name,"_pcalist_",type1,".rds")
        if(!file.exists(pcafile)){
            pcalist<-hvg_pca(expr)
            saveRDS(pcalist,pcafile)
        }else{
            pcalist<-readRDS(pcafile)
        }
        
        if(!file.exists(pcafile2)){
            pcalist2<-mixture_hvg_pca(expr,method_list=c("mv_lognc","logmv_lognc","scran_pos","seuratv1","mean_max_nc"))
            saveRDS(pcalist2,pcafile2)
        }else{
            pcalist2<-readRDS(pcafile2)
        }
        
        pcalist12=c(pcalist[[1]],pcalist2[[1]])
        if(is.null(rownames(pcalist12[[1]])[1])){
            rownames(pcalist12[[1]]) = 1:nrow(pcalist12[[1]])
        }
        for(i in 1:length(pcalist12)){
            rownames(pcalist12[[i]]) = rownames(pcalist12[[1]])
        }

        ev_res<-evaluate_hvg_discrete(pcalist12,cell_label)
        saveRDS(ev_res,paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/",dataset_name,"/final_",dataset_name,"_evaluate.rds"))
        print("Done!")
    }
}
if(datatype == 'cite'){
    ## CITESEQ
    
    dataset_name_range = c("cbmc_pbmc", #1
                           "cbmc8k",#2
                           "fetalBM",#3
                           "CD34",#4
                           "bmcite",#5
                           "seurat_cite",#6
                           "Sucovid")[id]  #7
    for (dataset_name in dataset_name_range) {
        cat(dataset_name,"\n")
        if (dataset_name == "cbmc_pbmc") {
            CBMCpath<-"/dcs04/hongkai/data/rzhao/pbmc/cbmc/"
            expr<-readRDS(paste0(CBMCpath,dataset_name,"_rna_filter.rds"))
            pro<-readRDS(paste0(CBMCpath,dataset_name,"_pro_filter.rds"))
        } else if (dataset_name == "cbmc8k") {
            CBMCpath<-"/dcs04/hongkai/data/rzhao/pbmc/cbmc/"
            expr<-readRDS(paste0(CBMCpath,dataset_name,"_rna_filter.rds"))
            pro<-readRDS(paste0(CBMCpath,dataset_name,"_pro_filter.rds"))
        } else if (dataset_name == "fetalBM") {
            ADT_file<-"/dcs04/hongkai/data/rzhao/pbmc/cellatlas/GSE166895_postQC_ADT_raw_FBM-MNCs.csv.gz"
            pro<-t(read.table(file =ADT_file, header = T, row.names=1,sep=",", as.is=T))
            colnames(pro)<-paste0(colnames(pro),'-1')
            RNA_file<-"/dcs04/hongkai/data/rzhao/pbmc/cellatlas/GSE166895_postQC_mRNAraw_FBM-MNCs.csv.gz"
            expr<-t(read.table(file =RNA_file, header = T, row.names=1,sep=",", as.is=T))
            expr<-expr[,colnames(expr)%in%colnames(pro)]
            pro<-pro[,colnames(expr)]
            
        } else if (dataset_name == "CD34") {
            cellatlas_path<-"/dcs04/hongkai/data/rzhao/pbmc/cellatlas"
            expr<-readRDS(file.path(cellatlas_path,"RNA_FL-FBM-CB.rds"))
            pro<-readRDS(file.path(cellatlas_path,"ADT_FL-FBM-CB.rds"))
        } else if (dataset_name == "bmcite") {
            
            expr<-readRDS("/dcs04/hongkai/data/rzhao/pbmc/bmcite_expr.rds")
            pro<-readRDS("/dcs04/hongkai/data/rzhao/pbmc/bmcite_pro.rds")
        } else if (dataset_name == "seurat_cite") {
            expr<-readRDS("/dcs04/hongkai/data/rzhao/pbmc/seurat_cite/pbmc_multimodal_count.rds")
            pro<-readRDS("/dcs04/hongkai/data/rzhao/pbmc/seurat_cite/pbmc_multimodal_pro.rds")
        } else if (dataset_name == "Sucovid") {
            expr <- readRDS("/users/rzhao/Sucovid/count.rds")
            pro<-readRDS("/users/rzhao/Sucovid/pro.rds")
            
        }
        #cat("expr: ",dim(as.matrix(expr)),"\n")
        pcafile<-paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/",dataset_name,"/final_",dataset_name,"_pcalist_",type0,".rds")
        pcafile2<-paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/",dataset_name,"/final_",dataset_name,"_pcalist_",type1,".rds")
        pcafile<-paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/",dataset_name,"/final_",dataset_name,"_pcalist_",type0,".rds")
        pcafile2<-paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/",dataset_name,"/final_",dataset_name,"_pcalist_",type1,".rds")
        if(!file.exists(pcafile)){
            pcalist<-hvg_pca(expr)
            saveRDS(pcalist,pcafile)
        }else{
            pcalist<-readRDS(pcafile)
        }
        
        if(!file.exists(pcafile2)){
            pcalist2<-mixture_hvg_pca(expr,method_list=c("mv_lognc","logmv_lognc","scran_pos","seuratv1","mean_max_nc"))
            saveRDS(pcalist2,pcafile2)
        }else{
            pcalist2<-readRDS(pcafile2)
        }
        
        pcalist12=c(pcalist[[1]],pcalist2[[1]])
        if(is.null(rownames(pcalist12[[1]])[1])){
            rownames(pcalist12[[1]]) = 1:nrow(pcalist12[[1]])
        }
        for(i in 1:length(pcalist12)){
            rownames(pcalist12[[i]]) = rownames(pcalist12[[1]])
        }
        #}
        ev_res<-evaluate_hvg_continuous(pcalist12,pro,input = "CITEseq")#,dataset_name=dataset_name)
        saveRDS(ev_res,paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/",dataset_name,"/final_",dataset_name,"_evaluate.rds"))
    }
}

if(datatype == 'multi'){
    ## multiomics
    
    dataset_name = c("pbmc3k",
                     "human_brain_3k",
                     "mouse_brain_fresh_5k",
                     "pbmc10k",
                     "lymphoma_14k")[id]
    
    ## Multiomics
    
    cat(dataset_name,"\n")
    if (dataset_name == "pbmc3k") {
        expr<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_rna.mat.rds"))
        pro<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_lsistdev.rds"))
    } else if (dataset_name == "human_brain_3k") {
        expr<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_rna.mat.rds"))
        pro<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_lsistdev.rds"))
    } else if (dataset_name == "mouse_brain_fresh_5k") {
        expr<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_rna.mat.rds"))
        pro<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_lsistdev.rds"))
    } else if (dataset_name == "pbmc10k") {
        expr<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_rna.mat.rds"))
        pro<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_lsistdev.rds"))
    } else if (dataset_name == "lymphoma_14k") {
        expr<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_rna.mat.rds"))
        pro<-readRDS(paste0("/dcs04/hongkai/data/rzhao/multiome/10x/",dataset_name,"/",dataset_name,"_lsistdev.rds"))
    }
    
    seurat.obj6<-CreateSeuratObject(expr,verbose = F)
    seurat.obj6<-NormalizeData(seurat.obj6,verbose = F)
    seurat.obj6<-FindVariableFeaturesMix(seurat.obj6,method.names = "seuratv3",nfeatures = nfeatures)
    #seurat.obj6<-FindVariableFeatures(seurat.obj6,nfeatures=nfeatures,verbose = F)
    seurat.obj6<-ScaleData(seurat.obj6,verbose = F)
    seurat.obj6<-RunPCA(seurat.obj6,npcs=30,verbose=F)
    
    rm(seurat.obj6)
    
    pcafile<-paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/",dataset_name,"/final_",dataset_name,"_pcalist_",type0,".rds")
    pcafile2<-paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/",dataset_name,"/final_",dataset_name,"_pcalist_",type1,".rds")
    
    if(!file.exists(pcafile)){
        pcalist<-hvg_pca(expr)
        saveRDS(pcalist,pcafile)
    }else{
        pcalist<-readRDS(pcafile)
    }
    
    if(!file.exists(pcafile2)){
        pcalist2<-mixture_hvg_pca(expr,method_list=c("mv_lognc","logmv_lognc","scran_pos","seuratv1","mean_max_nc"))
        saveRDS(pcalist2,pcafile2)
    }else{
        pcalist2<-readRDS(pcafile2)
    }
    
    
    pcalist12=c(pcalist[[1]],pcalist2[[1]])
    if(is.null(rownames(pcalist12[[1]])[1])){
        rownames(pcalist12[[1]]) = 1:nrow(pcalist12[[1]])
    }
    for(i in 1:length(pcalist12)){
        rownames(pcalist12[[i]]) = rownames(pcalist12[[1]])
    }
    #}
    ev_res<-evaluate_hvg_continuous(pcalist12,pro[-1,])#,dataset_name=dataset_name)
    saveRDS(ev_res,paste0("/dcs04/hongkai/data/rzhao/pbmc/explore_hvg/",dataset_name,"/final_",dataset_name,"_evaluate_jiuyao.rds"))
}
