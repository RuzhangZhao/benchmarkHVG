# Author: Ruzhang Zhao
# Dataset loading can be found  under
## https://github.com/RuzhangZhao/scRNAseqProcess/blob/main/DatatsetLoading.R

#' hvg_pca
#'
#' @details
#'
#' @param rna_mat input scRNA-seq count matrix
#' @param nfeatures Number of features to select as top variable features.
#'
#'
#' @return seurat.obj.pca:  a list of PCA.
#' @return var.seurat.obj:
#'
#' @import Matrix
#' @import scran
#' @import Seurat
#' @import FNN
#' @import Rfast
#' @importFrom scran modelGeneVar modelGeneVarByPoisson
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom scuttle logNormCounts
#' @export
#'
#' @examples
#' if(0){
#' simple_matrix<-matrix(1:2e4,nrow=4000,ncol=5)
#' rownames(simple_matrix)<-1:nrow(simple_matrix)
#' colnames(simple_matrix)<-1:ncol(simple_matrix)
#' simple_matrix_HVG<-FindVariableFeaturesMix(simple_matrix)
#' }
#'
hvg_pca<-function(rna_mat,
                  nfeatures = 2000){


    rna_mat_PFlog1pPF<-NormalizeData(rna_mat,scale.factor=mean(Matrix::colSums(rna_mat)),verbose=F)
    rna_mat_PFlog1pPF<-NormalizeData(rna_mat_PFlog1pPF,scale.factor=mean(Matrix::colSums(rna_mat_PFlog1pPF)),normalization.method = "RC",verbose=F)

    #seurat.obj0<-CreateSeuratObject(rna_mat)
    #rna_mat_PFlog1pPF<-t(t(rna_mat)/colSums(rna_mat))*mean(colSums(rna_mat))
    #rna_mat_PFlog1pPF<-log1p(rna_mat_PFlog1pPF)
    #rna_mat_PFlog1pPF<-t(t(rna_mat_PFlog1pPF)/colSums(rna_mat_PFlog1pPF))*mean(colSums(rna_mat_PFlog1pPF))

    seurat.obj.pca<-list()
    var.seurat.obj<-list()

    #1.random
    message("Method1: random")

    hvg<-sample(rownames(rna_mat)[which(rowSums(rna_mat)>0)],nfeatures)

    seurat.obj1<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj1<-NormalizeData(seurat.obj1,verbose = F)
    VariableFeatures(seurat.obj1)<-hvg
    seurat.obj1<-ScaleData(seurat.obj1,verbose = F)
    seurat.obj1<-RunPCA(seurat.obj1,npcs=30,verbose=F)

    seurat.obj.pca[[1]]<-seurat.obj1@reductions$pca@cell.embeddings
    var.seurat.obj[[1]]<-VariableFeatures(seurat.obj1)


    ######## 2. Most original var count ~ mean count
    message("Method2: mv_ct")
    seurat.obj2<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj2<-NormalizeData(seurat.obj2,verbose = F)

    #sce <- SingleCellExperiment(list(counts=rna_mat))
    #sce@assays@data$logcounts<-sce@assays@data$counts
    #dec <- modelGeneVar(sce)
    #top.hvgs <- getTopHVGs(dec, n=2000)
    #rm(sce)
    #VariableFeatures(seurat.obj2)<-top.hvgs
    seurat.obj2<-FindVariableFeaturesMix(seurat.obj2,method.names = "mv_ct",nfeatures = nfeatures)
    seurat.obj2<-ScaleData(seurat.obj2,verbose = F)
    seurat.obj2<-RunPCA(seurat.obj2,npcs=30,verbose=F)

    seurat.obj.pca[[2]]<-seurat.obj2@reductions$pca@cell.embeddings
    var.seurat.obj[[2]]<-VariableFeatures(seurat.obj2)
    rm(seurat.obj2)


    #################################################
    #3.original normalized selection: var(normalized)~mean(normalized)
    message("Method3: mv_nc")
    seurat.obj3<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj3<-NormalizeData(seurat.obj3,verbose = F)

    #rna_mat_norm<-seurat.obj3@assays$RNA@data
    #rna_mat_norm@x<-exp(rna_mat_norm@x)-1
    #sce <- SingleCellExperiment(list(counts=rna_mat_norm))
    #sce@assays@data$logcounts<-sce@assays@data$counts
    #dec <- modelGeneVar(sce)
    #top.hvgs <- getTopHVGs(dec, n=2000)
    #rm(sce)
    #VariableFeatures(seurat.obj3)<-top.hvgs

    seurat.obj3<-FindVariableFeaturesMix(seurat.obj3,method.names = "mv_nc",nfeatures = nfeatures)
    seurat.obj3<-ScaleData(seurat.obj3,verbose = F)
    seurat.obj3<-RunPCA(seurat.obj3,npcs=30,verbose=F)

    seurat.obj.pca[[3]]<-seurat.obj3@reductions$pca@cell.embeddings
    var.seurat.obj[[3]]<-VariableFeatures(seurat.obj3)
    rm(seurat.obj3)


    #################################################
    #4.scran selection: var(log1p(normalized))~ mean(log1p(normalized))
    message("Method4: scran")
    ## scran

    seurat.obj4<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj4<-NormalizeData(seurat.obj4,verbose = F)


    #library(scran)
    #sce <- SingleCellExperiment(list(counts=rna_mat))
    #sce <- logNormCounts(sce)
    #dec <- modelGeneVar(sce)
    #top.hvgs2 <- getTopHVGs(dec, n=2000)
    #rm(sce)
    #rm(dec)
    #VariableFeatures(seurat.obj4)<-top.hvgs2

    seurat.obj4<-FindVariableFeaturesMix(seurat.obj4,method.names = "scran",nfeatures = nfeatures)
    seurat.obj4<-ScaleData(seurat.obj4,verbose = F)
    seurat.obj4<-RunPCA(seurat.obj4,npcs=30,verbose=F)

    seurat.obj.pca[[4]]<-seurat.obj4@reductions$pca@cell.embeddings
    var.seurat.obj[[4]]<-VariableFeatures(seurat.obj4)

    rm(seurat.obj4)



    #################################################
    #5.PFlog1pPF selection: var(PFlog1pPF)~mean(PFlog1pPF)
    message("Method5: mv_PFlogPF")

    seurat.obj5<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj5<-NormalizeData(seurat.obj5,verbose = F)

    #library(scran)
    #sce <- SingleCellExperiment(list(counts=rna_mat_PFlog1pPF))
    #sce@assays@data$logcounts<-sce@assays@data$counts
    #dec <- modelGeneVar(sce)
    #top.hvgs <- getTopHVGs(dec, n=2000)
    #rm(sce)
    #VariableFeatures(seurat.obj5)<-top.hvgs

    seurat.obj5<-FindVariableFeaturesMix(seurat.obj5,method.names = "mv_PFlogPF",nfeatures = nfeatures)
    seurat.obj5<-ScaleData(seurat.obj5,verbose = F)
    seurat.obj5<-RunPCA(seurat.obj5,npcs=30,verbose=F)

    seurat.obj.pca[[5]]<-seurat.obj5@reductions$pca@cell.embeddings
    var.seurat.obj[[5]]<-VariableFeatures(seurat.obj5)

    rm(seurat.obj5)


    #################################################
    #6.seurat selection: log(var(count))~log(mean(count))
    message("Method6: seuratv3")
    seurat.obj6<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj6<-NormalizeData(seurat.obj6,verbose = F)
    seurat.obj6<-FindVariableFeaturesMix(seurat.obj6,method.names = "seuratv3",nfeatures = nfeatures)
    #seurat.obj6<-FindVariableFeatures(seurat.obj6,nfeatures=nfeatures,verbose = F)
    seurat.obj6<-ScaleData(seurat.obj6,verbose = F)
    seurat.obj6<-RunPCA(seurat.obj6,npcs=30,verbose=F)

    seurat.obj.pca[[6]]<-seurat.obj6@reductions$pca@cell.embeddings
    var.seurat.obj[[6]]<-VariableFeatures(seurat.obj6)
    rm(seurat.obj6)

    #################################################
    #7.normalized selection: log(var(normalized))~log(mean(normalized))
    message("Method7: logmv_nc")
    seurat.obj7<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj7<-NormalizeData(seurat.obj7,verbose = F)

    #seurat.obj7@assays$RNA@counts<-seurat.obj7@assays$RNA@data
    #seurat.obj7@assays$RNA@counts@x<-exp(seurat.obj7@assays$RNA@counts@x)-1
    #seurat.obj7<-FindVariableFeatures(seurat.obj7,nfeatures=nfeatures,verbose = F)

    seurat.obj7<-FindVariableFeaturesMix(seurat.obj7,method.names = "logmv_nc",nfeatures = nfeatures)
    seurat.obj7<-ScaleData(seurat.obj7,verbose = F)
    seurat.obj7<-RunPCA(seurat.obj7,npcs=30,verbose=F)

    seurat.obj.pca[[7]]<-seurat.obj7@reductions$pca@cell.embeddings
    var.seurat.obj[[7]]<-VariableFeatures(seurat.obj7)
    rm(seurat.obj7)


    #################################################
    #8.log-normalized selection: log(var(log1p(normalized)))~log(mean(log1p(normalized)))
    message("Method8: logmv_lognc")

    seurat.obj8<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj8<-NormalizeData(seurat.obj8,verbose = F)

    #seurat.obj8@assays$RNA@counts<-seurat.obj8@assays$RNA@data
    #seurat.obj8<-FindVariableFeatures(seurat.obj8,nfeatures=nfeatures,verbose = F)

    seurat.obj8<-FindVariableFeaturesMix(seurat.obj8,method.names = "logmv_lognc",nfeatures = nfeatures)
    seurat.obj8<-ScaleData(seurat.obj8,verbose = F)
    seurat.obj8<-RunPCA(seurat.obj8,npcs=30,verbose=F)

    seurat.obj.pca[[8]]<-seurat.obj8@reductions$pca@cell.embeddings
    var.seurat.obj[[8]]<-VariableFeatures(seurat.obj8)
    rm(seurat.obj8)


    #################################################
    #9.log-PFlog1pPF selection : log(var(PFlog1pPF))~log(mean(PFlog1pPF))
    message("Method9: logmv_PFlogPF")

    seurat.obj9<-CreateSeuratObject(rna_mat,verbose = F)
    #tmp<-seurat.obj9@assays$RNA@counts
    #seurat.obj9@assays$RNA@counts<-rna_mat_PFlog1pPF
    #seurat.obj9<-FindVariableFeatures(seurat.obj9,nfeatures=nfeatures,verbose = F)
    #seurat.obj9@assays$RNA@counts<-tmp

    seurat.obj9<-NormalizeData(seurat.obj9,verbose = F)
    seurat.obj9<-FindVariableFeaturesMix(seurat.obj9,method.names = "logmv_PFlogPF",nfeatures = nfeatures)
    seurat.obj9<-ScaleData(seurat.obj9,verbose = F)
    seurat.obj9<-RunPCA(seurat.obj9,npcs=30,verbose=F)

    seurat.obj.pca[[9]]<-seurat.obj9@reductions$pca@cell.embeddings
    var.seurat.obj[[9]]<-VariableFeatures(seurat.obj9)
    rm(seurat.obj9)


    #################################################
    #10.scran selection: possion: var(log1p(normalized))~ mean(log1p(normalized))
    message("Method10: scran_pos")
    ## scran

    seurat.obj10<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj10<-NormalizeData(seurat.obj10,verbose = F)

    #library(scran)
    #sce <- SingleCellExperiment(list(counts=rna_mat))
    #sce <- logNormCounts(sce)
    #dec<- modelGeneVarByPoisson(sce)
    #top.hvgs3 <- getTopHVGs(dec, n=2000)
    #rm(sce)
    #rm(dec)
    #VariableFeatures(seurat.obj10)<-top.hvgs3

    seurat.obj10<-FindVariableFeaturesMix(seurat.obj10,method.names = "scran_pos",nfeatures = nfeatures)
    seurat.obj10<-ScaleData(seurat.obj10,verbose = F)
    seurat.obj10<-RunPCA(seurat.obj10,npcs=30,verbose=F)

    seurat.obj.pca[[10]]<-seurat.obj10@reductions$pca@cell.embeddings
    var.seurat.obj[[10]]<-VariableFeatures(seurat.obj10)

    rm(seurat.obj10)



    #################################################
    #11.dispersion
    message("Method11: dispersion count")
    seurat.obj11<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj11<-NormalizeData(seurat.obj11,verbose = F)
    seurat.obj11<-FindVariableFeaturesMix(seurat.obj11,method.names = "disp_ct",nfeatures = nfeatures,verbose = F)
    seurat.obj11<-ScaleData(seurat.obj11,verbose = F)
    seurat.obj11<-RunPCA(seurat.obj11,npcs=30,verbose=F)

    seurat.obj.pca[[11]]<-seurat.obj11@reductions$pca@cell.embeddings
    var.seurat.obj[[11]]<-VariableFeatures(seurat.obj11)

    rm(seurat.obj11)

    #################################################
    #12.Seurat v1
    message("Method12: Seurat v1, disp_nc")
    seurat.obj12<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj12<-NormalizeData(seurat.obj12,verbose = F)
    seurat.obj12<-FindVariableFeatures(seurat.obj12,nfeatures=nfeatures,verbose = F,selection.method = "disp")
    seurat.obj12<-ScaleData(seurat.obj12,verbose = F)
    seurat.obj12<-RunPCA(seurat.obj12,npcs=30,verbose=F)

    seurat.obj.pca[[12]]<-seurat.obj12@reductions$pca@cell.embeddings
    var.seurat.obj[[12]]<-VariableFeatures(seurat.obj12)

    rm(seurat.obj12)

    #################################################
    #13.seurat
    message("Method13: seuratv2")
    ## seurat mvp
    seurat.obj13<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj13<-NormalizeData(seurat.obj13,verbose = F)
    seurat.obj13<-FindVariableFeatures(seurat.obj13,nfeatures=nfeatures,verbose = F,selection.method = "mvp")
    seurat.obj13<-ScaleData(seurat.obj13,verbose = F)
    seurat.obj13<-RunPCA(seurat.obj13,npcs=30,verbose=F)

    seurat.obj.pca[[13]]<-seurat.obj13@reductions$pca@cell.embeddings
    var.seurat.obj[[13]]<-VariableFeatures(seurat.obj13)

    rm(seurat.obj13)

    #################################################
    #14.dispersion lognc
    message("Method14: dispersion lognc")
    seurat.obj14<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj14<-NormalizeData(seurat.obj14,verbose = F)
    seurat.obj14<-FindVariableFeaturesMix(seurat.obj14,method.names = "disp_lognc",nfeatures = nfeatures,verbose = F)
    seurat.obj14<-ScaleData(seurat.obj14,verbose = F)
    seurat.obj14<-RunPCA(seurat.obj14,npcs=30,verbose=F)

    seurat.obj.pca[[14]]<-seurat.obj14@reductions$pca@cell.embeddings
    var.seurat.obj[[14]]<-VariableFeatures(seurat.obj14)

    rm(seurat.obj14)

    #################################################
    #15.dispersion PFlogPF
    message("Method15: dispersion PFlogPF")
    seurat.obj15<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj15<-NormalizeData(seurat.obj15,verbose = F)
    seurat.obj15<-FindVariableFeaturesMix(seurat.obj15,method.names = "disp_PFlogPF",nfeatures = nfeatures,verbose = F)
    seurat.obj15<-ScaleData(seurat.obj15,verbose = F)
    seurat.obj15<-RunPCA(seurat.obj15,npcs=30,verbose=F)

    seurat.obj.pca[[15]]<-seurat.obj15@reductions$pca@cell.embeddings
    var.seurat.obj[[15]]<-VariableFeatures(seurat.obj15)

    rm(seurat.obj15)


    #################################################
    #16.cell ranger selection (dispersion based, scanpy):
    message("Method16: cell ranger")
    sc = import("scanpy")
    seurat.obj16<-CreateSeuratObject(rna_mat,verbose = F)
    scanpy_rna=sc$AnnData( t((seurat.obj16@assays$RNA@counts)))
    scanpy_rna$obs_names = as.character(1:ncol(rna_mat))
    scanpy_rna$var_names = rownames(rna_mat)
    #sc$pp$filter_cells(scanpy_rna, min_genes=200)
    sc$pp$filter_genes(scanpy_rna, min_cells=1)
    sc$pp$normalize_total(scanpy_rna, target_sum=1e4)
    sc$pp$log1p(scanpy_rna)
    sc$pp$highly_variable_genes(scanpy_rna,n_top_genes=as.integer(2000),flavor='cell_ranger')
    #[‘seurat’, ‘cell_ranger’, ‘seurat_v3’]
    top_scanpy_cell_ranger<-scanpy_rna$var_names[scanpy_rna$var$highly_variable]
    top_scanpy_cell_ranger<-sapply(1:length(top_scanpy_cell_ranger), function(i){
      top_scanpy_cell_ranger[i-1]
    })

    seurat.obj16<-NormalizeData(seurat.obj16,verbose = F)
    VariableFeatures(seurat.obj16)<-top_scanpy_cell_ranger
    seurat.obj16<-ScaleData(seurat.obj16,verbose = F)
    seurat.obj16<-RunPCA(seurat.obj16,npcs=30,verbose=F)

    seurat.obj.pca[[16]]<-seurat.obj16@reductions$pca@cell.embeddings
    var.seurat.obj[[16]]<-VariableFeatures(seurat.obj16)
    rm(seurat.obj16)
    rm(scanpy_rna)


    #################################################
    #17.highly expressed genes
    message("Method17: mean_max_ct")
    seurat.obj17<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj17<-NormalizeData(seurat.obj17,verbose = F)
    seurat.obj17<-FindVariableFeaturesMix(seurat.obj17,method.names = "mean_max_ct",nfeatures = nfeatures,verbose = F)
    seurat.obj17<-ScaleData(seurat.obj17,verbose = F)
    seurat.obj17<-RunPCA(seurat.obj17,npcs=30,verbose=F)

    seurat.obj.pca[[17]]<-seurat.obj17@reductions$pca@cell.embeddings
    var.seurat.obj[[17]]<-VariableFeatures(seurat.obj17)

    rm(seurat.obj17)

    #################################################
    #18.highly expressed genes normalized count
    message("Method18: mean_max_nc")
    seurat.obj18<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj18<-NormalizeData(seurat.obj18,verbose = F)
    seurat.obj18<-FindVariableFeaturesMix(seurat.obj18,method.names = "mean_max_nc",nfeatures = nfeatures,verbose = F)
    seurat.obj18<-ScaleData(seurat.obj18,verbose = F)
    seurat.obj18<-RunPCA(seurat.obj18,npcs=30,verbose=F)

    seurat.obj.pca[[18]]<-seurat.obj18@reductions$pca@cell.embeddings
    var.seurat.obj[[18]]<-VariableFeatures(seurat.obj18)

    rm(seurat.obj18)


    #################################################
    #19.highly expressed genes log normalized count
    message("Method19: mean_max_lognc")
    seurat.obj19<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj19<-NormalizeData(seurat.obj19,verbose = F)
    seurat.obj19<-FindVariableFeaturesMix(seurat.obj19,method.names = "mean_max_lognc",nfeatures = nfeatures,verbose = F)
    seurat.obj19<-ScaleData(seurat.obj19,verbose = F)
    seurat.obj19<-RunPCA(seurat.obj19,npcs=30,verbose=F)

    seurat.obj.pca[[19]]<-seurat.obj19@reductions$pca@cell.embeddings
    var.seurat.obj[[19]]<-VariableFeatures(seurat.obj19)

    rm(seurat.obj19)

    #################################################
    #20.highly expressed genes PFlogPF
    message("Method20: mean_max_PFlogPF")
    seurat.obj20<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj20<-NormalizeData(seurat.obj20,verbose = F)
    seurat.obj20<-FindVariableFeaturesMix(seurat.obj20,method.names = "mean_max_PFlogPF",nfeatures = nfeatures,verbose = F)
    seurat.obj20<-ScaleData(seurat.obj20,verbose = F)
    seurat.obj20<-RunPCA(seurat.obj20,npcs=30,verbose=F)

    seurat.obj.pca[[20]]<-seurat.obj20@reductions$pca@cell.embeddings
    var.seurat.obj[[20]]<-VariableFeatures(seurat.obj20)

    rm(seurat.obj20)

    #################################################
    #21.Seurat SCT
    message("Method21: SCT")
    seurat.obj21<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj21<-SCTransform(seurat.obj21,variable.features.n=2000,verbose = F)
    seurat.obj21<-RunPCA(seurat.obj21,npcs=30,verbose=F)

    seurat.obj.pca[[21]]<-seurat.obj21@reductions$pca@cell.embeddings
    var.seurat.obj[[21]]<-VariableFeatures(seurat.obj21)

    list("seurat.obj.pca"=seurat.obj.pca,
         "var.seurat.obj"=var.seurat.obj)

}


# saved as paste0(save_path,dataset_name,"_pcalist_mixture.rds")
mixture_hvg_pca<-function(rna_mat,
                          nfeatures = 2000){

    seurat.obj.pca<-list()
    var.seurat.obj<-list()

    #seurat.obj0<-CreateSeuratObject(rna_mat,verbose = F)
    ##rna_mat<-seurat.obj0@assays$RNA@counts
    #rm(seurat.obj0)
    #rna_mat_PFlog1pPF<-t(t(rna_mat)/colSums(rna_mat))*mean(colSums(rna_mat))
    #rna_mat_PFlog1pPF<-log1p(rna_mat_PFlog1pPF)
    #rna_mat_PFlog1pPF<-t(t(rna_mat_PFlog1pPF)/colSums(rna_mat_PFlog1pPF))*mean(colSums(rna_mat_PFlog1pPF))

    rna_mat_PFlog1pPF<-NormalizeData(rna_mat,scale.factor=mean(colSums(rna_mat)))
    rna_mat_PFlog1pPF<-NormalizeData(rna_mat_PFlog1pPF,scale.factor=mean(colSums(rna_mat_PFlog1pPF)),normalization.method = "RC")


    sce <- SingleCellExperiment(list(counts=rna_mat))
    sce <- logNormCounts(sce)
    dec <- modelGeneVar(sce)
    dec.var <- dec@listData$bio
    dec.keep <- !is.na(dec.var) & dec.var > 0
    v_scran<-dec.var
    v_scran[!dec.keep]<-0
    rm(sce)
    rm(dec)

    seurat.obj1<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj1<-NormalizeData(seurat.obj1,verbose = F)
    seurat.obj1<-FindVariableFeatures(seurat.obj1,nfeatures=nfeatures,verbose = F)
    v_seuratv3<-seurat.obj1@assays$RNA@meta.features$vst.variance.standardized

    v_maxct<-rowMeans(seurat.obj1@assays$RNA@data)
    rm(seurat.obj1)
    seurat.obj11<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj11<-NormalizeData(seurat.obj11,verbose = F)
    seurat.obj11<-FindVariableFeatures(seurat.obj11,nfeatures=nfeatures,verbose = F,selection.method = "disp")
    v_disp<-seurat.obj11@assays$RNA@meta.features$mvp.dispersion
    v_disp[is.na(v_disp)]<-0
    v_disp[v_disp<0]<-0
    rm(seurat.obj11)

    seurat.obj10<-CreateSeuratObject(rna_mat,verbose = F)
    seurat.obj10<-NormalizeData(seurat.obj10,verbose = F)
    rna_mat_norm<-seurat.obj10@assays$RNA@data
    rna_mat_norm@x<-exp(rna_mat_norm@x)-1
    #library(scran)
    sce <- SingleCellExperiment(list(counts=rna_mat_norm))
    sce@assays@data$logcounts<-sce@assays@data$counts
    dec <- modelGeneVar(sce)
    top.hvgs <- getTopHVGs(dec, n=2000)
    dec.var <- dec@listData$bio
    dec.keep <- !is.na(dec.var) & dec.var > 0
    v_scran_nc<-dec.var
    v_scran_nc[!dec.keep]<-0
    rm(sce)
    rm(dec)

    sce <- SingleCellExperiment(list(counts=rna_mat))
    sce <- logNormCounts(sce)
    dec<- modelGeneVarByPoisson(sce)
    dec.var <- dec@listData$bio
    dec.keep <- !is.na(dec.var) & dec.var > 0
    v_scran_pos<-dec.var
    v_scran_pos[!dec.keep]<-0
    rm(sce)
    rm(dec)

    mixture_rank<-function(input_lst){
        input_lst_order<-list()
        for(i in 1:length(input_lst)){
            input_lst_order[[i]]<-order(order(input_lst[[i]],decreasing = T))
        }

        apply(matrix(unlist(input_lst_order),
                     ncol=length(input_lst_order),byrow = FALSE),1,FUN = min)
    }
    mixture_mat<-cbind(v_scran_nc,v_scran,v_scran_pos,v_disp,v_maxct,v_seuratv3)
    mixture_index_list<-list(c(1,2),c(1,3),c(1,4),c(1,5),c(2,3),c(2,4),c(2,5),c(3,4),c(3,5),c(4,5),
                          c(1,2,3),c(1,2,4),c(1,2,5),c(1,3,4),c(1,3,5),c(1,4,5),
                          c(2,3,4),c(2,3,5),c(2,4,5),c(3,4,5),
                          c(1,2,3,4),c(1,2,3,5),c(1,2,4,5),c(1,3,4,5),c(2,3,4,5),
                          c(1,2,3,4,5),c(1,2,3,4,5,6))

    hvg<-list()

    for(kk in 1:length(mixture_index_list)){
        cur_mixture_index<-mixture_index_list[[kk]]
        cur_mixture_list<-list()
        cur_mixture_list_count<-1
        for( i in cur_mixture_index){
            cur_mixture_list[[cur_mixture_list_count]]<-mixture_mat[,i]
            cur_mixture_list_count<-cur_mixture_list_count+1
        }
        hvg[[kk]]<-order(mixture_rank(cur_mixture_list))[1:nfeatures]
    }

    print("Method_Mixture")
    rna_mat_norm<-NormalizeData(rna_mat,verbose = FALSE)

    for( i in 1:length(hvg)){
        message(mixture_index_list[i])
        cur_index<-i
        rna_mat_norm_hvg<-rna_mat_norm[hvg[[i]],]
        rna_mat_scale<-ScaleData(rna_mat_norm_hvg,verbose = F)
        suppressWarnings(rna_mat_pca <- RunPCA(
            object = rna_mat_scale,
            features = rownames(rna_mat_scale),
            npcs = 30,
            verbose = FALSE)@cell.embeddings)
        seurat.obj.pca[[cur_index]]<-rna_mat_pca
        var.seurat.obj[[cur_index]]<-rownames(rna_mat)[hvg[[i]]]
    }
    newList<-list("seurat.obj.pca"=seurat.obj.pca,
                  "var.seurat.obj"=var.seurat.obj)
    print("Finish Methods Mixture!")
    return(newList)
}


