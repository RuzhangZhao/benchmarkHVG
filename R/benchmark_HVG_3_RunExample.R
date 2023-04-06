# Author: Ruzhang Zhao
## Example: CITEseq
## bmcite
if(0){
library(mixhvg)
library(SeuratData)
bm <- LoadData(ds = "bmcite")
expr<-bm@assays$RNA@counts
pro<-bm@assays$ADT@data
# baseline methods
pcalist<-hvg_pca(expr)
pcalist<-pcalist$seurat.obj.pca
# mixture methods
mixpcalist<-mixture_hvg_pca(expr)
mixpcalist<-mixpcalist$seurat.obj.pca
cite_pca_eval<-evaluate_hvg_continuous(pcalist=pcalist,pro=pro,input="CITEseq")
cite_mixpca_eval<-evaluate_hvg_continuous(pcalist=mipcalist,pro=pro,input="CITEseq")

## Example: MultiomeATAC
## pbmc3k_multi
## Assume we have processed pbmc3k_multi
## and get pbmc3k_multi_rna_mat and pbmc3k_multi_lsi
# baseline methods
pcalist<-hvg_pca(pbmc3k_multi_rna_mat)
pcalist<-pcalist$seurat.obj.pca
# mixture methods
mixpcalist<-mixture_hvg_pca(pbmc3k_multi_rna_mat)
mixpcalist<-mixpcalist$seurat.obj.pca
multi_pca_eval<-evaluate_hvg_continuous(pcalist=pcalist,pro=pbmc3k_multi_lsi,input="MultiomeATAC")
multi_mixpca_eval<-evaluate_hvg_continuous(pcalist=mixpcalist,pro=pbmc3k_multi_lsi,input="MultiomeATAC")

## Example: Cell Sorting
## duo4_pbmc
library(DuoClustering2018)
sce_names<-c("Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
sn<-sce_names[1]
sce<-do.call(paste0("sce_full_",sn),list())
expr<-sce@assays$data@listData$counts
cell_label<-sce$phenoid
# baseline methods
pcalist<-hvg_pca(expr)
pcalist<-pcalist$seurat.obj.pca
# mixture methods
mixpcalist<-mixture_hvg_pca(expr)
mixpcalist<-mixpcalist$seurat.obj.pca
sort_pca_eval<-evaluate_hvg_discrete(pcalist=pcalist,label=cell_label)
sort_mixpca_eval<-evaluate_hvg_discrete(pcalist=mixpcalist,label=cell_label)
}
