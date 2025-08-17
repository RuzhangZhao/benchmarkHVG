# Author: Ruzhang Zhao
## Example: CITEseq
## bmcite
library(mixhvg)
library(SeuratData)
bm <- LoadData(ds = "bmcite")
expr<-bm@assays$RNA@counts
pro<-bm@assays$ADT@data
# baseline methods
outputs<-rcti_pca(expr)

## Example: MultiomeATAC
## pbmc3k_multi
## Assume we have processed pbmc3k_multi
## and get pbmc3k_multi_rna_mat and pbmc3k_multi_lsi
# baseline methods
outputs<-rcti_pca(pbmc3k_multi_rna_mat)

## Example: Cell Sorting
## duo4_pbmc
library(DuoClustering2018)
sce_names<-c("Zhengmix4eq","Zhengmix4uneq","Zhengmix8eq")
sn<-sce_names[1]
sce<-do.call(paste0("sce_full_",sn),list())
expr<-sce@assays$data@listData$counts
cell_label<-sce$phenoid
# baseline methods
outputs<-rcti_pca(expr)

