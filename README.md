The files are for benchmarking highly variable gene(HVG) selection methods.

**Benchmark_HVG_1_Methods.R** contains the processing steps for all the baseline methods and methods mixture. Refer to https://github.com/RuzhangZhao/mixhvg for the R package.

**Benchmark_HVG_2_Evaluation_Criteria.R**: contains the all the evaluation criteria for cell sorting, CITEseq, and MultiomeATAC, separately.

**Benchmark_HVG_3_RunExample.R**: contains examples for processing HVG selection methods and running evaluation.

```R
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
```



```R
## Example: CITEseq
## bmcite

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

```



Download pbmc file from [pbmc_rna.rds](https://github.com/RuzhangZhao/pbmc3k/blob/main/pbmc3k_rna.rds) and [pbmc_lsi.rds](https://github.com/RuzhangZhao/pbmc3k/blob/main/pbmc3k_lsi.rds)

```R
## Example: MultiomeATAC
## pbmc3k_multi
## Assume we have processed pbmc3k_multi
## and get pbmc3k_multi_rna_mat and pbmc3k_multi_lsi
# baseline methods
pbmc3k_multi_rna_mat<-readRDS("pbmc_rna.rds")
pbmc3k_multi_lsi<-readRDS("pbmc_lsi.rds")[1:4,]
pcalist<-hvg_pca(pbmc3k_multi_rna_mat)
pcalist<-pcalist$seurat.obj.pca
# mixture methods
mixpcalist<-mixture_hvg_pca(pbmc3k_multi_rna_mat)
mixpcalist<-mixpcalist$seurat.obj.pca
multi_pca_eval<-evaluate_hvg_continuous(pcalist=pcalist,pro=pbmc3k_multi_lsi,input="MultiomeATAC")
multi_mixpca_eval<-evaluate_hvg_continuous(pcalist=mixpcalist,pro=pbmc3k_multi_lsi,input="MultiomeATAC")

```

