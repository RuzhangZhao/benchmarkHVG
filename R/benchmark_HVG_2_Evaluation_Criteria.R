# Author: Ruzhang Zhao
# Evaluate with cell sorting

#' Criterion 1: ARI
#' @details
#'
#' @param
#'
#'
#' @return
#'
#' @import Seurat
#' @importFrom mclust adjustedRandIndex
#' @export
#'
ARILouvain<-function(
        embedding,
        cell_label,
        maxit = 25
){
    N_label<-length(unique(cell_label))
    snn_<- FindNeighbors(object = embedding,
                         nn.method = "rann",
                         verbose = F)$snn
    cluster_label <- FindClusters(snn_,
                                  resolution = 0,
                                  verbose = F)[[1]]
    N_cluster0<-length(unique(cluster_label))
    if(N_cluster0>=N_label){
        print("Resolution 0 still larger than true")
        return(c(N_cluster0,
                 adjustedRandIndex(cluster_label,cell_label)))
    }

    snn_<- FindNeighbors(object = embedding,
                         nn.method = "rann",
                         verbose = F)$snn
    cluster_label <- FindClusters(snn_,
                                  resolution = 1,
                                  verbose = F)[[1]]
    N_cluster1<-length(unique(cluster_label))
    if(N_cluster1 <= N_label ){
        return(c(N_cluster1,
                 adjustedRandIndex(cluster_label,cell_label)))
    }else if (N_cluster1 > N_label){
        a<-0
        b<-1
    }
    keepsign<-TRUE
    iter_<-0
    while (keepsign) {
        c<-(a+b)/2
        snn_<- FindNeighbors(object = embedding,
                             nn.method = "rann",
                             verbose = F)$snn
        cluster_label <- FindClusters(snn_,
                                      resolution = c,
                                      verbose = F)[[1]]
        N_clusterc<-length(unique(cluster_label))
        if( N_clusterc == N_label){
            keepsign<-FALSE
            return(c(N_clusterc,
                     adjustedRandIndex(cluster_label,cell_label)))
        }else if (N_clusterc < N_label){
            a<-c
        }else{
            b<-c
        }
        iter_<-iter_+1
        if(iter_ > maxit){
            keepsign<-FALSE
        }
    }
    return(c(N_clusterc,
             adjustedRandIndex(cluster_label,cell_label)))
}

#' Criterion 2: Within VS Between Cell Type Variance Ratio
#' @details
#'
#' @param embedding
#' @param cell_label
#'
#'
#' @return
#'
#' @importFrom Rfast colVars
#' @export
#'
within_between_var_ratio<-function(embedding,cell_label){
    cell_label <- as.character(cell_label)
    # sorted unique cluster label
    label_index<-sort(unique(as.character(cell_label)))
    if(sum(label_index == "") >0){
        label_index<-label_index[-which(label_index == "")]
    }
    # number of labels
    N_label<-length(label_index)
    within<-sapply(1:N_label, function(i){
        index_i<-which(cell_label == label_index[i])
        (length(index_i)-1)*colVars(embedding[index_i,])
    })
    all<-(sum(nrow(embedding)-1)*colVars(embedding))
    sum(rowSums(within))/(sum(all)- sum(rowSums(within)))
}

#' Criterion 3: Nearest Neighbor Accuracy
#' @details
#'
#' @param embedding
#' @param cell_label
#' @param k
#' @param cutoff
#'
#'
#' @return
#'
#' @importFrom caret createDataPartition
#' @importFrom FNN knn
#' @export
#'
knn_accuracy<-function(embedding,cell_label,k = 3,cutoff = 3000){
    if (length(cell_label) > cutoff){
        sample_index<-createDataPartition(cell_label,p = cutoff/length(cell_label))$Resample1
    }else{
        sample_index<-1:length(cell_label)
    }

    N_cell<-nrow(embedding)
    pre_label<-sapply(sample_index, function(cell){
        train_label<-cell_label[-cell]
        a<-FNN::knn(train = embedding[-cell,],
                    test = embedding[cell,],
                    cl = train_label,k = k,prob = T)
        a[1]
    })
    sum(pre_label == cell_label[sample_index])/length(sample_index)
}

#' Criterion 4: Within VS Between Cell Type Distance Ratio
#' @details
#'
#' @param embedding
#' @param cell_label
#'
#'
#' @return
#'
#' @export
#'
within_between_dist_ratio<-function(embedding,cell_label){
    cell_label <- as.character(cell_label)
    # sorted unique cluster label
    label_index<-sort(unique(as.character(cell_label)))
    if(sum(label_index == "") >0){
        label_index<-label_index[-which(label_index == "")]
    }
    cell_label_unique<-unique(cell_label)
    ncell = length(cell_label)
    label_dist = matrix(1,ncell,ncell)
    for (k in seq_along(cell_label_unique)) {
        id_k = cell_label == cell_label_unique[k]
        label_dist[id_k,id_k] = 0
    }
    pc_dist<-as.matrix(dist(embedding))
    sum(pc_dist*(1-label_dist))/sum(pc_dist*label_dist)
}
#' Comprehensive Evaluation for cell sorting.
#' @details
#'
#' @param pcalist
#' @param label
#'
#' @importFrom caret createDataPartition
#' @return var_ratio
#' @return ari
#' @return 3nn
#' @return dist_ratio
#'
#' @export
#'
evaluate_hvg_discrete<-function(pcalist,label){
    Num_method<-length(pcalist)
    Nosample<-FALSE
    if(length(label)<10000){Nosample<-TRUE}

    if(!Nosample){
        set.seed(1)
        index_sample_pca<-createDataPartition(label,p = min(1,10000/length(label)))$Resample1
        set.seed(2)
        index_sample_pca1<-createDataPartition(label,p = min(1,10000/length(label)))$Resample1
    }

    #################################################
    # Within Cell Type Variance / Between Cell Type Variance Ratio
    message("var_ratio")
    variance_ratio<-rep(NA,Num_method)
    for(i in 1:Num_method){
        if(Nosample){
            variance_ratio[i]<-within_between_var_ratio(pcalist[[i]],label)
        }else{
            variance_ratio[i]<-
                (within_between_var_ratio(pcalist[[i]][index_sample_pca,],label[index_sample_pca])+
                     within_between_var_ratio(pcalist[[i]][index_sample_pca1,],label[index_sample_pca1]))/2
        }
    }

    #################################################
    # ARI with Louvain clustering  same number
    message("ari_louvain")
    ari_list<-rep(NA,Num_method)
    for(i in 1:Num_method){
        if(Nosample){
            ari_list[i]<-ARILouvain(pcalist[[i]],label)[2]
        }else{
            ari_list[i]<-(ARILouvain(pcalist[[i]][index_sample_pca,],label[index_sample_pca])[2]
                          +ARILouvain(pcalist[[i]][index_sample_pca1,],label[index_sample_pca1])[2])/2
        }
    }

    #################################################
    # 3NN accuracy
    message("3nn")
    nn_acc<-rep(NA,Num_method)

    for(i in 1:Num_method){
        if(Nosample){
            nn_acc[i]<-knn_accuracy(pcalist[[i]],label)
        }else{
            nn_acc[i]<-knn_accuracy(pcalist[[i]][index_sample_pca,],label[index_sample_pca])
        }
    }

    #################################################
    # Distance Ratio
    dist_ratio<-rep(NA,Num_method)

    for(i in 1:Num_method){
        if(Nosample){
            dist_ratio[i]<-within_between_dist_ratio(pcalist[[i]],label)
        }else{
            dist_ratio[i]<-within_between_dist_ratio(pcalist[[i]][index_sample_pca,],label[index_sample_pca])
        }
    }
    return(list(
            "var_ratio"=variance_ratio,
            "ari"=ari_list,
            "3nn"=nn_acc,
            "dist_ratio"=dist_ratio))
}

# Evaluate with CITEseq & MultiomeATAC
#' Criterion 1: Variance Ratio
#' @details
#'
#' @param embedding
#' @param pro
#' @param resolution
#'
#' @import Seurat
#' @importFrom Rfast rowVars
#' @return
#'
#' @export
#'
within_between_var_ratio_continuous<-function(
        embedding,
        pro,
        resolution = 0.2){

    snn_<- FindNeighbors(object = embedding,
                         nn.method = "rann",
                         verbose = F)$snn
    cluster_label <- FindClusters(snn_,
                                  resolution = resolution,
                                  verbose = F)[[1]]

    cluster_label <- as.numeric(as.character(cluster_label))
    # sorted unique cluster label
    label_index<-sort(as.numeric(
        unique(as.character(cluster_label))))
    # number of cluster label
    N_label<-length(label_index)
    #print(paste0("Current Label Number is ",N_label))
    within<-sapply(1:N_label, function(i){
        index_i<-which(cluster_label == label_index[i])
        (length(index_i)-1)*rowVars(pro[,index_i])
    })
    all_var<-(ncol(pro)-1)*rowVars(pro)
    var_r<-rowSums(within)/(all_var- rowSums(within))
    mean(var_r)
}

#' Criterion 2: Nearest Neighbor Mean Square Error
#' @details
#'
#' @param embedding
#' @param pro
#' @param k
#' @param cutoff
#'
#' @importFrom FNN knn
#' @return
#'
#' @export
#'
knn_regression<-function(embedding,pro,k=3,
                         cutoff = 2000){
    if (nrow(embedding) > cutoff){
        sample_index<-sample(1:nrow(embedding),size = cutoff)
    }else{
        sample_index<-1:nrow(embedding)
    }
    pred_val<-sapply(sample_index, function(cell){
        pro1<-pro[,-cell]
        ind<-FNN::knn(train=embedding[-cell,],
                    test= embedding[cell,],
                    cl = rep(1,nrow(embedding)-1),k=k)
        rowMeans(pro1[,attr(ind,"nn.index")[1,]])
    })
    sum((pred_val - pro[,sample_index])^2)/nrow(pro)
}

#' Criterion 3: Nearest Neighbor Distance Ratio
#' @details
#'
#' @param embedding
#' @param pro
#' @param k
#' @param cutoff
#'
#' @importFrom FNN knn
#' @importFrom pdist pdist
#' @return
#'
#' @export
#'
knn_ratio<-function(embedding,pro,k = 100,
                    cutoff = 2000){
    if (nrow(embedding) > cutoff){
        sample_index<-sample(1:nrow(embedding),size = cutoff)
    }else{
        sample_index<-1:nrow(embedding)
    }
    N_cell<-length(sample_index)
    k <- min(k,ceiling(N_cell/2))
    pro<-t(pro)
    dist_ratio<-sapply(sample_index, function(cell){
        ind<-FNN::knn(train = embedding[-cell,],
                    test = embedding[cell,],
                    cl = rep(1,nrow(embedding)-1),
                    k = k,prob = T)
        nn_index<-c(1:nrow(embedding))[-cell][attr(ind,"nn.index")[1,]]
        r<-mean(pdist(pro[cell,],pro[nn_index,])@dist)/
            mean(pdist(pro[cell,],pro[-c(nn_index,cell),])@dist)
        r
    })
    mean(dist_ratio)
}

# Criterion 4: Distance Correlation
# just cor function


#' Comprehensive Evaluation for CITEseq & MultiomeATAC
#' Select input from "MultiomeATAC" or "CITEseq
#' @details
#'
#' @param pcalist
#' @param pro
#' @param input
#'
#' @import Seurat
#' @return var_ratio
#' @return knn_ratio
#' @return 3nn
#' @return dist_cor
#'
#' @export
#'
evaluate_hvg_continuous<-function(pcalist,pro,
                               input="MultiomeATAC" ){
    if (input == "CITEseq"){
        scale_pro<-CreateSeuratObject(pro,verbose=FALSE)
        scale_pro <- NormalizeData(scale_pro, normalization.method = "CLR", margin = 2,verbose=F)
        scale_pro <- ScaleData(scale_pro,verbose=FALSE)
        pro <- scale_pro@assays$RNA@scale.data
    }
    Num_method<-length(pcalist)
    Nosample<-FALSE
    if(Num_method>10000){
        set.seed(1)
        index_sample_pca<-sample(1:ncol(pro),size = 10000)
        set.seed(2)
        index_sample_pca1<-sample(1:ncol(pro),size = 10000)
    }else{
        index_sample_pca<-1:ncol(pro)
        Nosample<-TRUE
    }


    #################################################
    # Within Between Cluster Variance Ratio
    message("var_ratio")
    variance_ratio<-rep(NA,Num_method)
    for(i in 1:Num_method){
        if(Nosample){
            variance_ratio[i]<-within_between_var_ratio_continuous(pcalist[[i]][index_sample_pca,],pro[,index_sample_pca])
        }else{
            variance_ratio[i]<-(within_between_var_ratio_continuous(pcalist[[i]][index_sample_pca,],pro[,index_sample_pca])+
                                    within_between_var_ratio_continuous(pcalist[[i]][index_sample_pca1,],pro[,index_sample_pca1]))/2
        }
    }

    #################################################
    # Within Nearest Neighbor VS Out of Nearest Neighbor Distance Ratio
    message("knn_ratio")
    knnratio<-rep(NA,Num_method)
    for(i in 1:Num_method){
        knnratio[i]<-knn_ratio(pcalist[[i]][index_sample_pca,],pro[,index_sample_pca])
    }


    #################################################
    # 3NN Regression MSE
    message("3nn")
    nn_mse<-rep(NA,Num_method)
    for(i in 1:Num_method){
        nn_mse[i]<-knn_regression(pcalist[[i]][index_sample_pca,],pro[,index_sample_pca])
    }

    #################################################
    # Distance Correlation
    dist_cor<-rep(NA,Num_method)
    pro_dist<-dist(t(pro[,index_sample_pca]))
    for(i in 1:Num_method){
        pc_dist<-dist(pcalist[[i]][index_sample_pca,])
        dist_cor[i]<-cor(c(pro_dist),c(pc_dist))
    }
    newList<-list(
        "var_ratio"=variance_ratio,
        "knn_ratio"=knnratio,
        "3nn"=nn_mse,
        "dist_cor"=dist_cor)
    return(newList)
}

# The resolution used for evaluation in all datasets

resolutionlist = list()
resolutionlist[["bmcite"]] = 0.1
resolutionlist[["pbmc_cite"]] = 0.2
resolutionlist[["cbmc8k_cite"]] = 0.1
resolutionlist[["FLiver_cite"]] = 0.05
resolutionlist[["FBM_cite"]] = 0.1
resolutionlist[["seurat_cite"]] = 0.08
resolutionlist[["sucovid_cite"]] = 0.3
resolutionlist[["homo_brain3k"]] = 0.15
resolutionlist[["lymphoma"]] = 0.05
resolutionlist[["mus_brain5k"]] = 0.1
resolutionlist[["pbmc3k_multi"]] = 0.1
resolutionlist[["pbmc10k_multi"]] = 0.1
