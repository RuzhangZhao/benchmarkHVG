evaluate_f1_from_embedding <- function(embedding,
                                       cell_label,
                                       algorithm=1) {
  unique_names <- unique(cell_label)
  unique_names <- unique_names[unique_names != "Others"]
  snn_ <- FindNeighbors(
    object = embedding,
    nn.method = "rann",
    verbose = F
  )$snn
  resolution_range <- seq(0.1, 2, 0.1)
  res <- sapply(resolution_range, function(cur_resolution) {
    cluster_label <- FindClusters(snn_,
                                  resolution = cur_resolution,
                                  algorithm = algorithm,
                                  verbose = F
    )[[1]]
    x <- sapply(unique(cluster_label), function(i) {
      cluster_cell_label <- rep("Others", length(cell_label))
      cluster_cell_label[cluster_label == i] <- unique_names
      compute_f1(cell_label, cluster_cell_label, unique_names)
    })
    max(x)
  })
  res
}

evaluate_f1_from_label <- function(cluster_label,
                                   cell_label) {
  unique_names <- unique(cell_label)
  unique_names <- unique_names[unique_names != "Others"]
  
  x <- sapply(unique(cluster_label), function(i) {
    cluster_cell_label <- rep("Others", length(cell_label))
    cluster_cell_label[cluster_label == i] <- unique_names
    compute_f1(cell_label, cluster_cell_label, unique_names)
  })
  max(x)
  
}


compute_f1 <- function(true_labels, predictions, positive_label) {
  conf_matrix <- table(Predicted = predictions, Actual = true_labels)
  TP <- conf_matrix[positive_label, positive_label]
  FP <- sum(conf_matrix[positive_label, ]) - TP
  FN <- sum(conf_matrix[, positive_label]) - TP
  
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  f1_score <- 2 * (precision * recall) / max(precision + recall, 1e-10)
  
  return(f1_score)
}