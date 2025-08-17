myendr = function (dat, config = endr_defs, ...) 
{
  config <- EDGE:::check_config(config, ...)
  # dat <- input_check(dat)
  if (config$n_dm > ncol(dat)) {
    endr_error("'n_dm' must be less than the number of columns of the data")
  }
  if (is.na(config$a) | is.na(config$b)) {
    config[c("a", "b")] <- EDGE:::ab_params(config$spread, config$min_dist)
  }
  n_wl <- config$n_wl
  n_dm <- config$n_dm
  n_comps <- config$n_comps
  seed <- config$seed
  nk <- config$n_neigs
  H <- config$H
  if (dim(dat)[1] > H) {
    H <- 1017881
  }
  amat <- affinity(dat, n_wl = n_wl, n_dm = n_dm, nk = nk, 
                   H = H, seed = seed)
  amat[, 3] <- amat[, 3]/(2 * n_wl)
  if (config$opt) {
    amat[amat[, 1] == amat[, 2], 3] <- 0
    colnames(amat) <- c("from", "to", "value")
    coor_knn <- EDGE:::make_coor(amat, rownames(dat), nrow(dat))
    amat_spec <- EDGE:::spectral_eigen(coor_knn, n_comps)
    opt_embedding <- EDGE:::endr_embed(coor_knn, amat_spec, config)
  }
  else {
    colnames(amat) <- c("from", "to", "value")
    coor_knn <- EDGE:::make_coor(amat, rownames(dat), nrow(dat))
    amat_spec <- EDGE:::spectral_eigen(coor_knn, n_comps)
    opt_embedding <- amat_spec
  }
  out_embedding <- EDGE:::center_embed(opt_embedding)
  class(out_embedding) <- "endr"
  out_embedding
}