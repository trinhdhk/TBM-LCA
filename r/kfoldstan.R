#functions slightly modified from: https://github.com/stan-dev/stancon_talks/blob/master/2017/Contributed-Talks/07_nicenboim/kfold.Rmd

#function to parrallelize all computations
#need at least two chains !!!
stan_kfold <- function(file, model, list_of_datas, chains, cores,...){
  library(pbmcapply)
  badRhat <- 1.1 # don't know why we need this?
  n_fold <- length(list_of_datas)
  if (missing(model)) model <- stan_model(file=file)
  # First parallelize all chains:
  sflist <- 
    pbmclapply(1:(n_fold*chains), mc.cores = cores, 
               function(i){
                 # Fold number:
                 k <- ceiling(i / chains)
                 s <- sampling(model, data = list_of_datas[[k]], 
                               chains = 1, chain_id = i,...)
                 return(s)
               })
  
  # Then merge the K * chains to create K stanfits:
  stanfit <- list()
  for(k in 1:n_fold){
    inchains <- (chains*k - (chains - 1)):(chains*k)
    #  Merge `chains` of each fold
    stanfit[[k]] <- sflist2stanfit(sflist[inchains])
  }  
  return(stanfit) 
}

#extract log-likelihoods of held-out data
extract_log_lik_K <- function(list_of_stanfits, list_of_holdout, ...){
  require(loo)
  K <- length(list_of_stanfits)
  list_of_log_liks <- plyr::llply(1:K, function(k){
    extract_log_lik(list_of_stanfits[[k]],...)
  })
  # `log_lik_heldout` will include the loglike of all the held out data of all the folds.
  # We define `log_lik_heldout` as a (samples x N_obs) matrix
  # (similar to each log_lik matrix)
  log_lik_heldout <- list_of_log_liks[[1]] * NA
  for(k in 1:K){
    log_lik <- list_of_log_liks[[k]]
    samples <- dim(log_lik)[1] 
    N_obs <- dim(log_lik)[2]
    # This is a matrix with the same size as log_lik_heldout
    # with 1 if the data was held out in the fold k
    heldout <- matrix(rep(list_of_holdout[[k]], each = samples), nrow = samples)
    # Sanity check that the previous log_lik is not being overwritten:
    if(any(!is.na(log_lik_heldout[heldout==1]))){
      warning("Heldout log_lik has been overwritten!!!!")
    }
    # We save here the log_lik of the fold k in the matrix:
    log_lik_heldout[heldout==1] <- log_lik[heldout==1]
  }
  return(log_lik_heldout)
}

#compute ELPD
kfold <- function(log_lik_heldout)  {
  library(matrixStats)
  logColMeansExp <- function(x) {
    # should be more stable than log(colMeans(exp(x)))
    S <- nrow(x)
    colLogSumExps(x) - log(S)
  }
  # See equation (20) of @VehtariEtAl2016
  pointwise <-  matrix(logColMeansExp(log_lik_heldout), ncol= 1)
  colnames(pointwise) <- "elpd"
  # See equation (21) of @VehtariEtAl2016
  elpd_kfold <- sum(pointwise)
  se_elpd_kfold <-  sqrt(ncol(log_lik_heldout) * var(pointwise))
  out <- list(
    pointwise = pointwise,
    elpd_kfold = elpd_kfold,
    se_elpd_kfold = se_elpd_kfold)
  #structure(out, class = "loo")
  return(out)
}