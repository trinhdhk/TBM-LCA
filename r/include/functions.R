# Function to change y label for bayesplots.
# Trinh Dong, 2020-12
# mcmc_plot: a plot from bayesplot
# ..., labs: new labels, in spread form or vector form.
# top_down: TRUE, whether the labels will be from top down or bottom up.
change_ylabs <- function(mcmc_plot, ..., labs = character(), top_down = TRUE){
  labs <- c(..., labs)
  params <- unique(mcmc_plot$data$parameter)
  stopifnot(length(labs) == length(params))
  mcmc_plot + ggplot2::scale_y_discrete(limits = if(top_down) rev(params) else params,
                                        labels = if(top_down) rev(labs) else labs)
}

# functions to create repeated kfold
repeated_kfold <- function(K = 10, N_rep = 4, N_obs, data, seed = 123, cores = 10){
  holdout <- vector("list", N_rep * K)
  keptin <- vector("list", N_rep * K)
  for(n in seq_len(N_rep)){
    set.seed(seed + n - 1)
    hh <- loo::kfold_split_random(K = K, N = N_obs)
    foldkept <- matrix(1, nrow = N_obs, ncol = K)
    for(i in 1:N_obs) foldkept[i, hh[i]] <- 0
    foldhold  <- 1- foldkept
    keptin[((n-1)*K+1):(n*K)] <- split(foldkept,rep(1:ncol(foldkept),each=nrow(foldkept)))
    holdout[((n-1)*K+1):(n*K)]  <- split(foldhold,rep(1:ncol(foldhold),each=nrow(foldhold)))
  }
  inputs <- pbmcapply::pbmclapply(seq_len(K * N_rep),
                                  mc.cores=cores,
                                  function(k) modifyList(data, list(keptin=keptin[[k]])))
  list(keptin=keptin, holdout=holdout, inputs = inputs)
}


#vb_kfold <- 

#functions slightly modified from: https://github.com/stan-dev/stancon_talks/blob/master/2017/Contributed-Talks/07_nicenboim/kfold.Rmd

#function to parrallelize all computations
#need at least two chains !!!
stan_kfold <- function(file, sampler, list_of_datas, include_paths=NULL, sample_dir = NULL, backend = c("cmdstanr", "rstan"), chains, cores, seed, pars, merge = TRUE,...){
  # library(pbmcapply)
  backend <- match.arg(backend)
  badRhat <- 1.1 # don't know why we need this?
  n_fold <- length(list_of_datas)
  if (missing(sampler)){
    model <- if (backend == "cmdstanr")
      cmdstanr::cmdstan_model(stan_file=file, include_paths = include_paths)
    else rstan::stan_model(file=file)
  } else {
    model <- sampler
    backend2 <- switch(class(sampler), stanmodel = "rstan", CmdStanModel = "cmdstanr")
    if (backend != backend2){
      cat("backend set to ", backend2)
      backend <- backend2
    }
  }
  
  # First parallelize all chains:
  future::plan(future::multisession, workers=cores)
  wd <- getwd()
  progressr::with_progress({
    p <- progressr::progressor(steps = n_fold*chains)
    #sflist <- vector("list", length(n_fold*chains))
    sflist <- furrr::future_map(seq_len(n_fold*chains), function(i){
    # for (i in seq_len(n_fold*chains))
      # sflist[[i]] <- future::future({
      setwd(wd)
      k <- ceiling(i / chains)
      if (backend == "rstan"){
        if (length(sample_dir)){
          sample_file <- file.path(sample_dir, paste0("data_", k, "_chain_", chains-(k*chains-i), ".csv"))
          sf <- rstan::sampling(model, data = list_of_datas[[k]], chains=1, seed = seed, 
                                chain_id = i, sample_file = sample_file, pars = pars,...)
        } else 
          
          sf <- rstan::sampling(model, data = list_of_datas[[k]], chains=1, seed = seed, 
                                chain_id = i, pars = pars, ...)
      }
        
      else {
        m <- model$clone()
        sf <- m$sample(data = list_of_datas[[k]],
                      chains = 1, seed = seed + i, chain_ids = i,
                      output_dir = sample_dir,...)
      }
      p()
      sf
      # }, seed=T, gc=T)
    }, .options = furrr::furrr_options(seed=TRUE))
  })
  # sflist = lapply(sflist, future::value)
  
  # Then merge the K * chains to create K stanfits:
  if (!merge) return(sflist) 
  if (backend=="rstan"){
    # browser()
    stanfit <- list()
    for(k in 1:n_fold){
      inchains <- (chains*k - (chains - 1)):(chains*k)
      #  Merge `chains` of each fold
      stanfit[[k]] <- rstan::sflist2stanfit(sflist[inchains])
    }  
    return(stanfit) 
  }
  
  else{
    stanfit <- cmdstanr::read_cmdstan_csv(list.files(output_dir, ".csv$", full.names = TRUE),
                                          variables = pars)
    
    stanfit
  }
}

#extract log-likelihoods of held-out data
extract_log_lik_K <- function(list_of_stanfits, list_of_holdout, ...){
  require(loo)
  K <- length(list_of_stanfits)
  holdout_matrix <- simplify2array(list_of_holdout)
  Nrep <- sum(holdout_matrix[1,])
  K <- ncol(holdout_matrix)/Nrep
  list_of_log_liks <- lapply(seq_len(Nrep*K), function(k){
    extract_log_lik(list_of_stanfits[[k]],...)
  })
  # `log_lik_heldout` will include the loglike of all the held out data of all the folds.
  # We define `log_lik_heldout` as a (samples x N_obs) matrix
  # (similar to each log_lik matrix)
  N_obs <- dim(list_of_log_liks[[1]])[2]
  # log_lik_heldout <- list_of_log_liks[[1]] * NA
  log_lik_heldout <-
    do.call(rbind,
      lapply(seq_len(Nrep),
             function(n){
               log_lik_heldout_n <- list_of_log_liks[[1]] * NA
               for (k in ((n-1)*K+1):(n*K)){
                 log_lik <- list_of_log_liks[[k]]
                 samples <- dim(log_lik)[1] 
                 N_obs <- dim(log_lik)[2]
                 # This is a matrix with the same size as log_lik_heldout
                 # with 1 if the data was held out in the fold k
                 heldout <- matrix(rep(list_of_holdout[[k]], each = samples), nrow = samples)
                 # Sanity check that the previous log_lik is not being overwritten:
                 if(any(!is.na(log_lik_heldout_n[heldout==1]))){
                   warning("Heldout log_lik has been overwritten!!!!")
                 }
                 # We save here the log_lik of the fold k in the matrix:
                 log_lik_heldout_n[heldout==1] <- log_lik[heldout==1]
               }
               log_lik_heldout_n
             })
    )
  
  attr(log_lik_heldout, "K") <- K
  return(log_lik_heldout)
}

#compute ELPD
kfold <- function(log_lik_heldout)  {
  library(matrixStats)
  logColMeansExp <- function(x) {
    # should be more stable than log(colMeans(exp(x)))
    S <- nrow(x)
    colLogSumExps(x, na.rm=TRUE) - log(S)
  }
  # See equation (20) of @VehtariEtAl2016
  pointwise <-  matrix(logColMeansExp(log_lik_heldout), ncol= 1)
  colnames(pointwise) <- "elpd"
  # See equation (21) of @VehtariEtAl2016
  elpd_kfold <- sum(pointwise)
  se_elpd_kfold <-  sqrt(ncol(log_lik_heldout) * var(pointwise))
  out <- list(
    pointwise = structure(pointwise, dimnames = list(NULL, "elpd_kfold")),
    estimates = structure(cbind(elpd_kfold, se_elpd_kfold), dimnames = list('elpd_kfold', c("Estimate", "SE")))
  )
  attr(out, "K") <- attr(log_lik_heldout, "K")
  class(out) <- c("kfold", "loo")
  return(out)
}

extract_K_fold <- function(list_of_stanfits, list_of_holdouts, pars = NULL, ...,  holdout=TRUE){
  K = length(list_of_holdouts)
  Nrep = sum(simplify2array(list_of_holdouts)[1,])
  holdout <- as.numeric(holdout)
  stopifnot(length(list_of_stanfits) == K)
  D = if (holdout) 1 else (K/Nrep - 1)
  par_extract_list <- lapply(list_of_stanfits,FUN = rstan::extract, pars=pars, ...)
  extract_pars <- names(par_extract_list[[1]])
  # browser()
  extract_holdout <- lapply(extract_pars, function(p) {
    extract_holdout_par <- array(NA, dim = c(dim(par_extract_list[[1]][[p]])[1] * Nrep * D, dim(par_extract_list[[1]][[p]])[-1]))
    # browser()
    for (n in seq_len(Nrep)){
      for (k in seq_len(K2 <- K/Nrep)){
        holdout_k <- list_of_holdouts[[(n-1)*K2+k]]
        for (d in seq_len(D)){
          # print((n-1)*K2+k)
          k_par <- par_extract_list[[(n-1)*K2+k]][[p]] 
          #browser()
          par_dims <- dim(k_par)
          # browser()
          if (length(par_dims) >= 2 && par_dims[min(length(par_dims),2)] == length(holdout_k)) {
            commas <- paste(rep(',', length(par_dims)-2), collapse = " ")
            # browser()
            call_string <- glue::glue('extract_holdout_par[(dim(par_extract_list[[1]][[p]])[1]*(n-1)+1):(dim(par_extract_list[[1]][[p]])[1]*n),
                                    holdout_k=={holdout} {commas}] <- k_par[,holdout_k=={holdout} {commas}]')
            eval(parse(text=call_string))
          } else {
            commas <- paste(rep(',', length(par_dims)-1), collapse = " ")
            call_string <-- glue::glue('extract_holdout_par[[(dim(par_extract_list[[1]][[p]])[1]*(n-1)+1):(dim(par_extract_list[[1]][[p]])[1]*n) {commas}] <- 
                                      k_par[{commas}]')
            eval(parse(text=call_string))
          }
        }
      }
    }
    
    extract_holdout_par
  })
  setNames(extract_holdout, extract_pars)
}

 

# Function to draw the calibration curve
calib_curve <- function(pred, obs, title = NULL, method="loess", se = TRUE, ...){
  require(ggplot2)
  ggplot(mapping=aes(x = pred)) + 
    geom_smooth(aes(y = as.numeric(obs)), color="red", se = se, method = method, ...) + 
    # geom_ribbon(aes(ymin = lb, ymax=ub), alpha=.5, fill=grey(.6)) +
    geom_line(aes(y = pred)) +
    geom_point(aes(y = as.numeric(obs))) +
    ylab("obs") +
    ggtitle(title)
}


# Repeated K-fold
repeated_stan_K_fold <- function(file, data, K, N, repetition, cores, ...){
  require(future)
  fold_sets <- 
    lapply(seq_len(repetition), function(r) {
      hh <- loo::kfold_split_random(K, N) 
      keptins <- matrix(1, nrow = N, ncol = K)
      for(i in 1:N) keptins[i, hh[i]] <- 0
      holdouts  <- 1 - keptins
      keptins <- split(keptins,rep(1:ncol(keptins),each=nrow(keptins)))
      holdouts  <- split(holdouts,rep(1:ncol(holdouts),each=nrow(holdouts)))
      datas <- lapply(seq_len(K), function(k) modifyList(data, list(keptin=keptins[[k]])))
      list(keptins = keptins, holdouts = holdouts, datas = datas)
    })
  # browser()
  progressr::with_progress({
    p <- progressr::progressor(steps = length(repetition))
    plan(multisession, workers=max(1, floor(cores/K)))
    for (r in seq_len(repetition)){
      mc_cores <- min(cores, K)
      dt = fold_sets[[r]]
      dt$fits %<-% stan_kfold(file = file, list_of_datas = dt$datas, cores = mc_cores, ...)
      foldsets[[r]] <- dt
      p()
    }
  })
  
  class(fold_sets) <- c("rep_kfold", "kfold")
  fold_sets
} 
