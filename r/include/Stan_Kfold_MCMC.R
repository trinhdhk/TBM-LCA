StanKfold_MCMC <- 
  R6::R6Class(
    "StanKfold_MCMC", 
    public = list(
      mc.cores = 1,
      initialize = function(K, N_Rep, inputs,
                            keptin, holdout,
                            stanfits = NULL, mc.cores = 1){
        N_chains <- length(stanfits[[1]]@stan_args)
        if (K * N_Rep != length(stanfits)) stop("Number of stanfits does not match with K and N_Rep!")
        for (i in seq_along(stanfits)) {
          fit <- stanfits[[i]]
          cat(crayon::green(glue::glue("Fit {i}\n")))
          check_all_diagnostics(fit, quiet = FALSE)
          if (length(fit@stan_args)!=N_chains) stop("Number of chains mismatched! List of stanfits is ill-formed!")
        }
        private$.inputs <- inputs
        private$.keptin <- keptin
        private$.holdout <- holdout
        private$.N_chains <- N_chains
        private$.K <- K
        private$.N_Rep <- N_Rep
        private$.stanfits <- stanfits
        self$mc.cores <- mc.cores
      },
      extract_log_lik = function(...){
        K <- private$.K
        N_Rep <- private$.N_Rep
        list_of_log_liks <- lapply(seq_len(N_Rep*K), function(k){
          loo::extract_log_lik(private$.stanfits[[k]],...)
        })
        list_of_holdout <- private$.holdout
        # `log_lik_heldout` will include the loglike of all the held out data of all the folds.
        # We define `log_lik_heldout` as a (samples x N_obs) matrix
        # (similar to each log_lik matrix)
        N_obs <- private$.N
        log_lik_heldout <-
          do.call(rbind,
                  lapply(seq_len(N_Rep),
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
        private$.log_lik <- log_lik_heldout
        return(log_lik_heldout)
      },
      estimate_elpd = function(){
        log_lik_heldout <- self$extract_log_lik()
        logColMeansExp <- function(x) {
          # should be more stable than log(colMeans(exp(x)))
          S <- nrow(x)
          matrixStats::colLogSumExps(x, na.rm=TRUE) - log(S)
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
        private$.elpd <- out
        return(out)
      },
      combine = function(){
        rstan::sflist2stanfit(private$.stanfits)
      },
      extract_kfold = function(pars, ...,mc.cores = self$mc.cores, group = c("test", "train", "all")){
        if (!length(private$.stanfits)) stop("No sample found. Run `sampling` first!")
        K <- private$.K
        N <- private$.N
        N_Rep <- private$.N_Rep
        group <- match.arg(group)
        if (group == "all") rstan::extract(rstan::sflist2stanfit(private$.stanfits))
        selection <- switch(group, 
                            train = private$.keptin,
                            test  = private$.holdout,
                            all   = lapply(private$.keptin, function(i) array(TRUE, dim = dim(private$.keptin[[1]]))))
        
        if (length(private$.stanfits) != K*N_Rep) 
          stop("Manipulated or empty samples!")
        
        raw_extracted <- pbmcapply::pbmclapply(private$.stanfits, FUN = rstan::extract, pars = pars, ..., mc.cores = mc.cores)
        N_Chains <- private$.N_chains
        N_Iters <- dim(raw_extracted[[1]][[1]])[1]
        pars_extracted <- names(raw_extracted[[1]]) #actual extracted pars might not agree with the pars in parameters
        # browser()
        holdout_extracted <- 
          lapply(
            pars_extracted, 
            FUN = function(p) {
              p_extracted <- lapply(raw_extracted, `[[`, p)
              p2 <- vector("list", length(p_extracted))
              selection_matrix <- lapply(selection, function(h) do.call(rbind, lapply(seq_len(N_Iters), function(.) h)))
              for (i in seq_along(p_extracted)){
                p_extracted[[i]][as.logical(selection_matrix[[i]])] <- NA 
                p2[[i]] <- apply(p_extracted[[i]], 1, na.omit)
              }
              abind::abind(p2, along = 1)
            })
        setNames(holdout_extracted, pars_extracted)
      },
      summary_kfold = function(pars, ..., group = c("test", "train", "all"), .custom_fn){
        default <- missing(.custom_fn)
        extracted_pars <- self$extract_kfold(pars, group = group, ...)
        if (default){
          summary <- sapply(extracted_pars,
                            apply, 2, function(l) data.frame(mean = mean(l), median = median(l), CrI2.5 = quantile(l, .025), CrI97.5 = quantile(l, .975)), 
                            simplify = FALSE)
          summary <- sapply(summary, dplyr::bind_rows, simplify = FALSE, USE.NAMES = TRUE)
          summary <- sapply(summary, `rownames<-`, NULL)
        } else summary <- sapply(extracted_pars, .custom_fn)
        
        summary
      },
      diagnostic_plot = function(predicted, observed, kind = c("calibration", "roc"), group = c("test", "train", "all"), 
                                 .estimation_fn = mean, title = NULL, ..., mc.cores = self$mc.cores){
        kind <- match.arg(kind)
        group <- match.arg(group)
        pred <- self$summary_kfold(predicted, group = group, mc.cores = mc.cores)
        Dx_plot <- switch(kind,
                          calibration = private$.calibration_plot(pred, obs, title=title,...),
                          roc = private$.roc_plot(pred, obs, title=title,...))
        Dx_plot
      }
    ),
    private = list(
      .K = integer(), .N_Rep = integer(), .inputs = NULL,
      .N_chains = integer(), .stanfits = NULL,
      .elpd = NULL, .log_lik = NULL, 
      .keptin = NULL, .holdout = NULL,
      .calibration_plot = function(pred, obs, title = NULL, method="loess", se = TRUE, ...){
        require(ggplot2)
        ggplot(mapping=aes(x = pred)) + 
          geom_smooth(aes(y = as.numeric(obs)), color="red", se = se, method = method, ...) + 
          geom_line(aes(y = pred)) +
          geom_point(aes(y = as.numeric(obs))) +
          xlab("Predicted") +
          ylab("Observed") +
          ggtitle(title)
      },
      .roc_plot = function(pred, obs, title = NULL, ...){
        classifierplots::roc_plot(pred, obs, ...) + 
          ggplot2::ggtitle(title)
      }
    ),
    active = list(
      log_lik = function(...){
        if (length(private$.log_lik)) return(private$.log_lik)
        self$extract_log_lik(...)
      },
      elpd = function(){
        if (length(private$.elpd)) return(private$.elpd)
        self$estimate_elpd()
      },
      draws = function(){
        abind::abind(private$.stanfits)
      }
    )
  )


# Check transitions that ended with a divergence
check_div <- function(fit, quiet=FALSE) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  N = length(divergent)
  
  if (!quiet) print(sprintf('%s of %s iterations ended with a divergence (%s%%)',
                            n, N, 100 * n / N))
  if (n > 0) {
    if (!quiet) print('  Try running with larger adapt_delta to remove the divergences')
    if (quiet) return(FALSE)
  } else {
    if (quiet) return(TRUE)
  }
}

# Check transitions that ended prematurely due to maximum tree depth limit
check_treedepth <- function(fit, max_depth = 10, quiet=FALSE) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  N = length(treedepths)
  
  if (!quiet)
    print(sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%)',
                  n, N, max_depth, 100 * n / N))
  
  if (n > 0) {
    if (!quiet) print('  Run again with max_treedepth set to a larger value to avoid saturation')
    if (quiet) return(FALSE)
  } else {
    if (quiet) return(TRUE)
  }
}

# Checks the energy fraction of missing information (E-FMI)
check_energy <- function(fit, quiet=FALSE) {
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  no_warning <- TRUE
  for (n in 1:length(sampler_params)) {
    energies = sampler_params[n][[1]][,'energy__']
    numer = sum(diff(energies)**2) / length(energies)
    denom = var(energies)
    if (numer / denom < 0.2) {
      if (!quiet) print(sprintf('Chain %s: E-FMI = %s', n, numer / denom))
      no_warning <- FALSE
    }
  }
  if (no_warning) {
    if (!quiet) print('E-FMI indicated no pathological behavior')
    if (quiet) return(TRUE)
  } else {
    if (!quiet) print('  E-FMI below 0.2 indicates you may need to reparameterize your model')
    if (quiet) return(FALSE)
  }
}

# Checks the effective sample size per iteration
check_n_eff <- function(fit, quiet=FALSE) {
  fit_summary <- rstan::summary(fit, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]
  
  iter <- dim(rstan:::extract(fit)[[1]])[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    if (is.nan(fit_summary[,'n_eff'][n])) {
      if (!quiet) print(sprintf('n_eff for parameter %s is NaN!',
                                rownames(fit_summary)[n]))
      no_warning <- FALSE
    } else {
      ratio <- fit_summary[,'n_eff'][n] / iter
      if (ratio < 0.001) {
        if (!quiet) print(sprintf('n_eff / iter for parameter %s is %s!',
                                  rownames(fit_summary)[n], ratio))
        no_warning <- FALSE
      }
    }
  }
  if (no_warning) {
    if (!quiet) print('n_eff / iter looks reasonable for all parameters')
    if (quiet) return(TRUE)
  }
  else {
    if (!quiet) print('  n_eff / iter below 0.001 indicates that the effective sample size has likely been overestimated')
    if (quiet) return(FALSE)
  }
}

# Checks the potential scale reduction factors
check_rhat <- function(fit, quiet=FALSE) {
  fit_summary <- rstan::summary(fit, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]
  
  no_warning <- TRUE
  for (n in 1:N) {
    rhat <- fit_summary[,'Rhat'][n]
    if (rhat > 1.1 || is.infinite(rhat) || is.nan(rhat)) {
      if (!quiet) print(sprintf('Rhat for parameter %s is %s!',
                                rownames(fit_summary)[n], rhat))
      no_warning <- FALSE
    }
  }
  if (no_warning) {
    if (!quiet) print('Rhat looks reasonable for all parameters')
    if (quiet) return(TRUE)
  } else {
    if (!quiet) print('  Rhat above 1.1 indicates that the chains very likely have not mixed')
    if (quiet) return(FALSE)
  }
}

check_all_diagnostics <- function(fit, quiet=FALSE) {
  if (!quiet) {
    check_n_eff(fit)
    check_rhat(fit)
    check_div(fit)
    check_treedepth(fit)
    check_energy(fit)
  } else {
    warning_code <- 0
    
    if (!check_n_eff(fit, quiet=TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 0))
    if (!check_rhat(fit, quiet=TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 1))
    if (!check_div(fit, quiet=TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 2))
    if (!check_treedepth(fit, quiet=TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 3))
    if (!check_energy(fit, quiet=TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 4))
    
    return(warning_code)
  }
}

parse_warning_code <- function(warning_code) {
  if (bitwAnd(warning_code, bitwShiftL(1, 0)))
    print("n_eff / iteration warning")
  if (bitwAnd(warning_code, bitwShiftL(1, 1)))
    print("rhat warning")
  if (bitwAnd(warning_code, bitwShiftL(1, 2)))
    print("divergence warning")
  if (bitwAnd(warning_code, bitwShiftL(1, 3)))
    print("treedepth warning")
  if (bitwAnd(warning_code, bitwShiftL(1, 4)))
    print("energy warning")
}

# Returns parameter arrays separated into divergent and non-divergent transitions
partition_div <- function(fit) {
  nom_params <- rstan:::extract(fit, permuted=FALSE)
  n_chains <- dim(nom_params)[2]
  params <- as.data.frame(do.call(rbind, lapply(1:n_chains, function(n) nom_params[,n,])))
  
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  params$divergent <- divergent
  
  div_params <- params[params$divergent == 1,]
  nondiv_params <- params[params$divergent == 0,]
  
  return(list(div_params, nondiv_params))
}
