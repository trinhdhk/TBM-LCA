#Function to create multi mcmc_intervals
mcmc_intervals_multi <-
  \(fits = list(),..., .width=.5){
    require(ggplot2)
    w <- .15 * length(fits)
    if (!length(names(fits))) names(fits) <- as.character(as.list(substitute(fits)[-1]))
    est  <- sapply(seq_along(fits),
                   \(i) {
                     bayesplot::mcmc_intervals_data(fits[[i]], ...) |> 
                       dplyr::mutate(Model = names(fits)[[i]])
                   }, simplify=FALSE, USE.NAMES = TRUE)
    est_dat <- dplyr::bind_rows(est) |> dplyr::mutate(Model = factor(Model, levels = names(fits)))
    ggplot(est_dat, aes(y=parameter, color=Model)) +
      geom_linerange(aes(xmin=ll, xmax=hh), position=position_dodge(width=w)) +
      geom_linerange(aes(xmin=l, xmax=h), size=1, position=position_dodge(width=w)) +
      geom_point(aes(x=m), size=1.2, position=position_dodge(width=w))
  }

#fix rstan check for date when recompiling
my_stan_model <- 
  function (file, model_name = "anon_model", model_code = "", 
            stanc_ret = NULL, boost_lib = NULL, eigen_lib = NULL, save_dso = TRUE, 
            verbose = FALSE, auto_write = rstan_options("auto_write"), 
            obfuscate_model_name = TRUE, allow_undefined = FALSE, allow_optimizations = FALSE, 
            standalone_functions = FALSE, use_opencl = FALSE, warn_pedantic = FALSE, 
            warn_uninitialized = FALSE, includes = NULL, isystem = c(if (!missing(file)) dirname(file), 
                                                                     getwd())) 
  {
    if (isTRUE(rstan_options("threads_per_chain") > 1L)) {
      Sys.setenv(STAN_NUM_THREADS = rstan_options("threads_per_chain"))
    }
    if (is.null(stanc_ret)) {
      model_name2 <- deparse(substitute(model_code))
      if (is.null(attr(model_code, "model_name2"))) 
        attr(model_code, "model_name2") <- model_name2
      if (missing(model_name)) 
        model_name <- NULL
      if (missing(file)) {
        tf <- tempfile()
        writeLines(model_code, con = tf)
        file <- file.path(dirname(tf), paste0(tools::md5sum(tf), 
                                              ".stan"))
        if (!file.exists(file)) 
          file.rename(from = tf, to = file)
        else file.remove(tf)
      }
      else file <- normalizePath(file)
      stanc_ret <- stanc(file = file, model_code = model_code, 
                         model_name = model_name, verbose = verbose, obfuscate_model_name = obfuscate_model_name, 
                         allow_undefined = allow_undefined, 
                         standalone_functions = standalone_functions, use_opencl = use_opencl, 
                         warn_pedantic = warn_pedantic, warn_uninitialized = warn_uninitialized, 
                         isystem = isystem)
      model_re <- "(^[[:alnum:]]{2,}.*$)|(^[A-E,G-S,U-Z,a-z].*$)|(^[F,T].+)"
      if (!is.null(model_name)) 
        if (!grepl(model_re, model_name)) 
          stop("model name must match ", model_re)
      S4_objects <- apropos(model_re, mode = "S4", ignore.case = FALSE)
      if (length(S4_objects) > 0) {
        e <- environment()
        stanfits <- sapply(mget(S4_objects, envir = e, inherits = TRUE), 
                           FUN = is, class2 = "stanfit")
        stanmodels <- sapply(mget(S4_objects, envir = e, 
                                  inherits = TRUE), FUN = is, class2 = "stanmodel")
        if (any(stanfits)) 
          for (i in names(which(stanfits))) {
            obj <- get_stanmodel(get(i, envir = e, inherits = TRUE))
            if (identical(obj@model_code[1], stanc_ret$model_code[1])) 
              return(obj)
          }
        if (any(stanmodels)) 
          for (i in names(which(stanmodels))) {
            obj <- get(i, envir = e, inherits = TRUE)
            if (identical(obj@model_code[1], stanc_ret$model_code[1])) 
              return(obj)
          }
      }
      mtime <- file.info(file)$mtime
      file.rds <- gsub("stan$", "rds", file)
      md5 <- tools::md5sum(file)
      if (!file.exists(file.rds)) {
        file.rds <- file.path(tempdir(), paste0(md5, ".rds"))
      }
      if (!file.exists(file.rds) || isTRUE((mtime.rds <- file.info(file.rds)$mtime) < 
                                           as.POSIXct(packageDescription("rstan")$Date)) || 
          !is(obj <- readRDS(file.rds), "stanmodel") || !is_sm_valid(obj) || 
          (!identical(stanc_ret$model_code, obj@model_code) && 
           is.null(message("hash mismatch so recompiling; make sure Stan code ends with a blank line"))) || 
          avoid_crash(obj@dso@.CXXDSOMISC$module) && is.null(message("recompiling to avoid crashing R session"))) {
      }
      else return(invisible(obj))
    }
    if (!is.list(stanc_ret)) {
      stop("stanc_ret needs to be the returned object from stanc.")
    }
    m <- match(c("cppcode", "model_name", "status"), names(stanc_ret))
    if (any(is.na(m))) {
      stop("stanc_ret does not have element `cppcode', `model_name', and `status'")
    }
    else {
      if (!stanc_ret$status) 
        stop("stanc_ret is not a successfully returned list from stanc")
    }
    if (.Platform$OS.type != "windows") {
      CXX <- get_CXX()
      if (!is.null(attr(CXX, "status")) || nchar(CXX) == 0) {
        WIKI <- "https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started"
        warning(paste("C++ compiler not found on system. If absent, see\n", 
                      WIKI))
      }
      else if (grepl("69", CXX, fixed = TRUE)) 
        warning("You may need to launch Xcode once to accept its license")
    }
    else CXX <- "g++"
    model_cppname <- stanc_ret$model_cppname
    model_name <- stanc_ret$model_name
    model_code <- stanc_ret$model_code
    model_cppcode <- stanc_ret$cppcode
    model_cppcode <- paste("#ifndef MODELS_HPP", "#define MODELS_HPP", 
                           "#define STAN__SERVICES__COMMAND_HPP", "#include <rstan/rstaninc.hpp>", 
                           model_cppcode, "#endif", sep = "\n")
    inc <- paste("#include <Rcpp.h>\n", "using namespace Rcpp;\n", 
                 if (is.null(includes)) 
                   model_cppcode
                 else sub("(class[[:space:]]+[A-Za-z_][A-Za-z0-9_]*[[:space:]]*)", 
                          paste(includes, "\\1"), model_cppcode), "\n", get_Rcpp_module_def_code(model_cppname), 
                 sep = "")
    if (verbose && interactive()) 
      cat("COMPILING THE C++ CODE FOR MODEL '", model_name, 
          "' NOW.\n", sep = "")
    if (verbose) 
      cat(system_info(), "\n")
    if (!is.null(boost_lib)) {
      old.boost_lib <- rstan_options(boost_lib = boost_lib)
      on.exit(rstan_options(boost_lib = old.boost_lib))
    }
    if (!file.exists(rstan_options("boost_lib"))) 
      stop("Boost not found; call install.packages('BH')")
    if (!is.null(eigen_lib)) {
      old.eigen_lib <- rstan_options(eigen_lib = eigen_lib)
      on.exit(rstan_options(eigen_lib = old.eigen_lib), add = TRUE)
    }
    if (!file.exists(rstan_options("eigen_lib"))) 
      stop("Eigen not found; call install.packages('RcppEigen')")
    dso <- cxxfunctionplus(signature(), body = paste(" return Rcpp::wrap(\"", 
                                                     model_name, "\");", sep = ""), includes = inc, plugin = "rstan", 
                           save_dso = save_dso | auto_write, module_name = paste("stan_fit4", 
                                                                                 model_cppname, "_mod", sep = ""), verbose = verbose)
    obj <- new("stanmodel", model_name = model_name, model_code = model_code, 
               dso = dso, mk_cppmodule = mk_cppmodule, model_cpp = list(model_cppname = model_cppname, 
                                                                        model_cppcode = model_cppcode))
    if (missing(file) || (file.access(dirname(file), mode = 2) != 
                          0) || !isTRUE(auto_write)) {
      tf <- tempfile()
      writeLines(model_code, con = tf)
      file <- file.path(tempdir(), paste0(tools::md5sum(tf), 
                                          ".stan"))
      if (!file.exists(file)) 
        file.rename(from = tf, to = file)
      else file.remove(tf)
      saveRDS(obj, file = gsub("stan$", "rds", file))
    }
    else if (isTRUE(auto_write)) {
      file <- gsub("stan$", "rds", file)
      if (file.exists(file)) {
        rds <- try(readRDS(file), silent = TRUE)
        if (!is(rds, "stanmodel")) 
          warning(rds, " exists but is not a 'stanmodel' so not overwriting")
        else saveRDS(obj, file = file)
      }
      else saveRDS(obj, file = file)
    }
    invisible(obj)
  }
environment(my_stan_model) <- environment(rstan::stan_model)

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
  if (K == 1) return(
    list(keptin = list(rep(1L, N_obs)), holdout = list(rep(0L, N_obs)), 
         inputs = list(modifyList(data, list(keptin = rep(1L, N_obs)))))
  )
  holdout <- vector("list", N_rep * K)
  keptin <- vector("list", N_rep * K)
  for(n in seq_len(N_rep)){
    set.seed(seed + n - 1)
    hh <- loo::kfold_split_random(K = K, N = N_obs)
    foldkept <- matrix(1, nrow = N_obs, ncol = K)
    for(i in 1L:N_obs) foldkept[i, hh[i]] <- 0
    foldhold  <- 1L - foldkept
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
normal_stan <- function(model, data, include_paths=NULL, sample_dir = NULL, backend, chains, cores, seed, pars,...){
  sample_file <- file.path(sample_dir, 'data_chain')
  if (backend == 'rstan') {
    sf <- rstan::sampling(model, data=data, sample_file = sample_file, chains=chains, cores=cores, seed=seed, pars=pars,... )
  } else {
    m <- model$clone()
    sf <- m$sample(data = data,
                   chains = 1, seed = seed,
                   output_dir = sample_dir,...)
  }
  
  sf
}

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
  if (length(list_of_datas) == 1L) 
    return(normal_stan(model, list_of_datas[[1]], include_path, sample_dir, backend,chains, cores, seed, pars,...))
  # First parallelize all chains:
  future::plan(future::multicore, workers=cores, gc=TRUE)
  wd <- getwd()
  progressr::handlers(progressr::handler_progress)
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
                      chains = 1, seed = seed, chain_ids = i,
                      output_dir = sample_dir,...)
      }
      p()
      sf
      # }, seed=T, gc=T)
    }, .options = furrr::furrr_options(seed=TRUE))
  }, enable = TRUE)
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
            call_string <- glue::glue('extract_holdout_par[(dim(par_extract_list[[1]][[p]])[1]*(n-1)*d+1):(dim(par_extract_list[[1]][[p]])[1]*n*d),
                                    holdout_k=={holdout} {commas}] <- k_par[,holdout_k=={holdout} {commas}]')
            eval(parse(text=call_string))
          } else {
            stop('Par(s) is/are not individual propert(y/ies). Please use the rstan::extract instead!')
            #commas <- paste(rep(',', length(par_dims)-1), collapse = " ")
            #call_string <-- glue::glue('extract_holdout_par[[(dim(par_extract_list[[1]][[p]])[1]*(n-1)+1):(dim(par_extract_list[[1]][[p]])[1]*n) {commas}] <- 
            #                          k_par[{commas}]')
            #eval(parse(text=call_string))
          }
        }
      }
    }
    
    extract_holdout_par
  })
  setNames(extract_holdout, extract_pars)
}

 # Thin histogram
thinhist_subplot <- function(x, normalize.factor = NULL, digits = 2, yscale = 1, zero_y = 0, plot = TRUE){
  round_factor = 10^digits
  sub = floor(min(x*round_factor))/round_factor
  sup = ceiling(max(x*round_factor))/round_factor
  breaks = seq(sub, sup, 10^-digits)
  x.breaks = hist(x, breaks, plot = FALSE)
  if (!length(normalize.factor)) normalize.factor <- max(x.breaks$count)
  counts = (x.breaks$count / normalize.factor) * yscale + zero_y
  require(ggplot2)
  p <- geom_segment(aes(y = zero_y, x = x.breaks$mids, yend = counts, xend = x.breaks$mids))
  if (!plot) return(p)
  ggplot() + p  +
    ylab("") + xlab("")
}

thinhist_subplot.binary <- function(x, y, normalize = TRUE, digits = 2, yscale = 1, plot = TRUE){
  round_factor = 10^digits
  sub = floor(min(x*round_factor))/round_factor
  sup = ceiling(max(x*round_factor))/round_factor
  breaks = seq(sub, sup, 10^-digits)
  x.breaks = hist(x, breaks, plot = FALSE)
  normalize.factor <- if (normalize) quantile(x.breaks$count, .9, na.rm=TRUE) else NULL
  p <-
    list(
      thinhist_subplot(x[y==0], normalize.factor=normalize.factor, digits = digits, yscale = yscale, zero_y=0, plot=FALSE),
      thinhist_subplot(x[y==1], normalize.factor=normalize.factor, digits = digits, yscale = -yscale, zero_y=1, plot=FALSE)  
    )
  if (!plot) return(p)
  ggplot() + p[[1]] + p[[2]] + 
    ylab("") + xlab("")
}

# Function to draw the calibration curve
calib_curve <- function(pred, pred_rep = NULL, obs, title = NULL, method=c("loess","splines"), se = TRUE, knots=3, hist = TRUE, hist.normalize = TRUE, ..., theme = ggplot2::theme_bw){
  require(ggplot2)
  maincolor <- "#CD113B"
  subcolor1 <- "#111111"
  subcolor2 <- "#999999"
  
  method <- match.arg(method) 
  # browser()
  p <- ggplot(mapping=aes(x = pred))
  add_pred <- function(pr, line = FALSE){
    pr <- pr
    if (method == "loess"){
      geom_smooth(aes(x = pr, y = as.numeric(obs)), color = subcolor1, fill = subcolor2, alpha=if (line) 1 else .5/length(pred_rep), size=line, se = !line, method = method,...)
    } else {
      stat_smooth(aes(x = pr, y = as.numeric(obs)), method="glm", formula=y~splines::ns(x,knots), color = subcolor1, fill = subcolor2, size=line, alpha=if (line) 1 else .6/length(pred_rep), se = !line, ...) 
    }
  }
  
  for (pr in pred_rep) p <- p + add_pred(pr)
  for (pr in pred_rep) p <- p + add_pred(pr, line = .3)
    
  
  if (method == "loess"){
	  p <- p + geom_smooth(aes( y = as.numeric(obs)), color=maincolor, fill = subcolor2, se = se, size=.9, alpha=.3, method = method, ...)
  }
  else {
    p <- p + stat_smooth(aes( y = as.numeric(obs)), fullrange = TRUE, method="glm", formula=y~splines::ns(x,knots), color = maincolor, fill = subcolor2, size=.9, alpha = .6/length(pred_rep), se = se,...) 
  }
  p <- p + 
    geom_line(aes(y = pred), color = "#000000") + 
    ylab("Observed") + xlab("Predicted")+
    scale_y_continuous(breaks = seq(0, 1, length.out = 5), limits = c(0,1), oob = scales::squish)+
    scale_x_continuous(breaks = seq(0, 1, length.out = 5))+
    ggtitle(title)
  
  if (hist){
    p0 <-  thinhist_subplot.binary(x = pred, y = as.numeric(obs), normalize = hist.normalize, yscale = .1, plot = FALSE)
    p <- p + p0[[1]] + p0[[2]]
  }
  
  p + theme()
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