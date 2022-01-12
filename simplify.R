#!/usr/bin/env Rscript

library(argparser)
misc <- new.env(hash = FALSE)
source("r/include/functions.R", local = misc)
rstan::rstan_options(javascript = FALSE, auto_write = TRUE)

args <- 
  arg_parser("TBM-LCA model simplifier - Written by Trinh Dong") |>
  add_argument("target", help = 'Target to simplify, must be an output of sample.R with --rep=1 and --fold=1') |>
  add_argument("--input", help = "input Rdata file", 
               default = "data/cleaned/data_input.Rdata", short = "-i") |>
  add_argument("--mode", help = "mode, either sampling, optimizing, or vb. Apart from sampling, only iteration number is tweakable", default = "sampling", short="-m") |>
  add_argument("--output-dir", help = "output folder", default = "outputs", short = "-o") |>
  add_argument("--output-file", help = "output filename, default to auto-inference", short = "-f") |>
  add_argument("--fold", help = "number of folds for cross-validation, ignored if not in a clean state and cache is available", default = 10, short = "-k") |>
  add_argument("--rep", help = "number of repetitions of K-fold procedure, ignored if not in a clean state and cache is available", default = 4, short = "-r") |>
  add_argument("--seed", help = "random seed, default to randomly select a number from 0 to 100,000", default = sample(0:100000, 1), short = "-S") |>
  add_argument("--clean-state", help = "ignore cache to run a clean state", flag = TRUE) |> 
  add_argument("--cache-dir", help = "cache directory, ignored if --no-cache flag is on", default = ".cache", short = "-d") |>
  add_argument("--no-cache", help = "disable caching and monitoring process, written to <cache-dir>/sampling/<model> ", flag = TRUE, short = "-n") |> 
  add_argument("--keep-cache", help = "disable clearing cache after running", flag = TRUE, short = "-l") |>
  add_argument("--chains", help = "number of chains per dataset", default = 3, short = "-N") |>
  add_argument("--cores", help = "number of cores", default = parallel::detectCores()-1, short = "-C") |>
  add_argument("--thin", help = "thining the chains", default = 1, short = "-T") |>
  add_argument("--iter", help = "number of iterations", default = 5000, short = "-I") |>
  add_argument("--warmup", help = "number of warm-up iterations", default = 2000, short = "-W") |> 
  add_argument("--init-r", help = "Stan init parameter", default = 1, short = "-R") |>
  add_argument("--adapt-delta", help = "Stan param", default = .8, short = "-D") |>
  add_argument("--penalty-family", help = "Family for the regularised prior, either 't', 'normal', or 'laplace'", default = "t", short = "-F") |>
  add_argument("--penalty-term", help = "Penalty scale term, if 0 then auto adaptation is performed, second term == -1 will make it follow the first", nargs = 1, default = 0, short = "-P") |>
  add_argument("--all-params", help = "Whether to expose imputation params, ignored for m0", flag = TRUE, short="-a") |>
  add_argument("--include-pars", help = "List of parameters to be extracted, default is based on the model", default = NA_character_, nargs = Inf) |>
  add_argument("--extra-x", help =  "Extra X", default = NA_character_, nargs = Inf) |>
  add_argument("--use-rstan-compiler", help = "DEBUG arg: by default, use a custom stan_model function, if TRUE, use the default rstan one", flag = TRUE) 

argparser <- parse_args(args)

create_folds <- 
  \(recipe, K, N, X_extra, seed, cache_file=NULL){
    recipe <- recipe
    inp <-
      with(recipe,
           list(
             N_all = nrow(data_19EI),
             nXc = ncol(Xc),
             nXd = ncol(Xd),
             nTd = ncol(Td),
             nTc = ncol(Tc),
             Y_Smear_all = data_19EI$csf_smear,
             Y_Mgit_all = data_19EI$csf_mgit,
             Y_Xpert_all = data_19EI$csf_xpert,
             obs_Smear_all = data_19EI$obs_smear,
             obs_Mgit_all = data_19EI$obs_mgit,
             obs_Xpert_all = data_19EI$obs_xpert,
             Xc_all = Xc,
             Xd_all = Xd,
             Td_all = Td,
             Tc_all = Tc,
             nX_extra = ncol(X_extra),
             X_extra_all = X_extra,
             obs_Xc_all = cbind(obs_Xc),
             obs_Xd_all = obs_Xd,
             obs_Td_all = obs_Td,
             obs_Tc_all = obs_Tc,
             penalty_family = penalty_family,
             penalty_term = penalty_term,
             mu_theta = mu_theta,
             mu_ztheta = mu_ztheta,
             sd_ztheta = sd_ztheta
           ))
    folds <- misc$repeated_kfold(inp, K = K, N_rep = N, N_obs = nrow(recipe$data_19EI), seed = seed)
    if (!is.null(cache_file)) saveRDS(folds, file = cache_file)
    folds
  }

argparser$model <- 's1'
results <- new.env(parent=emptyenv())
results$.META <- argparser

with(
  argparser,
  {
    
    my_stan_model <- if (use_rstan_compiler) rstan::stan_model else misc$my_stan_model
    if (is.na(output_file)) output_file <- paste(model,"RDS",sep=".")
    output_name <- fs::path_ext_remove(fs::path_file(output_file))
    
    if ((rep>1 || fold>1) && mode != 'sampling') 
      stop("Repeated K-Fold only support sampling.")
    
    penalty_family <- switch(penalty_family, 
                             "t" = 0,
                             "laplace" = 1,
                             "normal" = 2)
    if (is.na(penalty_family)) stop("Invalid penalty family. Must be either t, laplace, or normal.")
    # Load the recipe
    recipe <- new.env()
    load(input, envir = recipe)
    recipe$penalty_term = penalty_term
    recipe$penalty_family = penalty_family
    
    X_extra = if (!all(is.na(extra_x))) recipe$data_19EI[extra_x] else NULL
    
    # HIV missing case
    # if (hiv_missing == 0){
    #   recipe$obs_Xd[,1] = 1
    # }
    # if (hiv_missing == 1){
      recipe$obs_Xd[recipe$Td[,7]==0,1] = 1
      # recipe$obs_Xd[,1] = ifelse((recipe$Td[,7]==0)&(recipe$obs_Xd[,1]==0), 1, recipe$obs_Xd[,1])
    # }
    # if (hiv_missing == 2){
    #   recipe$Td[,7] = 1
    # }
    
    # If there are cache to use, create corresponding dirs
    if (!no_cache) {
      dir.create(file.path(cache_dir, "sampling"), recursive = TRUE, showWarnings = FALSE, mode = "770")
      dir.create(file.path(cache_dir, "folds"), recursive = TRUE, showWarnings = FALSE, mode = "770")
    }
    
    # Load cache or create new folds
    cache_file <- file.path(cache_dir, "folds", paste0(output_name, "_recipe.RDS"))
    has_cache <- file.exists(cache_file)
    
    prior_family_name <- 
      switch(as.character(penalty_family), 
             "0" = "Student t", "1" = "Laplace", "2" = "Normal")
    # Print out program information
    cli::cli_alert_info('Program starts at: {Sys.time()}')
    cli::cli_h1('Configuration:')
    cli::cli_ul()
    cli::cli_li('{.strong Model} {model}')
    cli::cli_li('Mode: {mode}')
    cli::cli_li('{fold} fold{?s} with {rep} repetition{?s}')
    cli::cli_li('{.strong Prior family:} {.field {prior_family_name}} with {.strong scales} = [{toString(ifelse(penalty_term==0, "~N[0,2.5]", penalty_term))}]')
    cli::cli_li('{.strong Random seed:} {.field {seed}}')
    cli::cli_li('Stan configurations:')
    ulid <- cli::cli_ul()
    cli::cli_li('Cores = {.field {cores}}')
    cli::cli_li('Chains = {.field {chains}}')
    cli::cli_li('Iters = {.field {iter}}')
    cli::cli_li('Warmup = {.field {warmup}}')
    cli::cli_li('adapt_delta = {.field {adapt_delta}}')
    cli::cli_li('init_r = {.field {init_r}}')
    cli::cli_end(ulid)
    cli::cli_end()
    cli::cli_h1('')
    
    cli::cli_alert("Extract posteriors from target")
    target_file <- file.path('outputs', paste0(target,'.RDS'))
    target_fit <- readRDS(target_file)
    if (!inherits(target_fit, 'stanfit')) target_fit <- target_fit$outputs 
    z_theta <- rstan::summary(target_fit, pars = c('z_theta'))$summary
    theta <- rstan::summary(target_fit, pars = c('theta'))$summary
    recipe$mu_ztheta <- z_theta[,'mean']
    recipe$sd_ztheta <- z_theta[,'sd']
    recipe$mu_theta <- theta[,'mean']
    
    if (clean_state){
      if (has_cache) {
        cli::cli_alert("Remove cache and create new folds")
        file.remove(cache_file)
      }
      folds <- create_folds(recipe=recipe, K=fold, N=rep, X_extra=X_extra, seed=seed, cache_file=cache_file)
    } else{
      if (has_cache) {
        cli::cli_alert("Cache file found. Use cache file.")
        folds <- readRDS(cache_file)
      } else {
        cli::cli_alert("No cache file found. Create new folds")
        folds <- create_folds(recipe=recipe, K=fold, N=rep, X_extra=X_extra, seed=seed, cache_file=cache_file)
      }
    }
    
    # Load the folds
    inputs <- folds$inputs
    
    # Remove the cache dir if exists
    if (!no_cache) {
      outdir <- file.path(cache_dir, "sampling", output_name)
      if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
      dir.create(outdir, showWarnings = F)
    } else outdir <- NULL
    
    cli::cli_alert('Compile model sampler')
    sampler <-"stan/s3.stan"
    sampler <- tryCatch(
      my_stan_model(sampler),
      error = function(e) {
        # Remove the cache file might fix the problem
        try(file.remove(file.path("stan", paste0(model, ".rds"))))
        try(file.remove(file.path("stan", paste0(model, ".RDS"))))
        my_stan_model(sampler)
      }
    )
    cli::cli_alert('Sample')
    pars <- c('a0', 'a', 'log_lik', 'ztheta_all','lambda')
    if (!all(is.na(include_pars))) pars <- include_pars
    if (fold == 1 && all_params) pars <- NA
    
    if (mode == "sampling"){
      results$outputs <- misc$stan_kfold(sampler = sampler,
                                         list_of_datas=inputs,
                                         backend = "rstan",
                                         chains = chains, cores = cores, 
                                         thin = thin, 
                                         merge = TRUE,
                                         control = list(adapt_delta=adapt_delta, max_treedepth=12), 
                                         init_r = init_r, seed = seed,
                                         sample_dir = outdir,
                                         pars = pars,
                                         iter=iter, warmup=warmup)
    } 
    
    if (mode == "optimizing"){
      results$outputs <- rstan::optimizing(object = sampler,
                                           verbose = TRUE,
                                           data = inputs[[1]],
                                           seed = seed,
                                           hessian = TRUE,
                                           draws = 1000,
                                           iter = iter)
    }
    
    if (mode == "vb"){
      results$outputs <- rstan::vb(object = sampler,
                                   data = inputs[[1]],
                                   seed = seed,
                                   iter = iter,
                                   elbo_samples = max(warmup, 100),
                                   adapt_iter = 200,
                                   tol_rel_obj = 1e-5)
    }
    
    # Clear cache
    if (!keep_cache){
      cli::cli_inform("Clean up cache")
      unlink(outdir, recursive = TRUE)
      file.remove(cache_file)
    } 
    
    # Save results
    cli::cli_alert("Save results")
    results$model_name <- model
    results$folds  <- folds
    if (mode != 'optimizing')
      results$.META$params <- if (fold == 1) results$outputs@model_pars else results$outputs[[1]]@model_pars
    saveRDS(results, file = file.path(output_dir, output_file))
    future::plan("sequential")
    cli::cli_alert_success('Sampling completed!')
  }
)