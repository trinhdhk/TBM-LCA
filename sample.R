#!/usr/bin/env Rscript

library(argparser)
library(magrittr, include.only = "%>%")
misc <- new.env(hash = FALSE)
source("r/include/functions.R", local = misc)
rstan::rstan_options(javascript = FALSE, auto_write = TRUE)

args <- 
  arg_parser("TBM-LCA model runner - Written by Trinh Dong") %>%
  add_argument("model", help = "model to run") %>%
  # add_argument("--help", help = "print this help page and exit", short = "-h", flag = TRUE)
  add_argument("--input", help = "input Rdata file", 
               default = "data/cleaned/data_input.Rdata", short = "-i") %>%
  add_argument("--sampler", help = "sampler file, default to stan/<model>.stan", short = "-s") %>% 
  add_argument("--mode", help = "mode, either sampling, optimizing, or vb. Apart from sampling, only iteration number is tweakable", default = "sampling", short="-m") %>%
  add_argument("--output-dir", help = "output folder", default = "outputs", short = "-o") %>%
  add_argument("--output-file", help = "output filename, default to auto-inference", short = "-f") %>%
  add_argument("--fold", help = "number of folds for cross-validation, ignored if not in a clean state and cache is available", default = 10, short = "-k") %>%
  add_argument("--rep", help = "number of repetitions of K-fold procedure, ignored if not in a clean state and cache is available", default = 4, short = "-r") %>%
  add_argument("--seed", help = "random seed", default = sample(0:100000, 1), short = "-S") %>%
  add_argument("--clean-state", help = "ignore cache to run a clean state", flag = TRUE) %>% 
  add_argument("--cache-dir", help = "cache directory, ignored if --no-cache flag is on", default = ".cache", short = "-d") %>%
  add_argument("--no-cache", help = "disable caching and monitoring process, written to <cache-dir>/sampling/<model> ", flag = TRUE, short = "-n") %>% 
  add_argument("--keep-cache", help = "disable clearing cache after running", flag = TRUE, short = "-l") %>%
  add_argument("--chains", help = "number of chains per dataset", default = 3, short = "-N") %>%
  add_argument("--cores", help = "number of cores", default = parallel::detectCores()-1, short = "-C") %>%
  add_argument("--thin", help = "thining the chains", default = 1, short = "-T") %>%
  add_argument("--iter", help = "number of iterations", default = 5000, short = "-I") %>%
  add_argument("--warmup", help = "number of warm-up iterations", default = 2000, short = "-W") %>% 
  add_argument("--init-r", help = "Stan init parameter", default = 1, short = "-R") %>%
  add_argument("--adapt-delta", help = "Stan param", default = .8, short = "-D") %>%
  add_argument("--penalty-family", help = "Family for the regularised prior, either 't', 'normal', or 'laplace'", default = "t", short = "-F") %>%
  add_argument("--penalty-term", help = "Penalty scale term, if 0 then auto adaptation is performed, second term == -1 will make it follow the first", nargs = 2, default = c(0,0), short = "-P") %>%
  add_argument("--all-params", help = "Whether to expose imputation params, ignored for m0", flag = TRUE, short="-a") %>%
  add_argument("--quad-Xc", help = "Add quadratic effect for continuous variable, number are the positions in Xc", default = NA, nargs = Inf, short = '-q') %>%
  add_argument("--pca-n-FA", help = "Number of latent factors for PCA, 0 for no PCA", default = 0, short = "-P") %>%
  add_argument("--re-b", help = "Position of continuous variable to be included as predictor in bacillary burden model [0=none, 1=age, 2=illness day, 3=blood glucose, 4=csf glucose, 5=csf lympho, 6=csf protein, 7=csf lacate, 8=csf neutro, 9=gcs, 10=csf eos, 11=csf rbc, 12+=quadratic variables if exists]", nargs = Inf, default = NA, short = "-B") %>%
  add_argument("--include-d", help = "Add fixed outer effect for the bacillary burden", flag=TRUE) %>%
  add_argument("--pos-a", help = "Position of positive coefficients [0=none, 1=hiv, 2=tb symptoms, 3=motor palsy, 4=nerve palsy, 5=past TB contact, 6=xray PTB, 7=xray MTB, 8=crypto, 9=age, 10=illness day, 11=blood glucose, 12=csf glucose, 13=csf lympho, 14=csf protein, 15=csf lacate, 16=csf neutro, 17=gcs, 18=csf eos, 19=csf rbc, 20+=quadratic variables if exists", nargs = Inf, default = NA, short = '-p') %>%
  add_argument("--neg-a", help = "Position of negative coefficients, follow the same conventions as pos-a", nargs = Inf, default = NA, short = '-g') %>%
  add_argument("--quad-RE", help = "Sensitivity analysis for quadratic effect of Random Effect", flag=TRUE, short = "-Q") %>%
  add_argument("--lifted-spc", help = "Wider prior for test specificities (variances all set to .7 (default is .7 for Xpert, .3 for Mgit and Smear).", flag = TRUE, short = "-U") %>%
  add_argument("--use-rstan-compiler", help = "DEBUG arg: by default, use a custom stan_model function, if TRUE, use the default rstan one", flag = TRUE) %>%
  add_argument("--include-pars", help = "List of parameters to be extracted, default is based on the model", default = NA_character_, nargs = Inf) %>%
  add_argument("--hiv-missing", help = "HIV missing treating condition, 0: all is 0, 1: MAR for suspected, 2: MAR for all ", default=1, nargs=1)
argparser <- parse_args(args)

# Functions to create folds
create_folds <- function(recipe, K, N, seed, cache_file=NULL, n_FA, B, lifted_spc, quad_RE){
  inp <- 
    with(recipe,
         list(
           N_all = nrow(data_19EI),
           unsure_spc = lifted_spc,
           nB = length(B),
           B = B,# as.array(7), #vector(),
           quad_RE = quad_RE,
           nFA = n_FA,
           nA_pos = length(pos_a),
           nA_neg = length(neg_a),
           A_pos = pos_a,
           A_neg = neg_a,
           nXc = ncol(Xc),
           nXd = ncol(Xd),
           nTd = ncol(Td),
           nTc = ncol(Tc),
           nD  = (include_d) * ncol(D),
           D_all = if (include_d) D else array(dim=c(nrow(data_19EI), 0)),
           Y_Smear_all = data_19EI$csf_smear,
           Y_Mgit_all = data_19EI$csf_mgit,
           Y_Xpert_all = data_19EI$csf_xpert,
           Y_NegCrypto_all = recipe$csf_NegCrypto,
           obs_Smear_all = data_19EI$obs_smear,
           obs_Mgit_all = data_19EI$obs_mgit,
           obs_Xpert_all = data_19EI$obs_xpert,
           Xc_all = Xc,
           Xd_all = Xd,
           Td_all = Td,
           Tc_all = Tc,
           obs_Xc_all = cbind(obs_Xc),
           obs_Xd_all = obs_Xd,
           obs_Td_all = obs_Td,
           obs_Tc_all = obs_Tc,
           nQ = length(quad_Xc),
           quad_Xc_idx = as.array(quad_Xc),
           penalty_family = penalty_family,
           penalty_term = penalty_term
         )
    )

  folds <- misc$repeated_kfold(inp, K = K, N_rep = N, N_obs = nrow(recipe$data_19EI), seed = seed)
  if (!is.null(cache_file)) saveRDS(folds, file = cache_file)
  folds
}


results <- new.env(parent=emptyenv())
results$.META <- argparser
with(
  argparser,
  {
    # Stan compiler
    my_stan_model <- if (use_rstan_compiler) rstan::stan_model else misc$my_stan_model
    # If no name for output file, do auto-inference
    model_no <- as.numeric(substr(model, 2,2))
    if (pca_n_FA == 0 && grepl('_pca', model)) stop('Number of latent factor for PCA must be positive!')
    is_pca <- grepl('_pca', model) || pca_n_FA > 0
    if (is_pca && !grepl('_pca', model)) model <- paste0(model, '_pca')
    if (is.na(output_file)) output_file <- paste(model,"RDS",sep=".")
    output_name <- fs::path_ext_remove(fs::path_file(output_file))
    stopifnot(!(length(intersect(pos_a, neg_a)) > 0
                && max(length(pos_a), length(neg_a)) > 0))
    if ((rep>1 || fold>1) && mode != 'sampling') 
      stop("Repeated K-Fold only support sampling.")
    # Coerce numeric
    # fold <- as.integer(fold)
    # rep  <- as.integer(rep)
    # chains <- as.integer(chains)
    # cores  <- as.integer(cores)
    # seed   <- as.integer(seed)
    # iter   <- as.integer(iter)
    # warmup <- as.integer(warmup)
    # init_r <- as.numeric(init_r)
    # penalty_term <- as.numeric(penalty_term)
    pos_a <- if (all(is.na(pos_a))) 0 else pos_a
    neg_a <- if (all(is.na(neg_a))) 0 else neg_a
    re_b <- if (all(is.na(re_b))) 0 else re_b
    if (penalty_term[1] < 0 || penalty_term[2] < -1) stop('Invalid penalty term(s)')
    n_FA <- if (any(!pca_n_FA %in% 0)) as.vector(pca_n_FA) else vector()
    B <- if (any(!re_b %in% 0)) as.array(as.numeric(re_b)) else vector()
    pos_a <- if (any(!pos_a %in% 0)) as.array(as.numeric(pos_a)) else vector()
    neg_a <- if (any(!neg_a %in% 0)) as.array(as.numeric(neg_a)) else vector()
    lifted_spc <- as.numeric(lifted_spc)
    quad_Xc <- as.integer(quad_Xc)
    if(any(is.na(quad_Xc))) quad_Xc <- NULL 
    
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
    # Remove crypto if model == m3e
    if (model == 'm3e') {
      recipe$csf_NegCrypto = !recipe$Xd[,8]
      recipe$Xd = recipe$Xd[,-8]
      recipe$obs_Xd = recipe$obs_Xd[,-8]
    }
    
    # HIV missing case
    if (hiv_missing == 0){
      recipe$obs_Xd[,1] = 1
    }
    if (hiv_missing == 1){
      recipe$obs_Xd[recipe$Td[,7]==0,1] = 1
      # recipe$obs_Xd[,1] = ifelse((recipe$Td[,7]==0)&(recipe$obs_Xd[,1]==0), 1, recipe$obs_Xd[,1])
    }
    if (hiv_missing == 2){
      recipe$Td[,7] = 1
    }
  
    
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
    writeLines('')
    cli::cli_inform('{.strong LCA Model Sampler}')
    cli::cli_alert_info('Program starts at: {Sys.time()}')
    cli::cli_h1('Configurations:')
    cli::cli_ul()
    cli::cli_li('{.strong Model} {model}')
    cli::cli_li('Mode: {mode}')
    cli::cli_li('{fold} fold{?s} with {rep} repetition{?s}')
    cli::cli_li('{.strong Prior family:} {.field {prior_family_name}} with {.strong scales} = [{toString(ifelse(penalty_term==0, "~N(0,2.5)", penalty_term))}]')
    cli::cli_li('{.strong Random seed:} {.field {seed}}')
    cli::cli_li('Stan configurations:')
    ulid <- cli::cli_ul()
    cli::cli_li('cores = {.field {cores}}')
    cli::cli_li('chains = {.field {chains}}')
    cli::cli_li('iter = {.field {iter}}')
    cli::cli_li('warmup = {.field {warmup}}')
    cli::cli_li('adapt_delta = {.field {adapt_delta}}')
    cli::cli_li('init_r = {.field {init_r}}')
    cli::cli_end(ulid)
    cli::cli_end()
    cli::cli_h1('')

    if (clean_state){
      if (has_cache) {
        cli::cli_alert("Remove cache and create new folds")
        file.remove(cache_file)
      }
      folds <- create_folds(recipe, fold, rep, seed, cache_file, n_FA, B, lifted_spc, quad_RE=m3_quadRE)
    } else{
      if (has_cache) {
        cli::cli_alert("Cache file found. Use cache file.")
        folds <- readRDS(cache_file)
      } else {
        cli::cli_alert("No cache file found. Create new folds")
        folds <- create_folds(recipe, fold, rep, seed, cache_file, n_FA, B, lifted_spc, quad_RE = m3_quadRE)
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
    if (is.na(sampler)) sampler <- file.path("stan", paste0(model, ".stan"))
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
    pars <- c("z_Smear", "z_Mgit", "z_Xpert", "z_theta",
              "log_lik", "p_Smear", "p_Mgit", "p_Xpert", "theta", "pairwise_corr")
    if (model_no > 0) pars <- c(pars, "a0", "a")
    # if (model != "m0kf" && any(penalty_term == 0)) pars <- c(pars, paste0("sp", c(1,2)[penalty_term==0]))
    if (model_no > 0 && any(penalty_term == 0)) pars <- c(pars, 'sp')
    if (model_no > 0 && all_params) pars <- c(pars, 
      "HIV_a0", "HIV_a", 
      "cs_a0", "cs_a", "L_Omega_cs",
      "mp_a0", "mp_a", "L_Omega_mp",
      "age_a0", "age_a", "age_sigma",
      "id_a0", "id_a", "id_sigma",
      if(model_no != 6) c( "L_Omega_csf", "L_sigma_csf") else c( 'mu_psi_csf', 'sigma_psi_csf','mu_lt_csf', 'sigma_lt_csf', 'psi0_csf','Q_csf'))
    if (!model_no %in% c(0,1)) pars <- c(pars, "b_RE", "b_HIV")
    if (model_no > 1) pars <- c(pars, "b")
    if (model_no == 4) pars <- c(pars, "b_FE")
    if (model == 'm3d') pars <- c(pars, 'a2')
    if (model == 'm3e') pars <- c(pars, 'z_NegCrypto')
    if (model == 'm4d') pars <- c(pars, "a2")
    if (is_pca) pars <- c(pars, "L_csf")
    # if (model_no == 5) pars <- c(pars, "b_cs")
    if (fold == 1 && all_params) pars <- NA
    if (include_d) pars <- c(pars, c('d'))
    if (grepl('missing$', model)) pars <- c(pars, 'z_obs') 
    if (grepl('missingXpert$', model)) pars <- c(pars, 'z_obs_Xpert') 
    if (!all(is.na(include_pars))) pars <- include_pars 
    
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
  })




