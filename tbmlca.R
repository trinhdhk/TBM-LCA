#!/usr/bin/env Rscript

library(argparser)
library(magrittr, include.only = "%>%")
misc <- new.env(hash = FALSE)
source("r/include/functions.R", local = misc)
rstan::rstan_options(javascript = FALSE, auto_write = TRUE)

args <- 
  arg_parser(paste("TBM-LCA model runner - Written by Trinh Dong, build", Sys.Date())) %>%
  add_argument("model", help = "model to run") %>%
  # add_argument("--help", help = "print this help page and exit", short = "-h", flag = TRUE)
  add_argument("--input", help = "input Rdata file", 
               default = "data/cleaned/data_input.Rdata", short = "-i") %>%
  add_argument("--sampler", help = "sampler file, default to stan/<model>.stan", short = "-s") %>% 
  add_argument("--output-dir", help = "output folder", default = "outputs", short = "-o") %>%
  add_argument("--output-file", help = "output filename, default to auto-inference", short = "-f") %>%
  add_argument("--fold", help = "number of folds for cross-validation, ignored if not in a clean state and cache is available", default = 10, short = "-k") %>%
  add_argument("--rep", help = "number of repetitions of K-fold procedure, ignored if not in a clean state and cache is available", default = 4, short = "-r") %>%
  add_argument("--seed", help = "random seed", default = 9802, short = "-S") %>%
  add_argument("--clean-state", help = "ignore cache to run a clean state", flag = TRUE) %>% 
  add_argument("--cache-dir", help = "cache directory, ignored if --no-cache flag is on", default = ".cache", short = "-d") %>%
  add_argument("--no-cache", help = "disable caching and monitoring process, written to <cache-dir>/sampling/<model> ", flag = TRUE, short = "-n") %>% 
  add_argument("--no-clear-cache", help = "disable clearing cache after running", flag = TRUE, short = "-l") %>%
  add_argument("--chains", help = "number of chains per dataset", default = 3, short = "-N") %>%
  add_argument("--cores", help = "number of cores", default = parallel::detectCores()-1, short = "-C") %>%
  add_argument("--thin", help = "thining the chains", default = 1, short = "-T") %>%
  add_argument("--iter", help = "number of iterations", default = 5000, short = "-I") %>%
  add_argument("--warmup", help = "number of warm-up iterations", default = 2000, short = "-W") %>% 
  add_argument("--init-r", help = "Stan init parameter", default = 1, short = "-R") %>%
  add_argument("--penalty-family", help = "Family for the regularised prior, either 't', 'normal', or 'laplace'", default = "t", short = "-F") %>%
  add_argument("--penalty-term", help = "Penalty scale term, if 0 then auto adaptation is performed", default = 0, short = "-P")

argparser <- parse_args(args)

# Functions to create folds
create_folds <- function(recipe, K, N, seed, cache_file=NULL){
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
           Xc_all = Xc,
           Xd_all = Xd,
           Td_all = Td,
           Tc_all = Tc,
           obs_Xc_all = obs_Xc,
           obs_Xd_all = obs_Xd,
           obs_Td_all = obs_Td,
           obs_Tc_all = obs_Tc,
           penalty_family = penalty_family,
           penalty_term = penalty_term
         )
    )
  
  folds <- misc$repeated_kfold(inp, K = K, N_rep = N, N_obs = nrow(recipe$data_19EI), seed = seed)
  if (!is.null(cache_file)) saveRDS(folds, file = cache_file)
  folds
}

with(
  argparser,
  {
    # If no name for output file, do auto-inference
    if (is.na(output_file)) output_file <- paste(model,"RDS",sep=".")
    
    # Coerce numerics
    fold <- as.integer(fold)
    rep  <- as.integer(rep)
    chains <- as.integer(chains)
    cores  <- as.integer(cores)
    seed   <- as.integer(seed)
    iter   <- as.integer(iter)
    warmup <- as.integer(warmup)
    init_r <- as.integer(init_r)
    penalty_term <- as.integer(penalty_term)
    penalty_family <- switch(penalty_family, 
      "t" = 0,
      "laplace" = 1,
      "normal" = 2)
    if (is.na(penalty_family)) stop("Invalid penalty family.")
    
    # Load the recipe
    recipe <- new.env()
    load(input, envir = recipe)
    recipe$penalty_term = penalty_term
    recipe$penalty_family = penalty_family
    
    # Create results env
    results <- new.env()
    
    # If there are cache to use, create corresponding dirs
    if (!no_cache) {
      dir.create(file.path(cache_dir, "sampling"), recursive = TRUE, showWarnings = FALSE, mode = "770")
      dir.create(file.path(cache_dir, "folds"), recursive = TRUE, showWarnings = FALSE, mode = "770")
    }
    
    # Load cache or create new folds
    cache_file <- file.path(cache_dir, "folds", paste0(model, ".RDS"))
    has_cache <- file.exists(cache_file)
    if (clean_state){
      if (has_cache) {
        writeLines(">> Remove cache and create new folds")
        file.remove(cache_file)
      }
      folds <- create_folds(recipe, fold, rep, seed, cache_file)
    } else{
      if (has_cache) {
        writeLines(">> Cache file found. Use cache file.")
        folds <- readRDS(cache_file)
      } else {
        writeLines(">> No cache file found. Create new folds")
        folds <- create_folds(recipe, fold, rep, seed, cache_file)
      }
    }
    
    # Load the folds
    inputs <- folds$inputs
    
    # Remove the cache dir if exists
    if (!no_cache) {
      outdir <- file.path(cache_dir, "sampling", model)
      if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
      dir.create(outdir, showWarnings = F)
    } else outdir <- NULL
    
    cat('>> Compile the sampler\n')
    if (is.na(sampler)) sampler <- file.path("stan", paste0(model, ".stan"))
    sampler <- tryCatch(
      rstan::stan_model(sampler),
      error = function(e) {
        # Remove the cache file might fix the problem
        file.remove(file.path("stan", paste0(model, ".rds")))
        rstan::stan_model(sampler)
      }
    )
    
    writeLines(">> Sample")
    pars <- c("z_Smear", "z_Mgit", "z_Xpert",
              "log_lik", "p_Smear", "p_Mgit", "p_Xpert", "theta")
    if (model != "m0kf") pars <- c(pars, "a0", "a")
    if (model %in% c("m2kf", "m3kf")) pars <- c(pars, "b_RE", "b_HIV")
    if (model != "m0kf" && penalty_term == 0) pars <- c(pars, "sp")
    results$outputs <- misc$stan_kfold(sampler = sampler,
                               list_of_datas=inputs,
                               backend = "rstan",
                               chains = chains, cores = cores, 
                               thin = thin, 
                               merge = TRUE,
                               control = list(adapt_delta=.75, max_treedepth=12), 
                               init_r = init_r, seed = seed,
                               sample_dir = outdir,
                               pars = pars,
                               iter=iter, warmup=warmup)
    
    # Clear cache
    if (!no_clear_cache){
      writeLines(">> Clean up cache")
      unlink(outdir, recursive = TRUE)
      file.remove(cache_file)
    } 
    
    # Save results
    writeLines(">> Save results")
    results$folds  <- folds
    results$seed   <- seed
    results$fold <- fold
    results$rep  <- rep
    results$chains <- chains
    results$iter   <- iter
    results$warmup <- warmup
    results$init_r <- init_r
    saveRDS(results, file = file.path(output_dir, output_file))
    future::plan("sequential")
  })




