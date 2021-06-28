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
  add_argument("--penalty-term", help = "Penalty scale term, if 0 then auto adaptation is performed", nargs = 2, default = c(0,0), short = "-P") %>%
  add_argument("--all-params", help = "Whether to expose imputation params, ignored for m0", flag = TRUE, short="-a") %>%
  add_argument("--n-FA", help = "For model 6 only: number of latent factor", default = 2, short = "-F")

argparser <- parse_args(args)

# Functions to create folds
create_folds <- function(recipe, K, N, seed, cache_file=NULL, n_FA){
  inp <- 
    with(recipe,
         list(
           N_all = nrow(data_19EI),
           nFA = n_FA,
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
        allow_undefined = allow_undefined, allow_optimizations = allow_optimizations, 
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

with(
  argparser,
  {
    # If no name for output file, do auto-inference
    if (is.na(output_file)) output_file <- paste(model,"RDS",sep=".")
    output_name <- fs::path_ext_remove(fs::path_file(output_file))
    # Coerce numerics
    fold <- as.integer(fold)
    rep  <- as.integer(rep)
    chains <- as.integer(chains)
    cores  <- as.integer(cores)
    seed   <- as.integer(seed)
    iter   <- as.integer(iter)
    warmup <- as.integer(warmup)
    init_r <- as.numeric(init_r)
    penalty_term <- as.integer(penalty_term)
    n_FA <- as.integer(`n-FA`)
    
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
    
    # Create results env
    results <- new.env()
    
    # If there are cache to use, create corresponding dirs
    if (!no_cache) {
      dir.create(file.path(cache_dir, "sampling"), recursive = TRUE, showWarnings = FALSE, mode = "770")
      dir.create(file.path(cache_dir, "folds"), recursive = TRUE, showWarnings = FALSE, mode = "770")
    }
    
    # Load cache or create new folds
    cache_file <- file.path(cache_dir, "folds", paste0(output_name, "_recipe.RDS"))
    has_cache <- file.exists(cache_file)
    writeLines(paste(">> Program starts at:", Sys.time()))
    writeLines(crayon::yellow("-----------------------------------------"))
    writeLines(crayon::yellow("Configuration:"))
    writeLines(crayon::yellow("- Model", model))
    writeLines(crayon::yellow("-", fold, "fold(s) with", rep, "repetition(s)"))
    writeLines(crayon::yellow("- Prior family:", switch(as.character(penalty_family), "0" = "Student t", "1" = "Laplace", "2" = "Normal"),
      "with variances = {", toString(ifelse(penalty_term==0, "~U[0,10]", penalty_term)), "}"))
    writeLines(crayon::yellow("- Seed =", seed))
    writeLines(crayon::yellow("- Stan config: Cores =", cores, "Chains =", chains, "Iter =", iter, "Warmup =", warmup, "Init_R =", init_r))
    if (model == "m6kf")
      writeLines(crayon::yellow("- Number of latent factors:", n_FA))
    writeLines(crayon::yellow("-----------------------------------------"))

    if (clean_state){
      if (has_cache) {
        writeLines(">> Remove cache and create new folds")
        file.remove(cache_file)
      }
      folds <- create_folds(recipe, fold, rep, seed, cache_file, n_FA)
    } else{
      if (has_cache) {
        writeLines(">> Cache file found. Use cache file.")
        folds <- readRDS(cache_file)
      } else {
        writeLines(">> No cache file found. Create new folds")
        folds <- create_folds(recipe, fold, rep, seed, cache_file, n_FA)
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
    
    cat('>> Compile the sampler\n')
    if (is.na(sampler)) sampler <- file.path("stan", paste0(model, ".stan"))
    sampler <- tryCatch(
      my_stan_model(sampler),
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
    if (model != "m0kf" && any(penalty_term == 0)) pars <- c(pars, paste0("sp", c(1,2)[penalty_term==0]))
    if (model != "m0kf" && all_params) pars <- c(pars, 
      "HIV_a0", "HIV_a", 
      "cs_a0", "cs_a", "L_Omega_cs",
      "mp_a0", "mp_a", "L_Omega_mp",
      "age_a0", "age_a", "age_sigma",
      "id_a0", "id_a", "id_sigma",
      if(model != "m6kf") c("csf_a0", "L_Omega_csf", "L_sigma_csf") else c( 'mu_psi_csf', 'sigma_psi_csf','mu_lt_csf', 'sigma_lt_csf', 'psi0_csf','Q_csf'))
    if (!model %in% c("m0kf", "m1kf")) pars <- c(pars, "b_RE", "b_HIV")
    if (model == "m4kf") pars <- c(pars, "b_XpertLevel")
    if (model == "m5kf") pars <- c(pars, "c_Xpert")
    if (model == "m6kf") pars <- c(pars, "U_csf_all", "L_csf")
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
    results$penalty <- list(family = penalty_family, term = penalty_term)
    saveRDS(results, file = file.path(output_dir, output_file))
    future::plan("sequential")
  })




