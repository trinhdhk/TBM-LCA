# Load libraries ----
library(magrittr)
library(data.table)
library(dplyr)

misc <- new.env(hash = FALSE)
source("r/include/functions.R", local = misc)
rstan::rstan_options(javascript = FALSE, autowrite = TRUE)

# Load the folds from cache, otherwise create new
m0_folds <- tryCatch({
  cat('Look for cache file.')
  if (!file.exists(".cache/folds/m0_folds.RDS")) stop('')
  readRDS(".cache/folds/m0_folds.RDS") 
  cat('>> Cache file found!\n')  
}, 
  error = function(e) {
    cat(">> No cache file found. Create new folds. \n")
    # Load data ----
    recipe <- new.env()
    load('data/cleaned/data_input.Rdata', envir = recipe) 
    
       # Prepare the input ----
    m0_input <-
      with(recipe,
           list(
             N_all = nrow(data_19EI),
             Y_Smear_all = data_19EI$csf_smear,
             Y_Mgit_all = data_19EI$csf_mgit,
             Y_Xpert_all = data_19EI$csf_xpert
           )
      )
    
    # Create the folds ----
    m0_folds <- misc$repeated_kfold(m0_input, K = 10, N_rep = 2, N_obs = nrow(recipe$data_19EI), seed = 612)
    saveRDS(m0_folds, file = '.cache/folds/m0_folds.RDS')
    m0_folds
  })

m0_inputs <- m0_folds$inputs
m0_outdir <- ".cache/sampling/m0kf"
if (dir.exists(m0_outdir)) unlink(m0_outdir, recursive = TRUE)
dir.create(m0_outdir, showWarnings = F)

# Compile the sampler ----
cat('Compile the sampler\n')
m0_sampler <- rstan::stan_model("stan/m0kf.stan")
cat('Sampling\n')
m0_outputs <- misc$stan_kfold(sampler = m0_sampler,
                              list_of_datas=m0_inputs,
                              backend = "rstan",
                              chains=3, cores=19, 
                              thin = 1, 
                              merge = TRUE,
                              control = list(adapt_delta=.75, max_treedepth=12), 
                              init_r = 1, seed=2693,
                              sample = m0_outdir,
                              pars = c("z_Smear", "z_Mgit", "z_Xpert",
                                       "log_lik", "p_Smear", "p_Mgit", "p_Xpert", "theta"),
                              iter=2000, warmup=1000)
saveRDS(m0_outputs, "outputs/m0kf.RDS")
unlink(m0_outdir, recursive = TRUE)
writeLines("Done!")
