# Load libraries ----
library(magrittr)
library(data.table)
library(dplyr)

misc <- new.env(hash = FALSE)
source("r/include/functions.R", local = misc)
rstan::rstan_options(javascript = FALSE, auto_write = TRUE)

# Load the folds from cache, otherwise create new
m3_folds <- tryCatch({
  cat('Look for cache file.')
  file <- ".cache/folds/m3_folds.RDS"
  if (!file.exists(file)) stop('')
  cat('>> Cache file found!\n') 
  readRDS(file) 
}, 
error = function(e) {
  cat(">> No cache file found. Create new folds. \n")
  # Load data ----
  recipe <- new.env()
  load('data/cleaned/data_input.Rdata', envir = recipe)
  
  # Prepare the input ----
  m3_input <-
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
           obs_Tc_all = obs_Tc
         )
    )
  
    # Create the folds ----
  m3_folds <- misc$repeated_kfold(m3_input, K = 10, N_rep = 2, N_obs = nrow(recipe$data_19EI), seed = 612)
  saveRDS(m3_folds, file = '.cache/folds/m3_folds.RDS')
  m3_folds
})

m3_inputs <- m3_folds$inputs
m3_outdir <- ".cache/sampling/m3kf"
if (dir.exists(m3_outdir)) unlink(m3_outdir, recursive = TRUE)
dir.create(m3_outdir, showWarnings = F)

# Compile the sampler ----
cat('Compile the sampler\n')
m3_sampler <- rstan::stan_model("stan/m3kf.stan")
cat('Sampling\n')
m3_outputs <- misc$stan_kfold(sampler = m3_sampler,
                              list_of_datas=m3_inputs,
                              backend = "rstan",
                              chains=3, cores=19, 
                              thin = 1, 
                              merge = TRUE,
                              control = list(adapt_delta=.75, max_treedepth=12), 
                              init_r = 1, seed=2693,
                              sample = m3_outdir,
                              pars = c("a0", "a", "b_RE", "b_HIV",
                                       "z_Smear", "z_Mgit", "z_Xpert",
                                       "log_lik", "p_Smear", "p_Mgit", "p_Xpert", "theta"),
                              iter=4000, warmup=2000)
saveRDS(m3_outputs, "outputs/m3kf.RDS")
unlink(m3_outdir, recursive = TRUE)
