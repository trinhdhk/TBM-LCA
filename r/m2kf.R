

misc <- new.env(hash = FALSE)
source("r/include/functions.R", local = misc)
rstan::rstan_options(javascript = FALSE, auto_write = TRUE)

# Load the folds from cache, otherwise create new
m2_folds <- tryCatch({
  cat('Look for cache file.')
  file <- ".cache/folds/m2_folds.RDS"
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
  m2_input <-
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
  m2_folds <- misc$repeated_kfold(m2_input, K = 10, N_rep = 4, N_obs = nrow(recipe$data_19EI), seed = 612)
  saveRDS(m2_folds, file = '.cache/folds/m2_folds.RDS')
  m2_folds
})

m2_inputs <- m2_folds$inputs
m2_outdir <- ".cache/sampling/m2kf"
if (dir.exists(m2_outdir)) unlink(m2_outdir, recursive = TRUE)
dir.create(m2_outdir, showWarnings = F)

# Compile the sampler ----
cat('Compile the sampler\n')
m2_sampler <- rstan::stan_model("stan/m2kf.stan")
cat('Sampling\n')
m2_outputs <- misc$stan_kfold(sampler = m2_sampler,
                              list_of_datas=m2_inputs,
                              backend = "rstan",
                              chains=3, cores=19, 
                              thin = 2, 
                              merge = TRUE,
                              control = list(adapt_delta=.75, max_treedepth=12), 
                              init_r = 1, seed=2693,
                              sample = m2_outdir,
                              pars = c("a0", "a", "b_RE", "b_HIV",
                                       "z_Smear", "z_Mgit", "z_Xpert",
                                       "log_lik", "p_Smear", "p_Mgit", "p_Xpert", "theta"),
                              iter=3000, warmup=2000)
unlink(m2_outdir, recursive = TRUE)
saveRDS(m2_outputs, "outputs/m2kf.RDS")
