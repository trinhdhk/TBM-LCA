# library(rstan)
library(bayesplot)
library(future)
# Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
rstan::rstan_options(auto_write = TRUE)

data_19EI_complete = readRDS('data/impute_19EI.RDS')

model <- rstan::stan_model('stan/model_0.stan')
model <- rstan::stan_model('stan/model_1.stan')
model <- rstan::stan_model('stan/model_2.stan')
model <- rstan::stan_model('stan/model_2.2.stan')
model <- rstan::stan_model('stan/model_3.stan')
model <- rstan::stan_model('stan/model_4.stan')

plan(multisession(workers = 19, gc=TRUE))
fits <- vector("list", length(data_19EI_complete))
for (i in seq_along(data_19EI_complete)){
  fits[[i]] <- future::future({
    options(mc.cores = 1)
    Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
    rstan_options(auto_write = TRUE)
    dat <- data_19EI_complete[[i]]
    
    X <- dat %$% as.matrix(cbind(
      hiv_stat,                             #1
      clin_symptoms,                        #2
      clin_motor_palsy,                     #3
      clin_nerve_palsy,                     #4
      xray_pul_tb,                          #5
      xray_miliary_tb                       #6
      (age-mean(age))/10,                   #7
      log_illness_day,                      #8
      clin_gcs,                             #9
      log2(glucose_ratio),                  #10
      log2(BLDGLU),                         #11
      log_lympho,                           #12
      log_protein,                          #13
      log_lactate,                          #14
      # ISDIABETE                           #15
    ))
    
    model_input_disc <- dat %$%
      list(
        N = nrow(dat), nX = ncol(X),
        Y_Smear = as.integer(csf_smear),
        Y_Mgit = as.integer(csf_mgit),
        Y_Xpert = as.integer(csf_xpert),
        X = X
      )
    
    rstan::sampling(model, 
                    data=model_input_disc, chain = 1,
                    # iter=40000, seed=2906+i, warmup=25000, #model1
                    iter=30000, seed=2906+i, warmup=15000,
                    control = list(max_treedepth = 18, adapt_delta=.7),
                    save_warmup = FALSE,
                    include=FALSE, 
                    pars = c("z_Smear_RE", "z_Mgit_RE", "z_Xpert_RE", "bac_load", "theta")) 
  }, packages = c('magrittr', 'rstan'))
  
  cat('Fitting imputated data', i, '\n')
}

fits_val <- lapply(fits, value)
fits_comb <- rstan::sflist2stanfit(fits_val)
bayesplot::mcmc_intervals(fits_comb, pars = dplyr::vars(starts_with('a')))
saveRDS(fits_comb, 'outputs/fitimpute2.rds')

log_lik_ <- loo::extract_log_lik(fits_comb, merge_chains = FALSE)
r_eff_ <- loo::relative_eff(exp(log_lik_))
loo_3h3a <- loo::loo(log_lik_, r_eff = r_eff_)

log_lik_ <- loo::extract_log_lik(mm, merge_chains = FALSE)
r_eff_ <- loo::relative_eff(exp(log_lik_), cores = 18)
loo_3h3a <- loo::loo(log_lik_, r_eff = r_eff_, moment_match = TRUE, cores=18)

qloo <- function(fit){
  log_lik_ <- loo::extract_log_lik(fit, merge_chains = FALSE)
  r_eff_ <- loo::relative_eff(exp(log_lik_))
  loo_3h3a <- loo::loo(log_lik_, r_eff = r_eff_)
  loo_3h3a
}

plan(multisession(workers = 19, gc=TRUE))
loos <- vector("list", length(fits_val))
for (i in seq_along(fits_val)){
  . <- fits_val[[i]]
  loos[[i]] <- future(qloo(.))
  rm(.)
}

loos <- lapply(loos, value)
elpds <- sapply(loos, function(l) exp(l$pointwise[,'elpd_loo']))
mean_elpds <- log(apply(elpds, 1, mean))

ps <- sapply(loos, function(l) exp(l$pointwise[,'p_loo']))
mean_ps <- log(apply(ps, 1, mean))

ll <- vector("list", length(data_19EI_complete))
for (i in seq_along(data_19EI_complete)){
  . <- fits_val[[i]]
  ll[[i]] <- future(loo::extract_log_lik(., merge_chains = FALSE))
} 
ll <- lapply(ll, value)
ll <- lapply(ll, exp)
# ll2 <- lapply(ll, abind::adrop, drop=2)
sumll <- Reduce(`+`, ll)
avgll <- log(sumll/133)
avgreff <- loo::relative_eff(avgll)
avgloo <- loo::loo(avgll, r_eff = avgreff)
