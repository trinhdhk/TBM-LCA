library(rstan)
library(bayesplot)
library(future)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
rstan_options(auto_write = TRUE)

data_19EI_complete = readRDS('data/impute_19EI.RDS')[1:126]

model <- rstan::stan_model('stan/model_disc_3h3a_19EI.stan')
model <- rstan::stan_model('stan/model_disc_4_19EI.stan')
model <- rstan::stan_model('stan/model_5_19EI.stan')

plan(multisession(workers = 19, gc=TRUE))
fits <- vector("list", length(data_19EI_complete))
# warn=c()
for (i in seq_along(data_19EI_complete)){
  fits[[i]] <- future::future({
    options(mc.cores = 1)
    Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
    rstan_options(auto_write = TRUE)
    dat <- data_19EI_complete[[i]]
    
    X <- dat %$% as.matrix(cbind(
      hiv_stat,                             #1
      (age-mean(age))/10,                   #2
      clin_symptoms,                        #3
      # log2(clin_illness_day), 
      log_illness_day,                      #4
      clin_nerve_palsy,                     #5
      clin_motor_palsy,                     #6
      ISDIABETE,                            #7
      clin_gcs,                             #8
      glucose_ratio,                        #9
      BLDGLU,                               #10
      # log2(csf_lympho+1),
      # log2(csf_protein),
      # log2(csf_lactate),
      log_lympho,                           #11
      log_protein,                          #12
      log_lactate,                          #13
      xray_pul_tb,                          #14
      xray_miliary_tb                       #15
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
                    iter=37000, seed=2906+i, warmup=25000,
                    control = list(max_treedepth = 18, adapt_delta=.72),
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
loos <- vector("list", length(data_19EI_complete))
for (i in seq_along(data_19EI_complete)){
  . <- fits_val[[i]]
  loos[[i]] <- future(qloo(.))
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
