rm(list = ls())
library(data.table)
library(magrittr)
library(rstan)
library(bayesplot)
options(mc.cores = 6)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
rstan_options(auto_write = TRUE)
source('r/fn_def.R')

data_19EI <- readRDS('data/cleaned/data_19EI.RDS')
# data_19EI$age <-2020+data_19EI$age 
data_19EI <- data_19EI[, c('hiv_stat', 'age', 'csf_smear', 'csf_mgit', 'csf_xpert',
                           'clin_illness_day', 'clin_symptoms', 'clin_gcs', 
                           'clin_nerve_palsy', 'clin_motor_palsy', 
                           'csf_clear', 'csf_protein', 'csf_lympho', 'glucose_ratio', 'csf_lactate',
                           'xray_miliary_tb', 'xray_pul_tb'), with=F]
data_19EI <- na.omit(data_19EI, cols = c('hiv_stat', 'age', 'csf_smear', 'csf_mgit', 'csf_xpert',
                                      'clin_illness_day', 'clin_symptoms', 'clin_gcs', 
                                      'clin_nerve_palsy', 'clin_motor_palsy', 
                                      'csf_clear', 'csf_protein', 'csf_lympho', 'glucose_ratio', 'csf_lactate',
                                      'xray_miliary_tb', 'xray_pul_tb'))

X <- data_19EI %$% as.matrix(cbind(
  hiv_stat,
  (age-mean(age))/10, 
  clin_symptoms,
  log2(clin_illness_day), 
  clin_nerve_palsy,
  clin_motor_palsy,
  ISDIABETE,
  clin_gcs,
  glucose_ratio,
  BLDGLU,
  log2(csf_lympho+1),
  log2(csf_protein),
  log2(csf_lactate),
  xray_pul_tb,
  xray_miliary_tb
))

Xp <- data_19EI %$% as.matrix(cbind(
  ISDIABETE
))



model_input_disc <- data_19EI %$%
  list(
    N = nrow(data_19EI), nX = ncol(X),
    Y_Smear = as.integer(csf_smear),
    Y_Mgit = as.integer(csf_mgit),
    Y_Xpert = as.integer(csf_xpert),
    # Y_Img = as.integer(xray_miliary_tb),
    X = X
    # diabetes = data_19EI$ISDIABETE
  )

model_3h3a_19EI <- load_model('disc_3h3a_19EI', data=model_input_disc, chain = 3, iter=20000, seed=2906, warmup=10000, 
                              control = list(max_treedepth = 18, adapt_delta=.8),
                              # init_r=.5,
                              include = FALSE,
                              pars = c("z_Smear_RE", "z_Mgit_RE", "z_Xpert_RE", "bac_load"))

# model_3h3a_19EI.model <- stan_model(file.path('stan', paste0('model_', 'disc_3h3a_19EI', '.stan')))
model <- rstan::stan_model('stan/model_disc_3h3a_19EI.stan')
model_3h3a_19EI <- sampling(model, data=model_input_disc, chain = 3, iter=10000, seed=2906, warmup=5000,
                   control = list(max_treedepth = 18, adapt_delta=.7),
                   include = FALSE,
                   pars = c("z_Smear_RE", "z_Mgit_RE", "z_Xpert_RE", "bac_load"))

e3h3a = rstan::extract(model_3h3a_19EI)

gge = ggmcmc::ggs(model_3h3a_19EI)
ggmcmc::ggs_running(gge, family="z_Xpert")
ggmcmc::ggs_pairs(gge, family="z_Xpert")
bayesplot::mcmc_areas(model_3h3a_19EI, pars=dplyr::vars(starts_with('a')))
bayesplot::mcmc_intervals(model_3h3a_19EI, pars=dplyr::vars(starts_with('a')))

posterior_ <- as.array(model_3h3a_19EI)
lp_ <- log_posterior(model_3h3a_19EI)
np_ <- nuts_params(model_3h3a_19EI)
color_scheme_set("darkgray")
mcmc_parcoord(posterior_, np = np_, pars = vars(-lp__))
mcmc_pairs(posterior_, np=np_, pars=c('b[1]','b_HIV'))


mcmc_intervals(posterior_, pars = vars(dplyr::starts_with('a'))) -> p_a
p_a + scale_y_discrete(labels = c("intercept", "HIV", "age", "sex", "illness day", "history contact tb", "cranial palsy", "motor palsy", "fever", "headache", "GCS", "glucose_ratio", "blood glucose", "csf lympho", "csf protein", "pulmonary tb", "miliary tb") %>% `[`(., length(.):1))

mcmc_intervals(posterior_, pars = vars(dplyr::starts_with("b"))) -> p_b

library(LaplacesDemon)
mcmc_intervals(posterior_, pars = vars(dplyr::starts_with("z")), transformations = "invlogit") -> p_z

log_lik_ <- loo::extract_log_lik(model_3h3a_19EI, merge_chains = FALSE)
r_eff_ <- loo::relative_eff(exp(log_lik_), cores = 3)
loo_3h3a <- loo::loo(log_lik_, r_eff = r_eff_, moment_match = TRUE)

# create a named list of draws for use with rstan methods
.rstan_relist <- function(x, skeleton) {
  out <- utils::relist(x, skeleton)
  for (i in seq_along(skeleton)) {
    dim(out[[i]]) <- dim(skeleton[[i]])
  }
  out
}

# rstan helper function to get dims of parameters right
.create_skeleton <- function(pars, dims) {
  out <- lapply(seq_along(pars), function(i) {
    len_dims <- length(dims[[i]])
    if (len_dims < 1) return(0)
    return(array(0, dim = dims[[i]]))
  })
  names(out) <- pars
  out
}

# extract original posterior draws
post_draws_stanfit <- function(x, ...) {
  as.matrix(x)
}

# compute a matrix of log-likelihood values for the ith observation
# matrix contains information about the number of MCMC chains
log_lik_i_stanfit <- function(x, i, parameter_name = "log_lik", ...) {
  loo::extract_log_lik(x, parameter_name, merge_chains = FALSE)[, , i]
}

# transform parameters to the unconstraint space
unconstrain_pars_stanfit <- function(x, pars, ...) {
  skeleton <- .create_skeleton(x@sim$pars_oi, x@par_dims[x@sim$pars_oi])
  upars <- apply(pars, 1, FUN = function(theta) {
    rstan::unconstrain_pars(x, .rstan_relist(theta, skeleton))
  })
  # for one parameter models
  if (is.null(dim(upars))) {
    dim(upars) <- c(1, length(upars))
  }
  t(upars)
}

# compute log_prob for each posterior draws on the unconstrained space
log_prob_upars_stanfit <- function(x, upars, ...) {
  apply(upars, 1, rstan::log_prob, object = x,
        adjust_transform = TRUE, gradient = FALSE)
}

# compute log_lik values based on the unconstrained parameters
log_lik_i_upars_stanfit <- function(x, upars, i, parameter_name = "log_lik",
                                    ...) {
  S <- nrow(upars)
  out <- numeric(S)
  for (s in seq_len(S)) {
    out[s] <- rstan::constrain_pars(x, upars = upars[s, ])[[parameter_name]][i]
  }
  out
}

lmm <- loo::loo_moment_match(x = model_3h3a_19EI,
                             loo = loo_3h3a,
                             cores = 2,
                             post_draws = post_draws_stanfit,
                             log_lik_i = log_lik_i_stanfit,
                             unconstrain_pars = unconstrain_pars_stanfit,
                             log_prob_upars = log_prob_upars_stanfit,
                             log_lik_i_upars = log_lik_i_upars_stanfit)
