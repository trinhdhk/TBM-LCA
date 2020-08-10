rm(list = ls())
library(data.table)
library(magrittr)
library(rstan)
library(bayesplot)
library(future)
# options(mc.cores = 18)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
rstan_options(auto_write = TRUE)
source('r/fn_def.R')

data_19EI <- readRDS('data/cleaned/data_19EI.RDS')
# data_19EI$age <-2020+data_19EI$age 
data_19EI <- na.omit(data_19EI, cols = c('csf_smear', 'csf_mgit', 'csf_xpert'))
data_19EI <- data_19EI[, c('hiv_stat', 'age', 'csf_smear', 'csf_mgit', 'csf_xpert',
                           'clin_illness_day', 'clin_gcs', 
                           'clin_nerve_palsy',
                           # 'clin_motor_palsy', 
                           'csf_clear', 'csf_protein', 'csf_lympho', 'glucose_ratio', 'csf_lactate',
                           'xray_miliary_tb', 'xray_pul_tb', 'ISDIABETE', 'BLDGLU', 'ISNSWEAT', 'ISCOUGH', 'ISWEIGHT',
                           'HEMIPLEGIA', 'PARAPLEGIA', 'TETRAPLEGIA'), with=F]
miss <- dplyr::summarise_all(data_19EI, ~sum(is.na(.x)))

data_19EI[, `:=`(log_illness_day = log2(clin_illness_day), clin_illness_day=NULL,
                 log_lympho = log2(csf_lympho+1),csf_lympho=NULL,
                 log_protein = log2(csf_protein),csf_protein=NULL,
                 log_lactate = log2(csf_lactate),csf_lactate=NULL)]

pred.mat <- mice::make.predictorMatrix(data_19EI)
pred.mat[3:5,] <- 0
pred.mat[14:19,] <- 0
pred.mat[14:16,14:16] <- 1
pred.mat[14:16,1] <- 1
pred.mat[17:19,] <- 0
pred.mat[17:19,c(1,12:13,17:19)] <- 1
pred.mat[17:19,1] <- 1
pred.mat[, 3:5] <- 1
data_19EI_imp <- mice::parlmice(data_19EI, n.core = 18, maxit=50, m=132, n.imp.core = 8, predictorMatrix = pred.mat, defaultMethod = c("pmm", "logreg", "polyreg", "polr"))
data_19EI_complete <- lapply(seq_len(data_19EI_imp$m), mice::complete, data = data_19EI_imp)
data_19EI_complete <- lapply(data_19EI_complete, function(dt)
  dplyr::mutate(dt,clin_symptoms = ISNSWEAT|ISCOUGH|ISWEIGHT, clin_motor_palsy = HEMIPLEGIA|PARAPLEGIA|TETRAPLEGIA))
model <- rstan::stan_model('stan/model_disc_3h3a_19EI.stan')
plan(multisession(workers = 18, gc=TRUE))

# model <- cmdstanr::cmdstan_model('stan/model_disc_3h3a_19EI.stan')

fits <- vector("list", 132) #length(data_19EI_complete))
# warn=c()
for (i in 1:132){#seq_along(data_19EI_complete)
  fits[[i]] <- future::future({
    options(mc.cores = 1)
    library(magrittr)
    library(rstan)
    # Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
    # dir.create(file.path('outputs', Sys.Date()))
    dat <- data_19EI_complete[[i]]
    
    X <- dat %$% as.matrix(cbind(
      hiv_stat,
      (age-mean(age))/10, 
      clin_symptoms,
      # log2(clin_illness_day), 
      log_illness_day,
      clin_nerve_palsy,
      clin_motor_palsy,
      ISDIABETE,
      clin_gcs,
      glucose_ratio,
      BLDGLU,
      # log2(csf_lympho+1),
      # log2(csf_protein),
      # log2(csf_lactate),
      log_lympho,
      log_protein,
      log_lactate,
      xray_pul_tb,
      xray_miliary_tb
    ))
    
    model_input_disc <- dat %$%
      list(
        N = nrow(dat), nX = ncol(X),
        Y_Smear = as.integer(csf_smear),
        Y_Mgit = as.integer(csf_mgit),
        Y_Xpert = as.integer(csf_xpert),
        X = X
      )
    
    # res <- model$sample(data=model_input_disc, chains = 3, seed = 2906, iter_warmup = 20000,
    #                     iter_sampling = 15000, max_treedepth = 18, adapt_delta=.8, output_dir = file.path('outputs', Sys.Date()))
    # cmdstanr::read_cmdstan_csv(res$output_files(), variables = c('a0', 'a', 'b_HIV', 'b_age', 'b_mil', 'b',
    #                                                              'z_Smear', 'z_Mgit', 'z_Xpert', 'RE', 'log_lik'))
    sampling(model, 
             data=model_input_disc, chain = 1,
             iter=21000, seed=2906+i, warmup=15000, 
             control = list(max_treedepth = 18, adapt_delta=.75),
             # sample_file = paste0('outputs/',Sys.Date(),'/model_',i), 
             include=FALSE, 
             pars = c("z_Smear_RE", "z_Mgit_RE", "z_Xpert_RE", "bac_load", "theta")) 
  })
  
  cat('Fitting imputated data', i, '\n')
}

fits_val <- lapply(fits, value)
fits_comb <- rstan::sflist2stanfit(fits_val)
bayesplot::mcmc_intervals(fits_comb, pars = dplyr::vars(starts_with('a')))
saveRDS(fits_comb, 'outputs/fitimpute2.rds')

log_lik_ <- loo::extract_log_lik(fits_comb, merge_chains = FALSE)
r_eff_ <- loo::relative_eff(exp(log_lik_), cores = 8)
loo_3h3a <- loo::loo(log_lik_, r_eff = r_eff_, moment_match = TRUE, cores=2)

log_lik_ <- loo::extract_log_lik(mm, merge_chains = FALSE)
r_eff_ <- loo::relative_eff(exp(log_lik_), cores = 18)
loo_3h3a <- loo::loo(log_lik_, r_eff = r_eff_, moment_match = TRUE, cores=18)

