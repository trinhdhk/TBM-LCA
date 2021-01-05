library(magrittr)
library(data.table)
library(dplyr)
source("r/include/functions.R")
data_19EI <- readRDS('data/cleaned/data_19EI.RDS')
data_19EI <- na.omit(data_19EI, cols = c('BLDGLU', 'csf_glucose', 'csf_lympho', 'csf_protein', 'csf_lactate'))

latentvar_input <- list(
  N = nrow(data_19EI),
  M = 6,
  K = 3,
  X = data_19EI %$% cbind(
    glu=sqrt(BLDGLU),                          #3    #9
    csfglu=sqrt(csf_glucose),
    gluratio = sqrt(glucose_ratio),
    csflym=log2(csf_lympho+1),                   #5    #11   
    csfpro=log2(csf_protein),                    #6    #12   
    csflac=log2(csf_lactate)
  )                    #7    #13  
)


rstan::rstan_options(auto_write=TRUE)
latentvar_sampler <- rstan::stan_model("stan/latentvar.stan")

latentvar_result <- rstan::sampling(latentvar_sampler, data=latentvar_input,chains=4,
                                    cores = 8, control=list(max_treedepth=15, adapt_delta=.82), iter=6000, warmup=1500)

rstan::loo(latentvar_result, moment_match=F) %>% plot
bayesplot::mcmc_hist_by_chain(latentvar_result, regex_pars = "^Y\\[1,")
rstan::rstan_options()