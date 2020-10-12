library(data.table)
library(magrittr)

data_19EI <- readRDS('data/cleaned/data_19EI.RDS')
data_19EI <- na.omit(data_19EI, cols = c('csf_smear', 'csf_mgit', 'csf_xpert'))
miss_summary <- dplyr::summarise_all(data_19EI, ~sum(is.na(.x)))



data_19EI_compl <- tidyr::replace_na(data_19EI, 
                               replace=as.list(structure(rep(-1, ncol(data_19EI)), names = names(data_19EI))))

#####################
# test.impute = rstan::stan_model('stan/test.impute.stan')
# 
# impute.input = list(
#   N = nrow(data_19EI),
#   HIV = data_19EI[, fcase(!is.na(hiv_stat), as.numeric(hiv_stat), default=0)],
#   cs = as.data.frame(data_19EI[, .(fcase(!is.na(ISCOUGH), as.numeric(ISCOUGH), default=0), 
#                      fcase(!is.na(ISNSWEAT),as.numeric(ISNSWEAT), default=0), 
#                      fcase(!is.na(ISWEIGHT),as.numeric(ISWEIGHT), default=0))]),
#   cbm_cs = fcase(!is.na(data_19EI$clin_symptoms),as.numeric(data_19EI$clin_symptoms),default=0),
#   obs_HIV = !is.na(data_19EI$hiv_stat),
#   obs_cs = as.data.frame(data_19EI[, .(as.numeric(!is.na(ISNSWEAT)), as.numeric(!is.na(ISCOUGH)), as.numeric(!is.na(ISWEIGHT)))]),
#   obs_cbm_cs = !is.na(data_19EI$clin_symptoms)
# )
# sample.impute = rstan::sampling(test.impute, data=impute.input, chains=4, cores=4, iter=10000, seed=1208)
# extract.impute = rstan::extract(sample.impute)
#####################

library(rstan)
options(mc.cores = 6)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
rstan_options(auto_write = TRUE)

Xd <- data_19EI %$% cbind(
  hiv_stat,                             #1
  clin_symptoms,                        #2
  clin_motor_palsy,                     #3
  clin_nerve_palsy,                     #4
  xray_pul_tb,                          #5
  xray_miliary_tb                       #6
)

Xc <- data_19EI %$% cbind(
  (age-mean(age, na.rm=TRUE))/10,       #1    #7
  log2(clin_illness_day),               #2    #8
  15-clin_gcs,                          #3    #9
  log2(BLDGLU),                         #4    #10 
  log2(csf_glucose),                    #5    #11
  log2(csf_lympho+1),                   #6    #12
  log2(csf_protein),                    #7    #13
  log2(csf_lactate)                     #8    #14
)

Td <- data_19EI %$% cbind(
  ISWEIGHT,                             #1
  ISNSWEAT,                             #2
  ISCOUGH,                              #3
  HEMIPLEGIA,                           #4
  PARAPLEGIA,                           #5
  TETRAPLEGIA,                          #6
  ISDIABETE                             #7
)

Tc <- data_19EI %$% cbind(
  4-GCSE,                               #1
  6-GCSM,                               #2
  5-GCSV                                #3
)

is.not.na <- Negate(is.na)
obs_Xd <- apply(Xd, 2, is.not.na) %>% apply(2, as.numeric)
obs_Xc <- apply(Xc, 2, is.not.na) %>% apply(2, as.numeric)
obs_Td <- apply(Td, 2, is.not.na) %>% apply(2, as.numeric)
obs_Tc <- apply(Tc, 2, is.not.na) %>% apply(2, as.numeric)

my.replace_na <- function(data, replace=0){
  apply(data, 2, tidyr::replace_na, replace=replace)
}
Xd <- my.replace_na(Xd)
Xc <- my.replace_na(Xc)
Td <- my.replace_na(Td)
Tc <- my.replace_na(Tc)

model_input <-
  list(
    N = nrow(data_19EI),
    nXc = ncol(Xc),
    nXd = ncol(Xd),
    nTc = ncol(Tc),
    nTd = ncol(Td),
    Y_Smear = data_19EI$csf_smear,
    Y_Mgit = data_19EI$csf_mgit,
    Y_Xpert = data_19EI$csf_xpert,
    Xc = Xc,
    Xd = Xd,
    Tc = Tc,
    Td = Td,
    obs_Xc = obs_Xc,
    obs_Xd = obs_Xd,
    obs_Tc = obs_Tc,
    obs_Td = obs_Td
  )
Sys.setenv(STAN_NUM_THREADS = 4)
model_2is <- rstan::stan_model("stan/model_2i.stan")
model_2i <- cmdstanr::cmdstan_model('stan/model_2i.stan', dir = 'stan', cpp_options = list(stan_threads = TRUE))
# model_2i <- cmdstanr::cmdstan_model('stan/model_2i.stan', dir = 'stan')
res_2i   <- model_2i$sample(data = model_input, seed = 29061993, output_dir = 'outputs',
                            chains = 3, parallel_chains = 3, threads_per_chain = 4,
                            iter_sampling = 10000, iter_warmup = 15000,
                            show_messages = FALSE,
                            adapt_delta = .8, max_treedepth = 18)
draw2i = res_2i$draws()
draw_2i = res_2i$draws(c("a0", "a", "b_HIV", "b", "z_Smear", "z_Mgit", "z_Xpert"))
bayesplot::mcmc_intervals(draw_2i, pars=vars(dplyr::starts_with("a")))
# 
# model_2iv <- cmdstanr::cmdstan_model('stan/model_2i_v.stan', dir = 'stan', cpp_options = list(stan_threads = TRUE))
# res_2iv   <- model_2iv$sample(data = model_input, seed = 1208, output_dir = 'outputs',
#                             chains = 3, parallel_chains = 3, threads_per_chain = 4, 
#                             iter_sampling = 15000, iter_warmup = 15000, 
#                             # show_messages = FALSE,
#                             adapt_delta = .8, max_treedepth = 18)
# 
# for (i in names(model_input)) assign(i, model_input[[i]])
# rstan::stan_rdump(names(model_input), file='in.R')
# model_2i$save_hpp_file(dir='stan')
