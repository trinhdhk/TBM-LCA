library(data.table)
library(magrittr)

data_19EI <- readRDS('data/cleaned/data_19EI.RDS')
data_19EI <- na.omit(data_19EI, cols = c('csf_smear', 'csf_mgit', 'csf_xpert'))
data_19EI[, clin_reducedconsciousness := (ISREDUCED %in% TRUE) | (clin_gcs < 15 %in% TRUE)]
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
  (15-clin_gcs),                        #3    #9
  sqrt(BLDGLU),                         #4    #10 
  sqrt(csf_glucose),                    #5    #11
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
moedl_2i_extI2 <- rstan::stan_model('stan/model_2i_exti.stan')
model_2iv <- cmdstanr::cmdstan_model('stan/model_2iv.stan', dir = 'stan', cpp_options = list(stan_threads = TRUE))
model_2i_extI <- cmdstanr::cmdstan_model('stan/model_2i_exti.stan', dir = 'stan')
model_0 <- cmdstanr::cmdstan_model('stan/model_0.stan', dir = 'stan')
model_2i_extIg <- cmdstanr::cmdstan_model('stan/model_2i_extig.stan', dir = 'stan')
model_2ih <- cmdstanr::cmdstan_model('stan/model_2ih.stan', dir = 'stan', cpp_options = list(stan_threads = TRUE))
# model_2i <- cmdstanr::cmdstan_model('stan/model_2i.stan', dir = 'stan')
res_2iv   <- model_2iv$sample(data = model_input, seed = 29061993, output_dir = 'outputs',
                            chains = 4, parallel_chains = 4, threads_per_chain = 4,
                            iter_sampling = 10000, iter_warmup = 15000,
                            show_messages = FALSE,
                            adapt_delta = .82, max_treedepth = 18)

res_2ih   <- model_2ih$sample(data = model_input, seed = 29061993, output_dir = 'outputs',
                              chains = 4, parallel_chains = 4, threads_per_chain = 4,
                              iter_sampling = 10000, iter_warmup = 15000,
                              show_messages = FALSE,
                              adapt_delta = .9, max_treedepth = 18)

res_2i_extI   <- model_2i_extI$sample(data = model_input, seed = 29061993, output_dir = 'outputs',
                                      chains = 4, parallel_chains = 4,
                                      iter_sampling = 1500, iter_warmup = 2000,
                                      show_messages = FALSE,
                                      adapt_delta = .88, max_treedepth = 18)
res_2i_extIg   <- model_2i_extIg$sample(data = model_input, seed = 29061993, output_dir = 'outputs',
                                        chains = 4, parallel_chains = 4,
                                        iter_sampling = 2000, iter_warmup = 2000,
                                        show_messages =FALSE,save_warmup = F,
                                        adapt_delta = .88, max_treedepth = 18)

res_0   <- model_0$sample(data = model_input, seed = 29061993, output_dir = 'outputs',
                                        chains = 4, parallel_chains = 4,
                                        iter_sampling = 5000, iter_warmup = 5000,
                                        show_messages = FALSE,
                                        adapt_delta = .8, max_treedepth = 18)
# draw2i = res_2i$draws()
draw_0=res_0$draws()
draw_2ih = res_2ih$draws(c("a0", "a", "b_HIV", "b", "z_Smear", "z_Mgit", "z_Xpert"))
draw_2i_extI = res_2i_extI$draws(c("a0", "a", "b_HIV", "b_RE", "z_Smear", "z_Mgit", "z_Xpert", "HIV_a0"))
draw_2i_extIg = res_2i_extIg$draws(c("a0", "a", "b_HIV", "b_RE", "z_Smear", "z_Mgit", "z_Xpert", "HIV_a0", "theta", "C"))
bayesplot::mcmc_intervals(draw_2i_extI, pars=vars(dplyr::starts_with("a")))
bayesplot::mcmc_intervals(draw_2i_extIg, pars=vars(dplyr::starts_with("a")))
bayesplot::mcmc_intervals(draw_2i_extI, pars=vars(dplyr::starts_with("z")),  prob_outer = .95, transformations = nimble::ilogit)
bayesplot::mcmc_intervals(draw_2i_extIg, pars=vars(dplyr::starts_with("z")), prob_outer = .95, transformations = nimble::ilogit)
saveRDS(res_2i, 'outputs/model_2iv1.RDS')

ylabs=rev(c("Intercept", "HIV", "TB Symptoms", "Motor Palsy", "Cranial Nerve Palsy",
        "Pulmonary TB", "Miliary TB",
        "Age/10", "Log(Illness days)", "GCS Loss (15-GCS)",
        "Sqrt(Blood Glucose)", "Sqrt(CSF Glucose)", 
        "Log(CSF Lympho)", "Log(CSF Protein)", "Log(CSF Lactate)", 
        "Sqrt(Glucose Ratio)", "Log(CSF Lympho):HIV"))
mcmc_intervals_data(draw_2i_extIg, pars=vars(dplyr::starts_with("a"))) -> dtin
bayesplot::mcmc_intervals(draw_2i_extIg, pars=vars(dplyr::starts_with("a")), prob_outer=.95)+
  ggplot2::scale_y_discrete(limits = unique(rev(dtin$parameter)), labels=unclass(ylabs))

res_2i_extIg$summary(c("a0", "a", "z_Smear", "z_Mgit", "z_Xpert")) -> summarytab
res_2i_extIg$summary(c("z_Smear", "z_Mgit", "z_Xpert"), mean, ~quantile(.x, c(.025, .975))) %>% 
  mutate(across(where(is.numeric), nimble::ilogit)) %>% 
  mutate(Summary = glue::glue("{round(mean,2)} ({round(`2.5%`,2)}, {round(`97.5%`,2)})"))

thetas <- res_2i_extIg$summary("theta")
ntheta <- seq_along(nrow(data_19EI))
lowprobs <- cbind(data_19EI$USUBJID, thetas$mean, Xd, Xc)[thetas$mean < .1 & data_19EI$csf_smear,]
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
