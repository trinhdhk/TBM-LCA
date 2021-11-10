library(rstan)
library(magrittr)
library(data.table)
data_19EI <- readRDS('data/cleaned/data_19EI.RDS')
data_19EI <- na.omit(data_19EI, cols = c('csf_smear', 'csf_mgit', 'csf_xpert'))
data_19EI[, clin_reducedconsciousness := (ISREDUCED %in% TRUE) | (clin_gcs < 15 %in% TRUE)]


options(mc.cores = 6)
# Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
rstan_options(auto_write = TRUE)

Xd <- data_19EI %$% cbind(
  hiv_stat,                             #1
  clin_symptoms,                        #2
  clin_motor_palsy,                     #3
  clin_nerve_palsy,                     #4
  xray_pul_tb,                          #5
  xray_miliary_tb,                      #6
  clin_reducedconsciousness             #7
)

Xc <- data_19EI %$% cbind(
  age=(age-mean(age, na.rm=TRUE))/10,          #1    #8
  id=log2(clin_illness_day),                   #2    #9
  glu=sqrt(BLDGLU),                            #3    #10 
  csfglu=sqrt(csf_glucose),                    #4    #11
  csflym=log2(csf_lympho+1),                   #5    #12
  csfpro=log2(csf_protein),                    #6    #13
  csflac=log2(csf_lactate)                     #7    #14
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

is.not.na <- Negate(is.na)
obs_Xd <- apply(Xd, 2, is.not.na) %>% apply(2, as.numeric)
obs_Xc <- apply(Xc, 2, is.not.na) %>% apply(2, as.numeric)
obs_Td <- apply(Td, 2, is.not.na) %>% apply(2, as.numeric)

my.replace_na <- function(data, replace=0){
  apply(data, 2, tidyr::replace_na, replace=replace)
}

Xd <- my.replace_na(Xd)
Xc <- my.replace_na(Xc)
Td <- my.replace_na(Td)

model_input <-
  list(
    N = nrow(data_19EI),
    nXc = ncol(Xc),
    nXd = ncol(Xd),
    nTd = ncol(Td),
    Y_Smear = data_19EI$csf_smear,
    Y_Mgit = data_19EI$csf_mgit,
    Y_Xpert = data_19EI$csf_xpert,
    Xc = Xc,
    Xd = Xd,
    Td = Td,
    obs_Xc = obs_Xc,
    obs_Xd = obs_Xd,
    obs_Td = obs_Td
  )

# model_input.boot <-
#   list(
#     N = nrow(data_19EI),
#     nXc = ncol(Xc),
#     nXd = ncol(Xd),
#     nTd = ncol(Td),
#     Y_Smear0 = data_19EI$csf_smear,
#     Y_Mgit0 = data_19EI$csf_mgit,
#     Y_Xpert0 = data_19EI$csf_xpert,
#     Xc0 = Xc,
#     Xd0 = Xd,
#     Td0 = Td,
#     obs_Xc0 = obs_Xc,
#     obs_Xd0 = obs_Xd,
#     obs_Td0 = obs_Td,
#     resample=T
#   )

model_2i_extI_dGCS <- cmdstanr::cmdstan_model('stan/model_2i_noGCS.stan', dir = 'stan')
model_2i_extI_dGCS <- cmdstanr::cmdstan_model('stan/m2i.stan', dir = 'stan')
model_2i_extI_dGCS_lasso <- cmdstanr::cmdstan_model('stan/m2ilasso.stan', dir = 'stan')
model_2i_extI_dGCS.boot <- cmdstanr::cmdstan_model('stan/model_2i_noGCS_boot.stan', dir = 'stan')

res_2i_extI_dGCS   <- model_2i_extI_dGCS$sample(data = model_input, seed = 29061913, output_dir = 'outputs',
                                      chains = 12, parallel_chains = 12,
                                      iter_sampling = 1500, iter_warmup = 1500,
                                      init = 1,
                                      show_messages = T, save_warmup = T,
                                      adapt_delta = .88, max_treedepth = 18)

res_2i_extI_dGCS_lasso   <- model_2i_extI_dGCS_lasso$sample(data = model_input, seed = 2906193, output_dir = 'outputs',
                                                chains = 12, parallel_chains = 12,
                                                iter_sampling = 2500, iter_warmup = 1500,
                                                init = 2,
                                                show_messages = T, save_warmup = TRUE,
                                                adapt_delta = .85, max_treedepth = 18)

res_2i_extI_dGCS_lasso_sum <- res_2i_extI_dGCS_lasso$summary()

draw_2i_extI_dGCS <- res_2i_extI_dGCS$draws(c("a0", "a", "b_HIV", "b_RE", "z_Smear", "z_Mgit", "z_Xpert", "HIV_a0"))
draw_2i_extI_dGCS_lasso <- res_2i_extI_dGCS_lasso$draws(c("a0", "a", "b_HIV", "b_RE", "z_Smear", "z_Mgit", "z_Xpert", "HIV_a0"))
bayesplot::mcmc_intervals(draw_2i_extI_dGCS, pars=vars(dplyr::starts_with("a")), prob_outer = .95)
bayesplot::mcmc_intervals(draw_2i_extI_dGCS_lasso, pars=dplyr::vars(dplyr::starts_with("a")))

thetas <- res_2i_extI_dGCS$summary("theta")
ntheta <- seq_along(nrow(data_19EI))
highprobs <- cbind(thetas$mean, Xc)[thetas$mean > .7,]
highprobs[,2] <- highprobs[,2]*10
highprobs[,c(3, 6:8)] <- 2^highprobs[,c(3,6:8)]
highprobs[,4:6] <- highprobs[,4:6]^2

lowprobs <- cbind(thetas$mean, Xc, Xd, Y_Smear = data_19EI$csf_smear,
                  Y_Mgit = data_19EI$csf_mgit,
                  Y_Xpert = data_19EI$csf_xpert)[thetas$mean < .1,]
lowprobs[,2] <- lowprobs[,2]*10
lowprobs[,c(3, 6:8)] <- 2^lowprobs[,c(3,6:8)]

model_2i_extI_dGCS_stan <- rstan::stan_model('stan/m2i.stan')
res_2i_extI_dGCS.2   <- rstan::sampling(model_2i_extI_dGCS_stan, chain=4, cores=4, seed = 2906,
                                        data = model_input, iter = 8000, warmup = 1500,
                                        control = list(adapt_delta = .88, max_treedepth = 18),
                                        init_r = 1,
                                        include = TRUE)

model_2i_extI_dGCS <- rstan::stan_model('stan/model_2i_noGCS.stan', model_name="model_2i_noGCS")
res_2i_extI_dGCS.2   <- rstan::sampling(model_2i_extI_dGCS, chain=3, cores=3, seed = 2906,
                                        data = model_input, iter = 18000, warmup = 2000,
                                        control = list(adapt_delta = .85, max_treedepth = 18),
                                        init_r = 1,
                                        include = TRUE)
pars=c("a0", "a", "b_HIV", "b_RE", "z_Smear", "z_Mgit", "z_Xpert", "HIV_a0", "log_lik"))

library(future)
future::plan(multisession, workers=6)
boot.fit = vector('list', 100)
for (i in 1:100){
  boot.fit[[i]] <- future({model_2i_extI_dGCS.boot$sample(data = model_input.boot, seed = 29061993, output_dir = 'outputs',
                                             chains = 3, parallel_chains = 3,
                                             iter_sampling = 1500, iter_warmup = 1200,
                                             show_messages = FALSE,
                                             adapt_delta = .88, max_treedepth = 18)})
}

res = rstan::stan('stan/model_2i_noGCS.stan',data = model_input, iter = 3000, warmup = 1000,
                  # include = TRUE,
                  # pars=c("a0", "a", "b_HIV", "b_RE", "z_Smear", "z_Mgit", "z_Xpert", "HIV_a0"),
                  cores = 4)


stancode <- 'data {real y_mean;} parameters {real y;} model {y ~ normal(y_mean,1);}'
mod <- stan_model(model_code = stancode, verbose = TRUE)
fit <- sampling(mod, data = list(y_mean = 0))
fit2 <- sampling(mod, data = list(y_mean = 5))