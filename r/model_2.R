# Model with random effect
# Author trinhdhk
# Version: 0.1.2004

# 0. Libs
rm(list = ls())
library(data.table)
library(magrittr)
library(rstan)
options(mc.cores = 6)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
rstan_options(auto_write = TRUE)

# 1. Load data
cleaned_data_dir <- 'data/cleaned'
data_23TB <- readRDS(file.path(cleaned_data_dir, 'data_23TB.RDS'))
data_05TB <- readRDS(file.path(cleaned_data_dir, 'data_05TB.RDS'))
maindt <- rbind(data_23TB,
                data_05TB[,names(data_23TB),with = FALSE])[,-1]
maindt[, lym_glu_ratio := (csf_lym_pct/glucose_ratio)/100]
maindt[, `:=`(
  clin_symptoms = as.numeric(clin_symptoms),
  clin_nerve_palsy = as.numeric(clin_nerve_palsy),
  clin_motor_palsy = as.numeric(clin_motor_palsy),
  csf_clear = as.numeric(csf_clear),
  xray_miliary_tb = as.numeric(xray_miliary_tb),
  xray_pul_tb = as.numeric(xray_pul_tb)
)]

# 2. Build input

maindtcomplete <- na.omit(maindt, cols = c('age', 'sex', 'bmi', 'csf_smear', 'csf_mgit', 'csf_xpert', 'lym_glu_ratio', #'img_score'
                                           'hiv_stat','clin_illness_day', 'clin_symptoms', 'clin_contact_tb', 'clin_gcs', 
                                           'clin_nerve_palsy', 'clin_motor_palsy', 'csf_clear',  'csf_lympho', #'csf_protein',
                                           'xray_miliary_tb', 'xray_pul_tb'))
X <- maindtcomplete %$% as.matrix(cbind(hiv_stat, age, sex, bmi, clin_illness_day, clin_contact_tb,
                                        clin_nerve_palsy, clin_motor_palsy, glucose_ratio, csf_lympho,
                                        xray_pul_tb))

model_input <- maindtcomplete %$% 
  list(
    N = nrow(maindtcomplete), nX = ncol(X),
    Y_Smear = csf_smear,
    Y_Mgit = csf_mgit,
    Y_Xpert = csf_xpert,
    Y_Img = xray_miliary_tb,
    Y_lnLymGlu = log(lym_glu_ratio+1),
    X = X
  )

model_2 <- stan_model('stan/model_2.stan')
result_2 <- sampling(model_2, data=model_input, chain=6, iter=8000, seed=128, warmup=1000, control = list(max_treedepth = 18, adapt_delta=.97))
extract_2 <- extract(result_2)
# plot((extract_1$theta), extract_1$a[,1])

library(bayesplot)
posterior_2 <- as.array(result_2)
np_2 <- nuts_params(result_2)
mcmc_parcoord(posterior_2, np = np_2, pars=dplyr::vars(-lp__))
mcmc_pairs(posterior_2, np = np_2, pars = c("z_Mgit[1]", 'z_Img[1]'),
           off_diag_args = list(size = 0.75))
