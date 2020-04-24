# Model without random effect
# Author trinhdhk
# Version: 0.1.2004

# 0. Libs
rm(list = ls())
library(data.table)
library(magrittr)
library(rstan)
options(mc.cores = 4)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')

# 1. Load data
cleaned_data_dir <- 'data/cleaned'
data_23TB <- readRDS(file.path(cleaned_data_dir, 'data_23TB.RDS'))
data_05TB <- readRDS(file.path(cleaned_data_dir, 'data_05TB.RDS'))
maindt <- rbind(data_23TB,
                data_05TB[,names(data_23TB),with = FALSE])[,-1]
maindt[, lym_glu_ratio := csf_lympho/glucose_ratio]
maindt[, `:=`(
  hiv_stat = hiv_stat=='Yes',
  clin_symptoms = as.numeric(clin_symptoms),
  clin_nerve_palsy = as.numeric(clin_nerve_palsy),
  clin_motor_palsy = as.numeric(clin_motor_palsy),
  csf_clear = as.numeric(csf_clear),
  xray_miliary_tb = as.numeric(xray_miliary_tb),
  xray_pul_tb = as.numeric(xray_pul_tb)
)]

# 2. Build input

maindtcomplete <- na.omit(maindt, cols = c('csf_smear', 'csf_mgit', 'csf_xpert', 'img_score', 'lym_glu_ratio', 
                         'hiv_stat','clin_illness_day', 'clin_symptoms', 'clin_contact_tb', 'clin_gcs', 
                         'clin_nerve_palsy', 'clin_motor_palsy', 'csf_clear',  'csf_lympho', 'csf_glucose',
                         'csf_protein', 'xray_miliary_tb', 'xray_pul_tb'))
model_input <- maindtcomplete %$% 
  list(
    N = nrow(maindtcomplete),
    Y_Smear = csf_smear,
    Y_Mgit = csf_mgit,
    Y_Xpert = csf_xpert,
    Y_Img = img_score+1,
    Y_lnLymGlu = log(lym_glu_ratio),
    X = as.matrix(cbind(hiv_stat, clin_illness_day, clin_symptoms, clin_contact_tb, clin_gcs, 
              clin_nerve_palsy, clin_motor_palsy, csf_clear,  csf_lympho, csf_glucose,
              csf_protein, xray_miliary_tb, xray_pul_tb))
  )

model_1 <- stan_model('stan/model_1.stan')
result_1 <- sampling(model_1, data=model_input, chain=1, iter=4000, seed=296, warmup=1000, control = list(max_treedepth = 18, adapt_delta=.9))
