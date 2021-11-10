# Load libs
rm(list = ls())
library(data.table)
library(magrittr)
library(rstan)
library(bayesplot)
options(mc.cores = 4)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
rstan_options(auto_write = TRUE)
source('r/fn_def.R')

data_23TB <- readRDS('data/cleaned/data_23TB.RDS')
N <- nrow(data_23TB)
data_23TBcomplete <- na.omit(data_23TB, cols = c('age', 'sex', 'bmi', 'csf_smear', 'csf_mgit', 'csf_xpert', 'img_score',
                                                 'clin_illness_day', 'clin_symptoms', 'clin_gcs', 
                                                 'clin_nerve_palsy', 'clin_motor_palsy', 'csf_lympho', 'glucose_ratio', 'csf_clear', 'csf_protein',
                                                 'xray_miliary_tb', 'xray_pul_tb'))
X <- data_23TBcomplete %$% as.matrix(cbind(age, 
                                           sex, 
                                           bmi, 
                                           clin_illness_day, #remove contact_tb due to many "unknown"s
                                           clin_nerve_palsy,
                                           clin_motor_palsy,
                                           csf_lympho,
                                           glucose_ratio,
                                           csf_protein,
                                           xray_pul_tb))
simul_input <- data_23TBcomplete %$%
  list(
    N = nrow(data_23TBcomplete), 
    nX = ncol(X),
    a0 = -.5,
    a = c(.01, 0, 0.001, .001, .02, 0, 0, -2, .5, .05),
    b_age = .2,
    b = c(.1, .5, .8),
    z_Smear = qnorm(c(.05, .7)), z_Mgit = qnorm(c(0, .6)), z_Xpert = qnorm(c(.01, .8)), z_Img = qnorm(c(.02, .5)), lambda_brainImg = c(0, 1),
    X=X
  )

simul_model_3h3a <-  load_model('simul_3h3a', data=simul_input, chain = 1, iter=5000, seed=1708, warmup=2000, 
                                control = list(max_treedepth = 18, adapt_delta=.99), init_r = .5)
es3h3a = extract(simul_model_3h3a)
i <- round(runif(1, 1, 3000))

model_input <- list(
  N = 50*nrow(data_23TBcomplete), 
  nX = ncol(X),
  Y_Smear = as.vector(es3h3a$Y_Smear[1:50,]),
  Y_Xpert = as.vector(es3h3a$Y_Xpert[1:50,]),
  Y_Mgit = as.vector(es3h3a$Y_Mgit[1:50,]),
  Y_Img = as.vector(es3h3a$Y_Img[1:50,]),
  Y_brainImg = as.vector(es3h3a$Y_brainImg[1:50,]),
  X=apply(X,2,rep,50)
)

model_3h3a <-  load_model('disc_3h3a', data=model_input, chain = 3, iter=50000, seed=1708, warmup=25000, 
                          control = list(max_treedepth = 18, adapt_delta=.99), init_r = .5)
ed3h3a = extract(model_3h3a)
