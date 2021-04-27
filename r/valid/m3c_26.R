library(magrittr)
library(data.table)
library(dplyr)

source("r/include/functions.R")
rstan::rstan_options(javascript=F)

data_26TB <- readRDS('data/cleaned/data_26TB.RDS')
mice_26TB <- mice::mice(data_26TB %>% select(-csf_lym_pct, -hiv_stat, -tube3, -USUBJID, -clin_symptoms, -clin_motor_palsy,-csf_smear, -csf_mgit, -csf_xpert), m = 10, maxit=50)

data_imp_26TB <- mice::complete(mice_26TB, 4) %>% mutate(
  csf_lym_pct = csf_lympho/csf_wbc, 
  hiv_stat = 0, 
  clin_symptoms = ISWEIGHT | ISNSWEAT | ISCOUGH,
  clin_motor_palsy = HEMIPLEGIA | TETRAPLEGIA | PARAPLEGIA,
  csf_xpert = data_26TB$csf_xpert,
  csf_mgit = data_26TB$csf_mgit,
  csf_smear = data_26TB$csf_smear
) %>% filter(!is.na(data_26TB$csf_lactate))
# data_26TB <- na.omit(data_26TB, cols = c('csf_smear', 'csf_mgit', 'csf_xpert', 'csf_lactate','xray_miliary_tb', 'xray_pul_tb'))

#data_26TB[,csf_lympho := csf_lym_pct*csf_wbc]
X <- data_imp_26TB %$% cbind(hiv_stat, clin_symptoms, clin_motor_palsy, clin_nerve_palsy, clin_contact_tb, xray_pul_tb, xray_miliary_tb,
                         age=(age-44)/10, id=log2(clin_illness_day+1), 
                         bldglu=sqrt(bld_glucose), csfglu = sqrt(csf_glucose), csflym=log2(csf_lympho+1),csfpro=log2(csf_protein),  csflac=log2(csf_lactate),sqrt(glucose_ratio),  gcs=log2(15-clin_gcs+1))
a0 <- rstan::summary(m3ckf_output, "a0")$summary[,"mean"]
a <- rstan::summary(m3ckf_output, "a")$summary[,"mean"]
b_RE <- rstan::summary(m3ckf_output, "b_RE")$summary[,"mean"]
b_HIV <- rstan::summary(m3ckf_output, "b_HIV")$summary[,"mean"]
z_Xpert <- rstan::summary(m3ckf_output)$summary[, "mean"]
z_Mgit <- rstan::summary(m3ckf_output, "z_Mgit")$summary[, "mean"]
z_Smear <- rstan::summary(m3ckf_output, "z_Smear")$summary[, "mean"]

n_rep = 10000
p_Xpert <- p_Mgit <- p_Smear <- matrix(numeric(), ncol = n_rep, nrow=nrow(data_imp_26TB))

z_theta <- a0 + X %*% a
theta <- nimble::ilogit(z_theta)

for (i in seq_len(n_rep)){
  RE <- rnorm(nrow(data_imp_26TB), 0, 1)
  z_Xpert_RE <- cbind(z_Xpert[1], z_Xpert[2] + b_RE[3] * (RE+ X[,1]*b_HIV))
  z_Mgit_RE <- cbind(z_Mgit[1], z_Mgit[2] + b_RE[2] * (RE+ X[,1]*b_HIV))
  z_Smear_RE <- cbind(z_Smear[1], z_Smear[2] + b_RE[1] * (RE+ X[,1] * b_HIV))
  
  p_Xpert_RE <- nimble::ilogit(z_Xpert_RE)
  p_Mgit_RE <- nimble::ilogit(z_Mgit_RE)
  p_Smear_RE <- nimble::ilogit(z_Smear_RE)
  
  
  p_Xpert[,i] <- (1-theta)*p_Xpert_RE[,1] + theta*p_Xpert_RE[,2]
  p_Mgit[,i] <- (1-theta)*p_Mgit_RE[,1] + theta*p_Mgit_RE[,2]
  p_Smear[,i] <- (1-theta)*p_Smear_RE[,1] + theta*p_Smear_RE[,2]
}

library(ggplot2); library(patchwork)
span=.75
calib_curve(rowMedians(p_Smear)[!is.na(rowMedians(p_Smear))],data_imp_26TB$csf_smear[!is.na(rowMedians(p_Smear))], "Smear", span=span)  + 
  calib_curve(rowMedians(p_Mgit)[!is.na(rowMedians(p_Mgit))]*3,data_imp_26TB$csf_mgit[!is.na(rowMedians(p_Mgit))], "Mgit", span=span)  +
  calib_curve(rowMedians(p_Xpert)[!is.na(rowMedians(p_Xpert))],data_imp_26TB$csf_xpert[!is.na(rowMedians(p_Xpert))], "Xpert", span=span) 

