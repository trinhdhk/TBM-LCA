library(magrittr); library(data.table); library(dplyr)
data_23TB <- readRDS("data/cleaned/data_23TB.RDS")
# load("zzzz.RDS")
source("r/include/functions.R")

data_23TB <- na.omit(data_23TB, cols = c('csf_smear', 'csf_mgit', 'csf_xpert', 'csf_lactate','xray_miliary_tb', 'xray_pul_tb'))
summarise(data_23TB %>% select(-diag), across(everything(), ~ sum(is.na(.x))))
data_23TB[, bldglu := csf_glucose / glucose_ratio]

X <- data_23TB %$% cbind(hiv_stat, clin_symptoms, clin_motor_palsy, clin_nerve_palsy, xray_pul_tb, xray_miliary_tb,
                         age=(age-42.77)/10, id=log2(clin_illness_day), 
                         bldglu=sqrt(bldglu), csfglu = sqrt(csf_glucose), csflym=log2(csf_lympho+1),csfpro=log2(csf_protein),  csflac=log2(csf_lactate),sqrt(glucose_ratio),  gcs=log2(15-clin_gcs+1))

a0 <- rstan::summary(ss2, "a0")$summary[,"mean"]
a <- rstan::summary(ss2, "a")$summary[,"mean"]
b_RE <- rstan::summary(ss2, "b_RE")$summary[,"mean"]
b_HIV <- rstan::summary(ss2, "b_HIV")$summary[,"mean"]
z_Xpert <- rstan::summary(ss2, "z_Xpert")$summary[, "mean"]
z_Mgit <- rstan::summary(ss2, "z_Mgit")$summary[, "mean"]
z_Smear <- rstan::summary(ss2, "z_Smear")$summary[, "mean"]

n_rep = 10000
p_Xpert <- p_Mgit <- p_Smear <- matrix(numeric(), ncol = n_rep, nrow=nrow(data_23TB))

z_theta <- a0 + X %*% a
theta <- nimble::ilogit(z_theta)

for (i in seq_len(n_rep)){
  RE <- rnorm(nrow(data_23TB), 0, 1)
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
span=.8
calib_curve(rowMeans(p_Smear),data_23TB$csf_smear, "Smear", span=span) + ylim(0, 1) + 
  calib_curve(rowMeans(p_Mgit),data_23TB$csf_mgit, "Mgit", span=span) + ylim(0, 1) +
  calib_curve(rowMeans(p_Xpert),data_23TB$csf_xpert, "Xpert", span=span) + ylim(0, 1)

