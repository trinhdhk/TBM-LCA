library(data.table)
library(magrittr)

data_19EI <- readRDS('data/cleaned/data_19EI.RDS')
data_19EI <- na.omit(data_19EI, cols = c('csf_smear', 'csf_mgit', 'csf_xpert'))
miss_summary <- dplyr::summarise_all(data_19EI, ~sum(is.na(.x)))


test.impute = rstan::stan_model('stan/test.impute.stan')

impute.input = list(
  N = nrow(data_19EI),
  HIV = data_19EI[, fcase(!is.na(hiv_stat), as.numeric(hiv_stat), default=0)],
  cs = as.data.frame(data_19EI[, .(fcase(!is.na(ISCOUGH), as.numeric(ISCOUGH), default=0), 
                     fcase(!is.na(ISNSWEAT),as.numeric(ISNSWEAT), default=0), 
                     fcase(!is.na(ISWEIGHT),as.numeric(ISWEIGHT), default=0))]),
  cbm_cs = fcase(!is.na(data_19EI$clin_symptoms),as.numeric(data_19EI$clin_symptoms),default=0),
  obs_HIV = !is.na(data_19EI$hiv_stat),
  obs_cs = as.data.frame(data_19EI[, .(as.numeric(!is.na(ISNSWEAT)), as.numeric(!is.na(ISCOUGH)), as.numeric(!is.na(ISWEIGHT)))]),
  obs_cbm_cs = !is.na(data_19EI$clin_symptoms)
)

sample.impute = rstan::sampling(test.impute, data=impute.input, chains=4, cores=4, iter=10000, seed=1208)
