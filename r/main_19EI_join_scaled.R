library(data.table)
library(magrittr)

data_19EI <- readRDS('data/cleaned/data_19EI.RDS')
data_19EI <- na.omit(data_19EI, cols = c('csf_smear', 'csf_mgit', 'csf_xpert'))

Xd <- data_19EI %$% cbind(
  hiv_stat,                             #1
  clin_symptoms,                        #2
  clin_motor_palsy,                     #3
  clin_nerve_palsy,                     #4
  xray_pul_tb,                          #5
  xray_miliary_tb                       #6
)

Xc <- data_19EI %$% cbind(
  age,                                  #1    #7
  log2(clin_illness_day),               #2    #8
  15-clin_gcs,                          #3    #9
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
    Xc = scale(Xc),
    Xd = Xd,
    Tc = scale(Tc),
    Td = Td,
    obs_Xc = obs_Xc,
    obs_Xd = obs_Xd,
    obs_Tc = obs_Tc,
    obs_Td = obs_Td
  )

model_2i_extIg_scale <- cmdstanr::cmdstan_model('stan/model_2i_extIg_scaled.stan', dir = 'stan')
res_2i_extIg_scale   <- model_2i_extIg_scale$sample(
  data = model_input, seed = 29061993, output_dir = 'outputs',
  chains = 4, parallel_chains = 4,
  iter_sampling = 2000, iter_warmup = 2000,
  show_messages =FALSE,save_warmup = T,
  adapt_delta = .85, max_treedepth = 18)
