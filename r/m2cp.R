library(magrittr)
library(data.table)
data_19EI <- readRDS('data/cleaned/data_19EI.RDS')
data_19EI <- na.omit(data_19EI, cols = c('csf_smear', 'csf_mgit', 'csf_xpert'))
data_19EI[222, "GCSV"] <- 5

Xd <- data_19EI %$% cbind(
  hiv_stat,                             #1
  clin_symptoms,                        #2
  clin_motor_palsy,                     #3
  clin_nerve_palsy,                     #4
  xray_pul_tb,                          #5
  xray_miliary_tb                       #6
)

Xc <- data_19EI %$% cbind(
  age=(age-mean(age, na.rm=TRUE))/10,          #1    #7
  id=log2(clin_illness_day),                   #2    #8
  glu=sqrt(BLDGLU),                            #3    #9 
  csfglu=sqrt(csf_glucose),                    #4    #10
  csflym=log2(csf_lympho+1),                   #5    #11   
  csfpro=log2(csf_protein),                    #6    #12   
  csflac=log2(csf_lactate),                    #7    #13   
  gcs=15-clin_gcs                              #8    #15             
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


Tc <- data_19EI %$% cbind(4-GCSE, 6-GCSM, 5-GCSV)

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

m2cp_input <-
  list(
    N = nrow(data_19EI),
    nXc = ncol(Xc),
    nXd = ncol(Xd),
    nTd = ncol(Td),
    nTc = ncol(Tc),
    Y_Smear = data_19EI$csf_smear,
    Y_Mgit = data_19EI$csf_mgit,
    Y_Xpert = data_19EI$csf_xpert,
    Xc = Xc,
    Xd = Xd,
    Td = Td,
    Tc = Tc,
    obs_Xc = obs_Xc,
    obs_Xd = obs_Xd,
    obs_Td = obs_Td,
    obs_Tc = obs_Tc
  )

m2cp_sampler <- cmdstanr::cmdstan_model('stan/m2cp.stan', dir = 'stan/bin')
dir.create("outputs/m2cp", showWarnings = F)
m2cp_result <- m2cp_sampler$sample(data = m2cp_input, seed = 43904, 
                                 output_dir = 'outputs/m2cp',
                                 chains = 10, parallel_chains = 10,
                                 iter_sampling = 1200, iter_warmup = 1200,
                                 init = 1,
                                 show_messages = T, save_warmup = T,
                                 adapt_delta = .85, max_treedepth = 18)

m2cp_result$cmdstan_diagnose()
m2cp_draw <- m2cp_result$draws(c("a0", "a", "b_HIV", "b_RE", "z_Smear", "z_Mgit", "z_Xpert", "HIV_a0", "X", "gcsv_imp"))
bayesplot::mcmc_intervals(m2cp_draw, pars=dplyr::vars(dplyr::starts_with("a")), prob_outer = .95)

######## 

m2ca_sampler <- cmdstanr::cmdstan_model('stan/m2ca.stan', dir = 'stan/bin')
dir.create("outputs/m2ca", showWarnings = F)
m2ca_result <- m2ca_sampler$sample(data = m2cp_input, seed = 43904, 
                                   output_dir = 'outputs/m2ca',
                                   chains = 10, parallel_chains = 10,
                                   iter_sampling = 1000, iter_warmup = 1000,
                                   init = 1,
                                   show_messages = T, save_warmup = T,
                                   adapt_delta = .85, max_treedepth = 18)

m2ca_result$cmdstan_diagnose()
m2ca_draw <- m2ca_result$draws(c("a0", "a", "b_HIV", "b_RE", "z_Smear", "z_Mgit", "z_Xpert", "HIV_a0", "X", "gcsv_imp"))
bayesplot::mcmc_intervals(m2ca_draw, pars=dplyr::vars(dplyr::starts_with("a")), prob_outer = .95)
