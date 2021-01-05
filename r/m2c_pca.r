library(magrittr)
library(data.table)
library(dplyr)
source("r/include/functions.R")
rstan::rstan_options(auto_write=T, javascript=F)
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

# Create k-fold
N <- nrow(data_19EI)
K <- 20
set.seed(6120)
hh <- loo::kfold_split_random(K, N) 
foldkept <- matrix(1, nrow = N, ncol = K)
for(i in 1:N) foldkept[i, hh[i]] <- 0
holdout  <- 1- foldkept
foldkept <- split(foldkept,rep(1:ncol(foldkept),each=nrow(foldkept)))
holdout  <- split(holdout,rep(1:ncol(holdout),each=nrow(holdout)))
m2cpcakf_input <-
  list(
    N_all = nrow(data_19EI),
    nXc = ncol(Xc),
    nXd = ncol(Xd),
    nTd = ncol(Td),
    nTc = ncol(Tc),
    Y_Smear_all = data_19EI$csf_smear,
    Y_Mgit_all = data_19EI$csf_mgit,
    Y_Xpert_all = data_19EI$csf_xpert,
    Xc_all = Xc,
    Xd_all = Xd,
    Td_all = Td,
    Tc_all = Tc,
    obs_Xc_all = obs_Xc,
    obs_Xd_all = obs_Xd,
    obs_Td_all = obs_Td,
    obs_Tc_all = obs_Tc,
    K_pca = 2
  )
m2cpcakf_inputs <- pbmcapply::pbmclapply(seq_len(K), mc.cores=10, function(k) modifyList(m2cpcakf_input, list(keptin=foldkept[[k]])))

dir.create("outputs/m2cpca_kf", showWarnings = F)

ss <- repeated_stan_K_fold(
  file="stan/m2cpca_kf.stan", data =m2cpcakf_input,
  K = 20, N = nrow(data_19EI),
  repetition = 10,
  include_paths="stan/",
  chains=3, cores=18, output_dir="outputs/m2cpca_kf",
  adapt_delta=.88, init=1, seed=290693, max_treedepth=12,
  thin=2,
  iter_sampling=1000, iter_warmup=1200
)
# ss <- stan_kfold(file="stan/m2cpca_kf.stan",list_of_datas=m2cpcakf_inputs,
#                  include_paths="stan/",
#                  chains=3, cores=16, output_dir="outputs/m2cpca_kf",
#                  adapt_delta=.85, init=1, seed=29063, max_treedepth=12,
#                  thin=2,
#                  iter_sampling=1000, iter_warmup=1000)
ss_sampler <- rstan::stan_model("stan/m2cpca_kf.stan")
sstest <- rstan::sampling(ss_sampler, data=m2cpcakf_inputs[[10]],
                          chains=3, cores=3, control = list(adapt_delta=.85, max_treedepth=12), init_r=1, 
                          seed=29063,
                          thin=1,
                          # pars = c("a0", "a", "b_RE", "b_HIV", "z_Smear", "z_Mgit", "z_Xpert",
                          #          "log_lik", "p_Smear", "p_Mgit", "p_Xpert"),
                          iter=2000, warmup=1000)

ss <- stan_kfold(sampler = ss_sampler, #file="stan/m2cpca_kf.stan",
                 list_of_datas=m2cpcakf_inputs,
                 backend = "rstan",
                 chains=3, cores=18,
                 control = list(adapt_delta=.85, max_treedepth=12), init_r=1, 
                 seed=29063,
                 thin=1,
                 pars = c("a0", "a", "b_RE", "b_HIV", "z_Smear", "z_Mgit", "z_Xpert",
                          "log_lik", "p_Smear", "p_Mgit", "p_Xpert"),
                 iter=2000, warmup=1000)

m2cpcakf_sampler <- cmdstanr::cmdstan_model(stan_file="stan/m2cpca_kf.stan", include_paths = "stan/")
m2cpcakf_result <- m2cpcakf_sampler$sample(data = m2cpcakf_inputs[[10]], seed = 43904,
                                 output_dir = 'outputs/m2cpca_kf',
                                 chains = 1, parallel_chains = 3,
                                 iter_sampling = 1000, iter_warmup = 1000,
                                 init = 1,
                                 show_messages = T, save_warmup = F,
                                 adapt_delta = .85)
ee <- extract_log_lik_K(ss,holdout)
kk <- kfold(ee)

p <- extract_K_fold(ss, holdout, pars=c("p_Smear", "p_Mgit", "p_Xpert"), cores=18)
p_summary <- sapply(p, apply, 2, function(l) data.frame(mean = mean(l), median = median(l), CI2.5 = quantile(l, .25), CI97.5 = quantile(l, .975)), simplify = FALSE, USE.NAMES = TRUE)
p_summary <- sapply(p_summary, dplyr::bind_rows, simplify = FALSE, USE.NAMES = TRUE)

library(ggplot2); library(patchwork)

span = .8
calib_curve(p_summary$p_Smear$mean,m2cpcakf_input$Y_Smear_all, "Smear", span=span) + 
  calib_curve(p_summary$p_Mgit$mean,m2cpcakf_input$Y_Mgit_all, "Mgit", span=span) +
  calib_curve(p_summary$p_Xpert$mean,m2cpcakf_input$Y_Xpert_all, "Xpert", span=span)
