library(magrittr)
library(data.table)
library(dplyr)
#source("r/include/functions.R")
#rstan::rstan_options(javascript=F)
#data_19EI <- readRDS('data/cleaned/data_19EI_3.RDS')
#data_19EI <- na.omit(data_19EI, cols = c('csf_smear', 'csf_mgit', 'csf_xpert'))
#data_19EI[222, "GCSV"] <- 5
#data_19EI[, clin_contact_tb := clin_contact_tb %in% TRUE]

source("r/include/functions.R")
rstan::rstan_options(javascript=F, autowrite=T)

data_19EI <- readRDS('data/cleaned/data_19EI_3.RDS')
# data_19EI <- na.omit(data_19EI, cols = c('csf_smear', 'csf_mgit'))
data_19EI[, `:=`(
  csf_smear = ifelse(is.na(csf_smear), 0, csf_smear),
  csf_mgit = ifelse(is.na(csf_mgit)  & !is.na(csf_smear), 0, csf_mgit),
  csf_xpert = ifelse(is.na(csf_xpert) & !is.na(csf_smear), 0, csf_xpert)
)]
data_19EI[is.na(GCSV) & GLASCOW == 15, GCSV := 5]
data_19EI[, clin_contact_tb := clin_contact_tb %in% TRUE]
data_19EI <- data_19EI[USUBJID!="003-335"]

Xd <- data_19EI %$% cbind(
  hiv_stat,                             #1
  clin_symptoms,                        #2
  clin_motor_palsy,                     #3
  clin_nerve_palsy,                     #4
  clin_contact_tb,
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


Tc <- data_19EI %$% cbind(4-GCSE, 6-GCSM, 5-GCSV, log(WHITE), log(LYMP))

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
K <- 10
set.seed(68730)
hh <- loo::kfold_split_random(K, N) 
foldkept <- matrix(1, nrow = N, ncol = K)
for(i in 1:N) foldkept[i, hh[i]] <- 0
holdout  <- 1- foldkept
foldkept <- split(foldkept,rep(1:ncol(foldkept),each=nrow(foldkept)))
holdout  <- split(holdout,rep(1:ncol(holdout),each=nrow(holdout)))
m2ckf_input <-
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
    obs_Tc_all = obs_Tc
  )
#m2ckf_inputs <- pbmcapply::pbmclapply(seq_len(K), mc.cores=10, function(k) modifyList(m2ckf_input, list(keptin=foldkept[[k]])))
m2ckf_folds <- repeated_kfold(m2ckf_input, K = 10, N_rep = 4, N_obs = nrow(data_19EI), seed = 6012)
m2ckf_inputs <- m2ckf_folds$inputs
dir.create("outputs/m2c_kf10", showWarnings = F)
m2ckf_sampler <- rstan::stan_model("stan/m2c_kf.stan")
m2ckf_outputs <- stan_kfold(file = "stan/m2c_kf.stan",
                            list_of_datas=m2ckf_inputs,
                            backend = "cmdstanr",
                            include_paths = "./stan", 
                            merge = F,
                            output_dir = 'outputs/m2c_kf10',
                            chains=3, cores=18, 
                            thin = 1,
                            adapt_delta=.7, max_depth=12, init = 1, seed=29693,
                            #pars = c("a0", "a", "b_RE", "b_HIV", "z_Smear", "z_Mgit", "z_Xpert",
                              #       "log_lik", "p_Smear", "p_Mgit", "p_Xpert", "theta"),
                            iter_sampling=1000, iter_warmup=1000, save_warmup=F)

m2ckf_outputs <- stan_kfold(sampler = m2ckf_sampler,
                 list_of_datas=m2ckf_inputs,
                 backend = "rstan",
                 chains=3, cores=19, 
                 thin = 1, 
                 merge = T,
                 control = list(adapt_delta=.75, max_treedepth=12), 
                 init_r = 1, seed=2693,
                 sample = "outputs/m2c_kf10/",
                 pars = c("a0", "a", "b_RE", "b_HIV", "z_Smear", "z_Mgit", "z_Xpert",
                          "log_lik", "p_Smear", "p_Mgit", "p_Xpert", "theta"),
                 iter=2000, warmup=1000)
# sstest=rstan::sampling(ss_sampler, d3ata=m2ckf_inputs[[1]], chains=2, cores=3)

holdout = m2ckf_folds$holdout
m2ckf_loglik <- extract_log_lik_K(m2ckf_outputs,holdout)
m2ckf_elpd <- kfold(m2ckf_loglik)
p <- extract_K_fold(m2ckf_outputs, holdout, pars=c("theta", "p_Smear", "p_Mgit", "p_Xpert"))
p_summary <- sapply(p, apply, 2, function(l) data.frame(mean = mean(l), median = median(l), CI2.5 = quantile(l, .25), CI97.5 = quantile(l, .975)), simplify = FALSE, USE.NAMES = TRUE)
p_summary <- sapply(p_summary, dplyr::bind_rows, simplify = FALSE, USE.NAMES = TRUE)

library(ggplot2); library(patchwork)

span = 1
calib_curve(p_summary$p_Smear$mean,m2ckf_input$Y_Smear_all, "Smear", span=span)  + 
  calib_curve(p_summary$p_Mgit$mean,m2ckf_input$Y_Mgit_all, "Mgit", span=span) +
  calib_curve(p_summary$p_Xpert$mean,m2ckf_input$Y_Xpert_all, "Xpert", span=span) 

ggplot(mapping=aes(x=p_summary$p_Smear$mean, y=as.numeric(m2ckf_input$Y_Smear_all))) + 
  stat_smooth(method="glm", formula=y~splines::ns(x,5)) + 
  geom_line(aes(y=p_summary$p_Smear$mean)) +
ggplot(mapping=aes(x=p_summary$p_Mgit$mean, y=as.numeric(m2ckf_input$Y_Mgit_all))) + 
  stat_smooth(method="glm", formula=y~splines::ns(x,5)) + 
  geom_line(aes(y=p_summary$p_Mgit$mean)) + 
ggplot(mapping=aes(x=p_summary$p_Xpert$mean, y=as.numeric(m2ckf_input$Y_Xpert_all))) + 
  stat_smooth(method="glm", formula=y~splines::ns(x,5)) + 
  geom_line(aes(y=p_summary$p_Xpert$mean))

p2 <- extract_K_fold(m2ckf_outputs, holdout, pars=c("theta", "p_Smear", "p_Mgit", "p_Xpert"), holdout=0)
p2_summary <- sapply(p2, apply, 2, function(l) data.frame(mean = mean(l), median = median(l, na.rm=T), CI2.5 = quantile(l, .25, na.rm=T), CI97.5 = quantile(l, .975, na.rm=T)), simplify = FALSE, USE.NAMES = TRUE)
p2_summary <- sapply(p2_summary, dplyr::bind_rows, simplify = FALSE, USE.NAMES = TRUE)


ggplot(mapping=aes(x=p2_summary$p_Smear$mean, y=as.numeric(m2ckf_input$Y_Smear_all))) + 
  stat_smooth(method="loess", formula=y~x, span=span) + 
  stat_smooth(mapping=aes(x=p_summary$p_Smear$mean), 
              method="loess", formula=y~x, color="red", fill="red", span=span) + 
  geom_line(aes(y=p2_summary$p_Smear$mean)) +
  # ylim(0,1)+
  ggplot(mapping=aes(x=p2_summary$p_Mgit$mean, y=as.numeric(m2ckf_input$Y_Mgit_all))) + 
  stat_smooth(method="loess", formula=y~x, span=span) + 
  stat_smooth(mapping=aes(x=p_summary$p_Mgit$mean), 
              method="loess", formula=y~x, color="red", fill="red", span=span) +
  geom_line(aes(y=p2_summary$p_Mgit$mean)) + 
  # ylim(0,1)+
  ggplot(mapping=aes(x=p2_summary$p_Xpert$mean, y=as.numeric(m2ckf_input$Y_Xpert_all))) + 
  stat_smooth(method="loess", formula=y~x, span=span) + 
  stat_smooth(mapping=aes(x=p_summary$p_Xpert$mean), 
              method="loess", formula=y~x, color="red", fill="red", span=span) +
  geom_line(aes(y=p2_summary$p_Xpert$mean))

m2ckf_output <- rstan::sflist2stanfit(m2ckf_outputs)
bayesplot::mcmc_intervals(m2ckf_output, reg="^[a,b]")
bayesplot::mcmc_intervals(m2ckf_output, reg="^z", transformations = nimble::ilogit)
