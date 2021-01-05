library(magrittr)
library(data.table)
library(dplyr)
source("r/include/functions.R")
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
set.seed(68270)
hh <- loo::kfold_split_random(K, N) 
foldkept <- matrix(1, nrow = N, ncol = K)
for(i in 1:N) foldkept[i, hh[i]] <- 0
holdout  <- 1- foldkept
foldkept <- split(foldkept,rep(1:ncol(foldkept),each=nrow(foldkept)))
holdout  <- split(holdout,rep(1:ncol(holdout),each=nrow(holdout)))
m2gkf_input <-
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
m2cgkf_inputs <- pbmcapply::pbmclapply(seq_len(K), mc.cores=10, function(k) modifyList(m2cgkf_input, list(keptin=foldkept[[k]])))

ss_sampler <- rstan::stan_model("stan/m2cg_kf.stan")
ss <- stan_kfold(sampler = ss_sampler,
                 list_of_datas=m2cgkf_inputs,
                 backend = "rstan",
                 chains=3, cores=19,
                 control = list(adapt_delta=.8, max_treedepth=13), init_r = 1, seed=290693,
                 pars = c("a0", "a", "b_RE", "b_HIV", "b_glu", "z_Smear", "z_Mgit", "z_Xpert",
                          "log_lik", "p_Smear", "p_Mgit", "p_Xpert", "theta"),
                 iter=1800, warmup=800)

ee <- extract_log_lik_K(ss,holdout)
kk <- kfold(ee)
p <- extract_K_fold(ss, holdout, pars=c("theta", "p_Smear", "p_Mgit", "p_Xpert"), cores=18)
p_summary <- sapply(p, apply, 2, function(l) data.frame(mean = mean(l), median = median(l), CI2.5 = quantile(l, .25), CI97.5 = quantile(l, .975)), simplify = FALSE, USE.NAMES = TRUE)
p_summary <- sapply(p_summary, dplyr::bind_rows, simplify = FALSE, USE.NAMES = TRUE)

library(ggplot2); library(patchwork)

span = .8
calib_curve(p_summary$p_Smear$mean,m2cgkf_input$Y_Smear_all, "Smear", span=span) + ylim(0, 1) + 
  calib_curve(p_summary$p_Mgit$mean,m2cgkf_input$Y_Mgit_all, "Mgit", span=span) + ylim(0, 1) +
  calib_curve(p_summary$p_Xpert$mean,m2cgkf_input$Y_Xpert_all, "Xpert", span=span) + ylim(0, 1)

ggplot(mapping=aes(x=p_summary$p_Mgit$mean, y=as.numeric(m2cgkf_input$Y_Mgit_all))) + 
  stat_smooth(method="glm", formula=y~splines::ns(x,3)) + 
  geom_line(aes(y=p_summary$p_Mgit$mean))

ggplot(mapping=aes(x=p_summary$p_Xpert$mean, y=as.numeric(m2cgkf_input$Y_Xpert_all))) + 
  stat_smooth(method="glm", formula=y~splines::ns(x,3)) + 
  geom_line(aes(y=p_summary$p_Xpert$mean))
