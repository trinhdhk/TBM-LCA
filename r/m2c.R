library(magrittr)
library(data.table)
library(dplyr)
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

m2c_input <-
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

m2c_sampler <- cmdstanr::cmdstan_model('stan/m2c.stan', dir = 'stan/bin')
dir.create("outputs/m2c", showWarnings = F)
m2c_result <- m2c_sampler$sample(data = m2c_input, seed = 261993, 
                                 output_dir = 'outputs/m2c',
                                 chains = 10, parallel_chains = 10,
                                 iter_sampling = 1000, iter_warmup = 1200,
                                 init = 1,
                                 show_messages = T, save_warmup = T,
                                 adapt_delta = .86, max_treedepth = 18)

m2c_result$cmdstan_diagnose()
m2c_draw <- m2c_result$draws(c("a0", "a", "b_HIV", "b_RE", "z_Smear", "z_Mgit", "z_Xpert", "HIV_a0", "p_Xpert", "p_Mgit", "p_Smear", "X"))

mcmc_labels <- function(pl, ..., from_top=TRUE){
  labs <- c(...)
  params <- pl$data$parameter
  stopifnot(length(labs) == length(params))
  pl + ggplot2::scale_y_discrete(limits = if(from_top) rev(params) else params,
                                 labels = if(from_top) rev(labs) else labs)
}

pl <- bayesplot::mcmc_intervals(m2c_draw, 
                          pars=dplyr::vars(dplyr::starts_with("a")), 
                          prob_outer = .95,
                          point_est = "mean") 
mcmc_labels(pl, "Intercept","HIV", "Symptoms", "Motor", "Nerve", "PulTB", "MilTB", 
            "Age", "TBDays", "sqrt(BldGlu)", "sqrt(CsfGlu)", "log2(CsfLymp)",
            "log2(CsfPro)", "log2(CsfLac)", "sqrt(GluRatio)", "GCS")

library(ggplot2)
Predict.Smear = cbind(theta=m2c_result$summary("theta")$median, pred=m2c_result$summary("p_Smear")$median, obs=m2c_input$Y_Smear)
Predict.Smear = Predict.Smear %>% 
  as.data.frame() %>% 
  mutate(ranges = cut(pred, seq(0, 1, .1)),
         decile = cut(pred, quantile(pred, seq(0,1,.1))))
Predict.summary.Smear = Predict.Smear %>% group_by(ranges) %>% summarise(obs = sum(obs)/n(), n=n()) %>% mutate(avg = seq(.05, .95, .1))
Predict.summary.Smear = Predict.Smear %>% group_by(decile) %>% summarise(obs = sum(as.numeric(obs)-1)/n(), n=n()) %>%
  mutate(avg = c(0.006611082,0.017512770,0.037648303, 0.068974318, 0.111360510, 0.164783250, 0.220714850, 0.305990050, 0.468795400, 0.762204700, 0.975592250))
ggplot(Predict.summary.Smear) + geom_line(aes(x=avg,y=avg), linetype=2) + geom_line(aes(x=avg, y=obs))

Predict.Xpert = cbind(theta=m2c_result$summary("theta")$mean, pred=m2c_result$summary("p_Xpert")$median, obs=m2c_input$Y_Xpert)
Predict.Xpert = Predict.Xpert %>% 
  as.data.frame() %>% 
  mutate(ranges = cut(pred, seq(0, 1, .1)))
Predict.summary.Xpert = Predict.Xpert %>% group_by(ranges) %>% summarise(obs = sum(obs)/n(), n=n()) %>% mutate(avg = seq(.05, .85, .1))
ggplot(Predict.summary.Xpert) + geom_line(aes(x=avg,y=avg), linetype=2) + geom_line(aes(x=avg, y=obs))

