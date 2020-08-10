library(nimble)
library(data.table)
library(magrittr)
modelcode_19EI <- nimbleCode({
  a0 ~ dnorm(0,1)
  b_age ~ dnorm(0,1)
  b_HIV ~ dnorm(0,1)
  b_mil ~ dnorm(0,1)
  for (i in 1:nX) a[i] ~ dnorm(0,1)
  for (i in 1:N){
    X[i,4] ~ dpois(d00 + d01*X[i,1] + d02*Xaux[i,1]+d03*Xaux[i,2]+d04*Xaux[i,3]+d05*Xaux[i,4])
    X[i,1] ~ dbern(d10 + d11*Xaux[i,1]+d13*Xaux[i,3])
    Xaux[i,1] ~ dbern(d20 + d21*X[i,1]+d22*Xaux[i,2]+d23*Xaux[i,3])
    Xaux[i,2] ~ dbern(d30 + d31*X[i,1]+d32*Xaux[i,1])
    Xaux[i,3] ~ dbern(d40 + d41*X[i,1]+d43*Xaux[i,2])
  }
  
  X[1:N,3] <- Xaux[1:N,1] | Xaux[1:N,2] | Xaux[1:N,3]
  X[1:N,4] <- log(X[1:N,4])
  
  # for (i in 1:N){
  #   for (j in 1:2){
  #     X[i,j] <- X[i,j]
  #   }
  #   X[i,3] <- Xaux[i,1] | Xaux[i,2] | Xaux[i,3]
  #   X[i,4] <- log(X[i,4])
  # 
  #   for (j in 5:nX){
  #     X[i,j] <- X[i,j]
  #   }
  # }
  
  
  # X[1:N,1:2] <- Xprime[1:N,1:2]
  # X[1:N,3] <- Xaux[1:N,1] | Xaux[1:N,2] | Xaux[1:N,3]
  # X[1:N,4] <- log(Xprime[1:N,4])
  # X[1:N,5:nX] <- Xprime[1:N,4:(nX-1)]
  
  z_Xpert[1] ~ dnorm(probit(.005), .7)
  z_Mgit[1] ~ dnorm(-3.023, .89)
  z_Smear[1] ~ dnorm(-3.023, .89)
  
  z_Xpert[2] ~ dnorm(probit(.593), .117);
  z_Mgit[2]  ~ dnorm(probit(.665), .217);
  z_Smear[2] ~ dnorm(probit(.786), .405);
  
  for (i in 1:3) b[i] ~ T(dt(0, 1, 1),0,)

  theta[] <- ilogit(rep(a0, N) + inprod(X[1:N, 1:nX], a[1:nX]))
  bac_load[] <- RE[1:N] + b_HIV * X[1:N,1] + b_age*X[1:N,2] + b_mil*X[1:N, nX]
  
  for (i in 1:N){
    RE[i] ~ T(dnorm(0,1),0,)
    C[i] ~ dbern(theta[i])
    
    p_Xpert[i] <- ilogit(z_Xpert[1])*(1-C[i]) + ilogit(z_Xpert[2] + b[3]*bac_load[i])*C[i]
    p_Mgit[i]  <- ilogit(z_Mgit[1])*(1-C[i])  + ilogit(z_Mgit[2] +  b[2]*bac_load[i])*C[i]
    p_Smear[i] <- ilogit(z_Smear[1])*(1-C[i]) + ilogit(z_Smear[2]+ b[1]*bac_load[i])*C[i]
  
    Y_Xpert[i] ~ dbern(p_Xpert[i]);
    Y_Mgit[i]  ~ dbern(p_Mgit[i]);
    Y_Smear[i] ~ dbern(p_Smear[i]);
  }
})

data_19EI <- readRDS('data/cleaned/data_19EI.RDS')
data_19EI <- na.omit(data_19EI, cols = c('hiv_stat', 'age', 'csf_smear', 'csf_mgit', 'csf_xpert',
                                         'clin_illness_day', 'clin_symptoms', 'clin_gcs', 
                                         'clin_nerve_palsy', 'clin_motor_palsy', 
                                         'csf_clear', 'csf_protein', 'csf_lympho', 'glucose_ratio', 'csf_lactate',
                                         'xray_miliary_tb', 'xray_pul_tb'))

X <- data_19EI %$% as.matrix(cbind(
  hiv_stat,
  (age-mean(age))/10, 
  clin_symptoms,
  # log2(clin_illness_day),
  clin_illness_day,
  clin_nerve_palsy,
  clin_motor_palsy,
  ISDIABETE,
  GLASCOW,
  glucose_ratio,
  BLDGLU,
  log2(csf_lympho+1),
  log2(csf_protein),
  log2(csf_lactate),
  xray_pul_tb,
  xray_miliary_tb
)) %>% unname

Xaux <- data_19EI %$% cbind(
  ISNSWEAT,
  ISWEIGHT,
  ISCOUGH,
  ISFEVER
)

const_19EI <-list(N = nrow(data_19EI), nX= ncol(X))
dt_19EI <- data_19EI %$% list(
  Y_Smear = as.integer(csf_smear),
  Y_Mgit = as.integer(csf_mgit),
  Y_Xpert = as.integer(csf_xpert),
  X=X,
  Xaux =Xaux
)

init_19EI <- function() list(a0 = runif(1,-1,1), a = runif(ncol(X), -1, 1),
                  b_HIV=runif(1,-1,1), b_age = runif(1,-1,1), b_mil = runif(1,-1,1), b = runif(3,0,2),
                  d01=runif(1,-1,1),d02=runif(1,-1,1),d03=runif(1,-1,1),d04=runif(1,-1,1),d05=runif(1,-1,1),
                  d11=runif(1,-1,1),d12=runif(1,-1,1),d13=runif(1,-1,1),d14=runif(1,-1,1),
                  d21=runif(1,-1,1),d22=runif(1,-1,1),d23=runif(1,-1,1),d24=runif(1,-1,1),
                  d31=runif(1,-1,1),d32=runif(1,-1,1),d33=runif(1,-1,1),d34=runif(1,-1,1),
                  d41=runif(1,-1,1),d42=runif(1,-1,1),d43=runif(1,-1,1),d44=runif(1,-1,1),
                  z_Smear = runif(2,-1,1),
                  z_Mgit=runif(2,-1,1),
                  z_Xpert=runif(2,-1,1),
                  RE = runif(nrow(data_19EI),0,2), 
                  C=rbinom(nrow(data_19EI),1,.5),
                  theta = runif(nrow(data_19EI), 0,1), bac_load = runif(nrow(data_19EI), 0, 2))

nimble_19EI <- nimbleModel(code=modelcode_19EI, constants = const_19EI, data=dt_19EI, init=init_19EI())
mcmcconf_19EI <- configureMCMC(nimble_19EI)
mcmc_19EI <- buildMCMC(mcmcconf_19EI)
model_19EI <- compileNimble(nimble_19EI, mcmc_19EI)
res_19EI <-  runMCMC(model_19EI$mcmc_19EI, niter=20000, nburnin = 8000, nchains=3)


dyesFoldFunction <- function(i){
  foldNodes_i <- c(paste0('Y_Smear[', i, ']'),paste0('Y_Mgit[', i, ']'),paste0('Y_Xpert[', i, ']'))
  return(foldNodes_i)
}

RMSElossFunction <- function(simulatedDataValues, actualDataValues){
  dataLength <- length(simulatedDataValues) # simulatedDataValues is a vector
  SSE <- 0
  for(i in 1:dataLength){
    SSE <- SSE + (simulatedDataValues[i] - actualDataValues[i])^2
  }
  MSE <- SSE / dataLength
  RMSE <- sqrt(MSE)
  return(RMSE)
}

EI_MCMCconfiguration <- configureMCMC(nimble_19EI)

crossValOutput <- runCrossValidate(MCMCconfiguration = EI_MCMCconfiguration,
                                   k = 5,
                                   foldFunction = 'random',
                                   lossFunction = RMSElossFunction,
                                   nCores = 15,
                                   MCMCcontrol = list(niter = 20000,
                                                      nburnin = 10000))

res_19EI <- nimbleMCMC(code=modelcode_19EI,constants = const_19EI, data=dt_19EI, init=init_19EI,
                       monitors = c('C', 'a', 'a0', 'b_HIV', 'b_age', 'b_mil', 'b', 'z_Smear', 'z_Mgit', 'z_Xpert', 
                                    'theta'), niter=50000, nburnin = 20000, nchains=3, setSeed=2908, progressBar = T,
                       )

res_19EI <- nimbleMCMC(code=modelcode_19EI,constants = const_19EI, data=dt_19EI, 
                       init=init_19EI, niter=20000, nburnin = 10000, nchains=1, setSeed=2908, progressBar = T,
                       monitors = 'logProbsY_Mgit')
nimbleOptions("showCompilerOutput")
com_19EI <- compileNimble()