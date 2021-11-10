# Master execution
# Author: trinhdhk
# Version 0.1.2005

# Load libs
rm(list = ls())
library(data.table)
library(magrittr)
library(rstan)
library(bayesplot)
options(mc.cores = 4)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
rstan_options(auto_write = TRUE)
source('r/fn_def.R')

# Load data
source('r/load_data.R')

# Load model
## model 1: without random effect
model_disc_1 <- load_model('disc_1', data=model_input_disc, chain = 4, iter=40000, seed=2906, warmup=10000, 
                           control = list(max_treedepth = 18, adapt_delta=.99))

model_disc_1b <- load_model('disc_1b', data=model_input_disc, chain = 4, iter=30000, seed=1208, warmup=5000, 
                           control = list(max_treedepth = 18, adapt_delta=.99))

model_cont_1 <- load_model('cont_1', data=model_input_cont, chain = 4, iter=40000, seed=1208, warmup=10000, 
                           control = list(max_treedepth = 18, adapt_delta=.99))

model_cont_1b <- load_model('cont_1b', data=model_input_cont, chain = 4, iter=40000, seed=1208, warmup=10000, 
                           control = list(max_treedepth = 18, adapt_delta=.99))
## model 2: with random effect ~ normal(0,5)

model_disc_2 <- load_model('disc_2', data=model_input_disc, chain = 4, iter=40000, seed=1208, warmup=10000, 
                           control = list(max_treedepth = 18, adapt_delta=.99))

model_disc_2h <- load_model('disc_2h', data=model_input_disc, chain = 4, iter=50000, seed=1208, warmup=20000, 
                           control = list(max_treedepth = 18, adapt_delta=.99))

model_disc_2h2 <- load_model('disc_2h2', data=model_input_disc, chain = 4, iter=50000, seed=1208, warmup=20000, 
                            control = list(max_treedepth = 18, adapt_delta=.99))

model_cont_2 <- load_model('cont_2', data=model_input_cont, chain = 4, iter=60000, seed=128, warmup=30000, 
                           control = list(max_treedepth = 18, adapt_delta=.99))

## model 3: with random effect ~ normal(0,s)
model_cont_3 <- load_model('cont_3', data=model_input_cont, chain = 4, iter=40000, seed=1208, warmup=10000, 
                           control = list(max_treedepth = 18, adapt_delta=.99)) #diverged

model_disc_3h <- load_model('disc_3h', data=model_input_disc, chain = 4, iter=60000, seed=1208, warmup=30000, 
                            control = list(max_treedepth = 18, adapt_delta=.99))

model_disc_3hc <- load_model('disc_3hc', data=model_input_disc, chain = 4, iter=60000, seed=1208, warmup=30000, 
                             control = list(max_treedepth = 18, adapt_delta=.995), init_r = .5)

model_disc_3h2 <- load_model('disc_3h2', data=model_input_disc, chain = 3, iter=60000, seed=1708, warmup=30000, 
                             control = list(max_treedepth = 18, adapt_delta=.999), init_r = .5)

library(future)
plan('multiprocess')


model_disc_3hc <- future({
  options(mc.cores = 4)
  Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
  
})

model_disc_3h2c <- future({
  load_model('disc_3h2c', data=model_input_disc, chain = 4, iter=60000, seed=1208, warmup=30000, 
             control = list(max_treedepth = 18, adapt_delta=.99))
})
