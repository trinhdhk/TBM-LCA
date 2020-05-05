# Master execution
# Author: trinhdhk
# Version 0.1.2005

# Load libs
rm(list = ls())
library(data.table)
library(magrittr)
library(rstan)
library(bayesplot)
options(mc.cores = 6)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
rstan_options(auto_write = TRUE)
source('r/fn_def.R')

# Load data
source('r/load_data.R')

# Load model
## model 1: without random effect
load_model(1)
## model 2: with random effect ~ normal(0,5)
load_model(2)