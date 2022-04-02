#!/usr/bin/env Rscript

library(argparser)
library(magrittr, include.only = "%>%")
misc <- new.env(hash = FALSE)
source("r/include/functions.R", local = misc)
rstan::rstan_options(javascript = FALSE, auto_write = TRUE)

args <- 
  arg_parser("TBM-LCA missing simulator") %>%
  add_argument("scenario", help = "scenario: 1: no prior info of sen, spc, and obs; 2: Good est but wide confidence intervals; 3: Good est of sen, spc, but bad for obs_rate; 4: Bad est of sen, spc, but good for obs_rate; 5: All good") %>%
  add_argument("--num-simulation", help="Number of simulation", short="-n") %>%
  add_argument("--cores", help="Number of cores", default = parallel::detectCores() -1, short="-c") %>%
  add_argument("--output-file", help="output-file", short="-o") %>%
  add_argument("--output-dir", help="output-dir", default="outputs", short="-d")

generate_Y = \(C, probs){
  Cs = sort(unique(C))
  sapply(C, \(c) rbinom(1,1, probs[Cs==c]))
}

create_data <- function(good_prior){
  # Sample size
  N = 600 
  # Create data with 3 predictors, X, X2, and X3
  X = rnorm(N, 0, 1)
  X2 = rnorm(N, 0, 1)
  X3 = rbinom(N,1, .3)
  
  # Create latent class
  probs = plogis(3*X+X2+5*X3-1)
  C = sapply(probs, \(p) rbinom(1,1,p))
  
  # Manifest variables
  Y1 = generate_Y(C, c(.1,.30))
  Y2 = generate_Y(C, c(.01,.50))
  Y3 = generate_Y(C, c(0.05, .8))
  
  # Simulate class-aware missing data for Y3. If obs==1, Y is observed
  obs_rate = c(0.1, .95)
  obs = generate_Y(C, obs_rate)

  list(N=N,nX=3, X=cbind(X,X2,X3),C=C, Y1=Y1, Y2=Y2, Y3=Y3, obs=obs, good_prior=good_prior)
}


sampler = misc$my_stan_model('stan/simulation.stan')

argparser <- parse_args(args)
results <- new.env(parent=emptyenv())
results$.META <- argparser

with(argparser, {
  output_file <- file.path(output_dir, paste0(output_file,".RDS"))
  
  good_prior <- switch(scenario, "1" = c(0, 0), "2" = c(1, 1), "3" = c(1, 0), "4" = c(0, 1), "5" = c(2, 2))
  num_simulation <- as.integer(num_simulation)
  results$datas <- lapply(seq_len(num_simulation),
                  function(.) create_data(good_prior))
  
  future::plan(future::multicore, workers=cores, gc=TRUE)
  
  cli::cli_alert("Running simulations")
  progressr::handlers(progressr::handler_progress)
  progressr::with_progress({
    p <- progressr::progressor(steps = num_simulation)
    
    results$outputs <- 
      furrr::future_map(seq_len(num_simulation),
                                  function(i){
                                    # cli::cli_progress_update(.envir = results)
                                    sf <- rstan::sampling(
                                      sampler,
                                      data = results$datas[[i]],
                                      chains=2,
                                      cores=1,
                                      iter=6000,
                                      warmup=2000,
                                      init_r=1,
                                      control=list(adapt_delta=.69, max_treedepth=12)
                                    ) 
                                    p()
                                    sf
                                  }, .options = furrr::furrr_options(seed=TRUE))
  }, enable=TRUE)
  

  cli::cli_alert("Saving simulations")
  
  saveRDS(results, output_file)
  cli::cli_alert_success("Sampling completed!")
})
