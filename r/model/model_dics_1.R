# Model with dichomous CSF lab
# Author trinhdhk
# Version: 0.1.2005

library(rstan)
model_disc_1 <- stan_model('stan/model_disc_1.stan')
result_disc_1 <- sampling(model_disc_1, data=model_input_cont, chain=6, iter=30000, seed=128, warmup=5000, control = list(max_treedepth = 18, adapt_delta=.99))
saveRDS(result_disc_1, file='outputs/result_disc_1.RDS')
rm(model_disc_1)