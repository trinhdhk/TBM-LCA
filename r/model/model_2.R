# Model with random effect
# Author trinhdhk
# Version: 0.1.2005

library(rstan)
model_2 <- stan_model('stan/model_2.stan')
result_2 <- sampling(model_2, data=model_input, chain=6, iter=18000, seed=128, warmup=1000, control = list(max_treedepth = 18, adapt_delta=.99))
extract_2 <- extract(result_2)
# plot((extract_1$theta), extract_1$a[,1])
saveRDS(result_2, file='outputs/result2.RDS')
# library(bayesplot)
# posterior_2 <- as.array(result_2)
# np_2 <- nuts_params(result_2)
# mcmc_parcoord(posterior_2, np = np_2, pars=dplyr::vars(-lp__))
# mcmc_pairs(posterior_2, np = np_2, pars = c("z_Mgit[1]", 'z_Img[1]'),
#            off_diag_args = list(size = 0.75))
