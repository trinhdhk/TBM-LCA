# Model without random effect
# Author trinhdhk
# Version: 0.1.2004

model_1 <- stan_model('stan/model_1.stan')
result_1 <- sampling(model_1, data=model_input, chain=6, iter=8000, seed=128, warmup=1000, control = list(max_treedepth = 18, adapt_delta=.97))
extract_1 <- extract(result_1)
# plot((extract_1$theta), extract_1$a[,1])

library(bayesplot)
posterior_1 <- as.array(result_1)
np_1 <- nuts_params(result_1)
mcmc_parcoord(posterior_1, np = np_1, pars=dplyr::vars(-lp__))
mcmc_pairs(posterior_1, np = np_1, pars = c("z_Mgit[1]", 'z_Img[1]'),
           off_diag_args = list(size = 0.75))
