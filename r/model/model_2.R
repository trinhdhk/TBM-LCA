# Model with random effect
# Author trinhdhk
# Version: 0.1.2005

library(rstan)
model_2 <- stan_model('stan/model_2.stan')
result_2 <- sampling(model_2, data=model_input, chain=5, iter=30000, seed=128, warmup=5000, control = list(max_treedepth = 18, adapt_delta=.99))
extract_2 <- extract(result_2)
# plot((extract_1$theta), extract_1$a[,1])
saveRDS(result_2, file='outputs/result2.RDS')

