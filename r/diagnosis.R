library(bayesplot)
posterior_2 <- as.array(result_2)
np_2 <- nuts_params(result_2)
mcmc_parcoord(posterior_2, np = np_2, pars=dplyr::vars(-lp__))
mcmc_pairs(posterior_2, np = np_2, pars = c("z_Mgit[1]", 'z_Mgit[2]', 'a[8]'),
           off_diag_args = list(size = 0.75))
