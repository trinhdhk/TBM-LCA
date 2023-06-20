library(ggplot2)
# Script to calculate the relative TBM risk after test result
# Author: Trinh Dong
# Email: trinhdhk@oucru.org
##############################################################

#m3=readRDS('path/to/fit_object_with_r1k1.RDS')
thetas=rstan::extract(m3$outputs, pars=c('theta','Y0_theta', 'Y00_theta', 'Y000_theta'))
S = list(
  a = thetas$Y0_theta[,1,] / thetas$theta,
  b = thetas$Y0_theta[,3,] / thetas$theta,
  c = thetas$Y00_theta[,1,] / thetas$theta,
  d = thetas$Y00_theta[,3,] / thetas$theta,
  e = thetas$Y000_theta / thetas$theta,
  hiv = rstan::extract(m3$outputs, pars='X')$X[,,1]
)

saveRDS(S, 'export/post_test_fraction.RDS')
