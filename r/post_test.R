library(ggplot2)

library(future)
m3=m3t
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

# theta_dt = data.frame(
#   usubjid = data_19EI$USUBJID,
#   hiv = data_19EI$hiv_stat %in% TRUE
# )
# 
# cli::cli_progress_bar ("adding stuff", total=6400, clear=TRUE)
# # future::plan(future::multiprocess, workers=19)
# 
# 
# for (i in seq_len(6400)){
#   # future({
#     no_test <- glue::glue("no_test.{i}")
#     smear_neg <- glue::glue("smear_neg.{i}")
#     mgit_neg <- glue::glue("mgit_neg.{i}")
#     xpert_neg <- glue::glue("xpert_neg.{i}")
#     smear_mgit_neg <- glue::glue("smear_mgit_neg.{i}")
#     smear_xpert_neg <- glue::glue("smear_xpert_neg.{i}")
#     mgit_xpert_neg <- glue::glue("mgit_xpert_neg.{i}")
#     all_neg <- glue::glue("all_neg.{i}")
#     
#     theta_dt[[no_test]] = thetas$theta[i,]
#     theta_dt[[smear_neg]] = thetas$Y0_theta[i,1,]
#     theta_dt[[mgit_neg]] = thetas$Y0_theta[i,2,]
#     theta_dt[[xpert_neg]] = thetas$Y0_theta[i,3,]
#     theta_dt[[smear_mgit_neg]] = thetas$Y00_theta[i,1,]
#     theta_dt[[smear_xpert_neg]] = thetas$Y00_theta[i,2,]
#     theta_dt[[mgit_xpert_neg]] = thetas$Y00_theta[i,3,]
#     theta_dt[[all_neg]] = thetas$Y000_theta[i,]
#   # })
#   cli::cli_progress_update()
# }
# 
# theta_dt = data.frame(
#   usubjid = data_19EI$USUBJID,
#   hiv = data_19EI$hiv_stat %in% TRUE,
#   no_test = apply(thetas$theta,2,mean),
#   smear_neg = apply(thetas$Y0_theta[,1,],2,mean),
#   mgit_neg = apply(thetas$Y0_theta[,2,],2,mean),
#   xpert_neg = apply(thetas$Y0_theta[,2,],2,mean),
#   smear_mgit_neg = apply(thetas$Y00_theta[,1,],2,mean),
#   smear_xpert_neg = apply(thetas$Y00_theta[,2,],2,mean),
#   mgit_xpert_neg = apply(thetas$Y00_theta[,3,],2,mean),
#   all_neg = apply(thetas$Y000_theta,2,mean)
# ) |>
#   tidyr::pivot_longer(
#     cols = -c(usubjid, hiv),
#     # names_pattern = "(.*)\\.(\\d*)" ,
#     names_to = c('time'),
#     values_to = 'value'
#   ) |> 
#   mutate(time = factor(time, 
#                        levels=c('no_test','smear_neg','mgit_neg', 'xpert_neg', 'smear_mgit_neg', 
#                                 "smear_xpert_neg", "mgit_xpert_neg", 'all_neg'))) 
# 
# theta_dt_long = 
#   theta_dt |>
#   tidyr::pivot_longer(
#     cols = -c(usubjid, hiv),
#     names_pattern = "(.*)\\.(\\d*)" ,
#     names_to = c('time', 'chain'),
#     values_to = 'value'
#   ) |> 
#   mutate(time = factor(time, 
#                        levels=c('no_test','smear_neg','mgit_neg', 'xpert_neg', 'smear_mgit_neg', 
#                                 "smear_xpert_neg", "mgit_xpert_neg", 'all_neg')),
#          chain = as.numeric(chain)) 
# 
# ggplot(data=theta_dt_long, aes(x=time,y=value,color=hiv)) +
#   geom_boxplot()
# 
# fit=list()
# fit$fit2 = lme4::lmer(log2(value)~time*hiv + (1|usubjid), theta_dt)
# fit$fit1= lme4::lmer(log2(value)~time +hiv+ (1|usubjid), theta_dt)
# # fit = lm(log2(value)~time+hiv, data=theta_dt)
# saveRDS(fit, 'export/post_test.RDS')
