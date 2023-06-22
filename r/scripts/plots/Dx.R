# Plot the Stan diagnostic metrics
# m3 is the full model, s is the simplified model

mcse = rstan::stan_mcse(m3$outputs)
rhat = rstan::stan_rhat(m3$outputs)
ess = rstan::stan_ess(m3$outputs)
saveRDS(list(mcse=ggplotGrob(mcse), rhat=ggplotGrob(rhat), ess=ggplotGrob(ess)), 'outputs/metrics/m3_diag.RDS')

mcse = rstan::stan_mcse(s)
rhat = rstan::stan_rhat(s)
ess = rstan::stan_ess(s)
saveRDS(list(mcse=ggplotGrob(mcse), rhat=ggplotGrob(rhat), ess=ggplotGrob(ess)), 'outputs/metrics/s_diag.RDS')
