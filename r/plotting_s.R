
# s_m3_t0= rstan::sflist2stanfit(s_m3_t0_2$outputs)
# s_m3 = s_m3$outputs |> rstan::sflist2stanfit()
a = rstan::extract(s_m3, pars=c('a0', 'a'))
a = cbind(a$a0, a$a)
colnames(a) = c('a0', paste0('a[',1:(ncol(a)-1),']'))
scale = purrr::transpose(scale_Xc)
a[,10] = -a[,10] / scale$`scaled:scale`$gcs
a[,9] = a[,9] /scale$`scaled:scale`$id

a_plot_s2 = td.misc::mcmc_intervals_multi(list(a),
                                     pars = c('a0',
                                              paste0('a[',c(1:5, 9, 6:8, 10:11),']')),
                                     multi_point_est = FALSE,
                                     point_est = c('median'),
                                     point_size = 2.5,
                                     prob_outer = .95) |>
  td.misc::change_ylabs(
    'Intercept (centred)',
    'HIV +',
    'TB-suggested symptoms',
    'Focal neurological deficit',
    'Cranial nerve palsy',
    'Past noticed TB contact',
    'Glasgow coma score',
    'Pulmonary TB/X-ray',
    'Miliary TB/X-Ray',
    # '*log<sub>2</sup>* (Age)',
    '*log<sub>2</sup>* (Days from onset)',
    'Headache',
    'Pyschosis'
  ) +  
  geom_vline(aes(xintercept=0), color=grey(.5), alpha=.5) +
  theme_bw() +
  scale_color_discrete(type=RColorBrewer::brewer.pal(3,'Set1')) +
  ggplot2::theme(axis.text.y = ggtext::element_markdown(), legend.position = 'none') +
  xlab('') + ylab('')

# a_plot_s2$data$m = bayesplot::mcmc_intervals_data(a,
#                                                   pars = c('a0',
#                                                            paste0('a[',c(1:5, 9, 6:8, 10:11),']')),
#                                                   transformations = plogis, 
#                                                   prob_outer=.95, point_est='mean') %>%
#   mutate(across(c(ll:hh), ~ qlogis(.x))) %>% .$m
a_plot_s = a_plot_s2 +
  scale_x_continuous(
    name = 'TBM odd ratio',
    breaks = log(c(.001, .01,.05, .1, .25, .5, 1, 2, 5, 10, 100)),
    labels = c(.001, .01,.05, .1, .25, .5, 1, 2, 5, 10, 100))

intercept_translation = unlist(scale$`scaled:center`[c('id', 'gcs')]) * cbind(a[,9], -a[,10])
intercept_translation[,2] = intercept_translation[,2] - (15 * (-a[,10]))
a[,1] = a[,1] - apply(intercept_translation,1,sum)
a_plot_s3 = td.misc::mcmc_intervals_multi(list(a),
                                          pars = c('a0',
                                                   paste0('a[',c(1:5, 9, 6:8, 10:11),']')),
                                          point_est = 'median',
                                          point_size = 2.5,
                                          prob_outer = .95) |>
  td.misc::change_ylabs(
    'Intercept',
    'HIV +',
    'TB-suggested symptoms',
    'Focal neurological deficit',
    'Cranial nerve palsy',
    'Past noticed TB contact',
    'Glasgow coma score',
    'Pulmonary TB/X-ray',
    'Miliary TB/X-Ray',
    # '*log<sub>2</sup>* (Age)',
    '*log<sub>2</sup>* (Days from onset)',
    'Headache',
    'Pyschosis'
  ) +  
  geom_vline(aes(xintercept=0), color=grey(.5), alpha=.5) +
  theme_bw() +
  scale_color_discrete(type=RColorBrewer::brewer.pal(3,'Set1')) +
  ggplot2::theme(axis.text.y = ggtext::element_markdown(), legend.position = 'none') +
  xlab('') + ylab('') +
  scale_x_continuous(
    name = 'TBM odd ratio',
    breaks = log(c(.001, .01,.05, .1, .25, .5, 1, 2, 5, 10, 100)),
    labels = c(.001, .01,.05, .1, .25, .5, 1, 2, 5, 10, 100))


# a_plot_s3$data$m = bayesplot::mcmc_intervals_data(a,
#                                                   pars = c('a0',
#                                                            paste0('a[',c(1:5, 9, 6:8, 10:11),']')),
#                                                   transformations = plogis, 
#                                                   prob_outer=.95, point_est='mean') %>%
#   mutate(across(c(ll:hh), ~ qlogis(.x))) %>% .$m

saveRDS(list(a_plot_s=a_plot_s,a_plot_s2=a_plot_s2,a_plot_s3=a_plot_s3), 'export/a_plot_s.RDS')
# save(a_plot_s, file='export/a_plot_s.Rdata')

