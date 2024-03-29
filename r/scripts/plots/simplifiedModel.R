#Script to generate plots of coefficient for simplified model
# Author: Trinh Dong
# Email: trinhdhk@oucru.org
library(tidyverse)
library(ggplot2)
load('data/cleaned/data_input.Rdata')
# sm3 = s_m3$outputs |> rstan::sflist2stanfit()
a = rstan::extract(sm3, pars=c('a0', 'a'))
a = cbind(a$a0, a$a)
colnames(a) = c('a0', paste0('a[',1:(ncol(a)-1),']'))
scale = purrr::transpose(scale_Xc)
a[,ncol(a)] = -a[,ncol(a)] / scale$`scaled:scale`$gcs #10
a[,ncol(a)-1] = a[,ncol(a)-1] /scale$`scaled:scale`$id #9

a_plot_s2 = td.misc::mcmc_intervals_multi(list(a),
                                     pars = c('a0',
                                              paste0('a[',c(1:5,ncol(a)-1, 6:7, ncol(a)-2,8,9,10,11),']')),
                                              # paste0('a[',c(1:5, 7, 6, 8:9),']')),
                                     multi_point_est = FALSE,
                                     point_est = c('median'),
                                     point_size = 2.5,
                                     prob_outer = .95) |>
  td.misc::change_ylabs(
    'Intercept',
    'HIV +',
    'TB-suggestive symptoms',
    'Focal neurological deficit',
    'Cranial nerve palsy',
    'Past noticed TB contact',
    'Glasgow coma score',
    'Pulmonary TB/X-ray',
    'Miliary TB/X-ray',
    # '*log<sub>2</sup>* (Age)',
    '*log<sub>2</sup>* (Symptom duration, days)',
    'Headache',
    'Fever',
    'Neck stiffness',
    'Pyschosis'
  ) +  
  geom_vline(aes(xintercept=0), color=grey(.5), alpha=.5) +
  theme_bw() +
  scale_color_discrete(type=RColorBrewer::brewer.pal(3,'Set2')) +
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
    name = 'TBM odds ratio',
    breaks = log(c(.001, .01,.05, .1, .25, .5, 1, 2, 5, 10, 20)),
    labels = c(.001, .01,.05, .1, .25, .5, 1, 2, 5, 10, 20))

intercept_translation = unlist(scale$`scaled:center`[c('id', 'gcs')]) * cbind(a[,ncol(a)-1], -a[,ncol(a)])
intercept_translation[,2] = intercept_translation[,2] - (15 * (-a[,ncol(a)]))
a[,1] = a[,1] - apply(intercept_translation,1,sum)
a_plot_s3 = td.misc::mcmc_intervals_multi(list(a),
                                          pars = c('a0',
                                                   # paste0('a[',c(1:5, 9, 6:8, 10:11),']')),
                                                   paste0('a[',c(1:5,ncol(a)-1, 6:7, ncol(a)-2,8,9,10,11),']')),
                                                   # paste0('a[',c(1:5, 7, 6, 8:9),']')),
                                          point_est = 'median',
                                          point_size = 2.5,
                                          prob_outer = .95) |>
  td.misc::change_ylabs(
    'Intercept',
    'HIV +',
    'TB-suggestive symptoms',
    'Focal neurological deficit',
    'Cranial nerve palsy',
    'Past noticed TB contact',
    'Glasgow coma score',
    'Pulmonary TB/X-ray',
    'Miliary TB/X-Ray',
    # '*log<sub>2</sup>* (Age)',
    '*log<sub>2</sup>* (Symptom duration, days)',
    'Headache',
    'Fever',
    'Neck stiffness',
    'Pyschosis'
  ) +  
  geom_vline(aes(xintercept=0), color=grey(.5), alpha=.5) +
  theme_bw() +
  scale_color_discrete(type=RColorBrewer::brewer.pal(3,'Set2')) +
  ggplot2::theme(axis.text.y = ggtext::element_markdown(), legend.position = 'none') +
  xlab('') + ylab('') +
  scale_x_continuous(
    name = 'TBM odds ratio',
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

