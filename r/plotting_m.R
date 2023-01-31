library(ggplot2)
a = rstan::extract(m3$outputs, pars=c('a0', 'a'))
a_mX = rstan::extract(m3m$outputs, pars=c('a0', 'a'))
a_mXu = rstan::extract(m3mu$outputs, pars=c('a0', 'a'))

a$a[,18] = -a$a[,18]
a = cbind(a$a0, a$a)

a_mX$a[,18] = -a_mX$a[,18]
a_mX = cbind(a_mX$a0, a_mX$a)

a_mXu$a[,18] = -a_mXu$a[,18]
a_mXu = cbind(a_mXu$a0, a_mXu$a)

colnames(a) <- colnames(a_mX) <- colnames(a_mXu) <-  c('a0', paste0('a[',1:(ncol(a)-1),']'))
labs = c(
  'Intercept', 
  'HIV +',
  'TB-suggestive symptoms',
  'Focal neurological deficit',
  'Cranial nerve palsy',
  'Past noticed TB contact',
  'Glasgow coma score',
  'Pulmonary TB/X-ray',
  'Miliary TB/X-Ray',
  '*log<sub>2</sub>* (Symptom duration, days)',
  '*log<sub>2</sub>* (Paired blood glucose)',
  '*log<sub>2</sub>* (CSF glucose)',
  '*log<sub>2</sub>* (CSF protein)',
  '*log<sub>2</sub>* (CSF lactate)',
  '*log<sub>10</sub>* (CSF lymphocyte)',
  '*log<sub>10</sub>* (CSF WBC)',
  '*log<sub>10</sub>* (CSF WBC)<sup>2</sup>',
  'CSF eosinophil > 0',
  '*log<sub>10</sub>* (CSF eosinophil)',
  '*log<sub>10</sub>* (CSF RBC)',
  'Evidence of cryptococcus',
  'Positive CSF Gram stain'
) 

cutpoints = c(-20,-16,-8,-4, -2,-1, log10(.5), log10(.9) ,0, log10(10/9),log10(2),1,2,4)
lab.cutpoints = c(paste0('1e', c(-20,-16,-8,-4)), .01, .1,.5, .9,1, 1.1, 2, 10 ,100, '1e4')

a_plot_m <-
  (td.misc::mcmc_intervals_multi(
    fits = list('Selected' = a, `Incomplete Xpert` = a_mX, `Incomplete Xpert + Wide prior` = a_mXu),
    pars = c('a0',
             paste0('a[',c(1:5,18,6:7,11:13,15,16,14,17,21,10,19,20,8,9),']')),
    multi_point_est = FALSE,
    point_est = c('median'),
    point_size=2,
    transformations = function(x) ifelse(x==0, 0, sign(x) * sqrt(abs(x))),
    prob_outer = .95) +
     geom_vline(aes(xintercept=0), color=grey(.5), alpha=.5) +
     scale_color_discrete(type=RColorBrewer::brewer.pal(name='Set1', n=3)[c(2,1,3)]) +
     scale_x_continuous(name='TBM odds ratio', 
                        breaks = (function(x) ifelse(x==0, 0, sign(x) * sqrt(abs(x))))(log(10^cutpoints)),
                        labels=lab.cutpoints) +
     ylab('')+
     theme_bw() +
     theme(axis.text.y = ggtext::element_markdown(), legend.position = 'bottom')) |>
     # guides(color=guide_legend(title="Scale"))) %>%
     td.misc::change_ylabs(
       labs = labs,
       top_down = T
     ) 

# a_plot_m_data_m1 = bayesplot::mcmc_intervals_data(a,
#                                                   pars = c('a0',
#                                                            paste0('a[',c(1:5,18,6:7,11:13,15,16,14,17,21,10,19,20,8,9),']')),
#                                                   transformations = plogis, 
#                                                   prob_outer=.95, point_est='mean') %>%
#   mutate(across(c(ll:hh), ~ sign(qlogis(.x)) * sqrt(abs(qlogis(.x))))) %>% .$m

# a_plot_m_data_m2 = bayesplot::mcmc_intervals_data(a_m,
#                                                   pars = c('a0',
#                                                            paste0('a[',c(1:5,18,6:7,11:13,15,16,14,17,21,10,19,20,8,9),']')),
#                                                   transformations = plogis, 
#                                                   prob_outer=.95, point_est='mean') %>%
#   mutate(across(c(ll:hh), ~ sign(qlogis(.x)) * sqrt(abs(qlogis(.x))))) %>% .$m

# a_plot_m$data$m = c(a_plot_m_data_m1, a_plot_m_data_m2)

b = rstan::extract(m3$outputs, pars=c('b_HIV', 'b'))
b_mX = rstan::extract(m3m$outputs, pars=c('b_HIV', 'b'))
b_mXu = rstan::extract(m3mu$outputs, pars=c('b_HIV', 'b'))

b$b[,6]=-b$b[,6]
b = cbind(b$b_HIV, b$b)

b_mX$b[,6]=-b_mX$b[,6]
b_mX = cbind(b_mX$b_HIV, b_mX$b)

b_mXu$b[,6]=-b_mXu$b[,6]
b_mXu = cbind(b_mXu$b_HIV, b_mXu$b)

colnames(b) <- colnames(b_mX) <- colnames(b_mXu) <- c('b_HIV', paste0('b[',1:(ncol(b)-1), ']'))

b_plot_m = td.misc::mcmc_intervals_multi(list('Selected' = b, `Incomplete Xpert` = b_mX, `Incomplete Xpert + Wide prior` = b_mXu),
                                       pars = c('b_HIV',  paste0('b[',c(1,6,4,5,2,3), ']')),
                                       point_est = 'mean',
                                       point_size=2,
                                       prob_outer = .95) |>
  td.misc::change_ylabs(
    'HIV +',
    'GCS',
    '*log<sub>2</sub>* (CSF glucose)',
    '*log<sub>2</sub>* (CSF protein)',
    '*log<sub>2</sub>* (CSF lactate)',
    '*log<sub>10</sub>* (CSF lymphocyte)',
    '*log<sub>10</sub>* (CSF WBC)',
    top_down = T
  ) +
  geom_vline(aes(xintercept=0), color=grey(.5), alpha=.5) +
  scale_color_discrete(type=RColorBrewer::brewer.pal(name='Set1', n=3)[c(2,1,3)]) + 
  guides(color=guide_legend(title="Model")) +
  ylab('') +
  scale_x_continuous(name='Standardised bacillary burden', breaks=c(-3,-2,-1,0,1,2,3))+
  theme_bw() + theme(axis.text.y = ggtext::element_markdown(family=NULL), legend.position = 'bottom')

z = rstan::extract(m3$outputs, pars=c('z_Smear', 'z_Mgit', 'z_Xpert'))
z_mX = rstan::extract(m3m$outputs, pars=c('z_Smear', 'z_Mgit', 'z_Xpert', 'z_obs'))
z_mXu= rstan::extract(m3mu$outputs, pars=c('z_Smear', 'z_Mgit', 'z_Xpert', 'z_obs'))
z = do.call(cbind, z)
z_mX = do.call(cbind, z_mX)
z_mXu = do.call(cbind, z_mXu)

colnames(z) <- paste0('z_',rep(c('Smear', 'Mgit', 'Xpert'),each=2), rep(c('[1]', '[2]'), 3))
colnames(z_mX) <- colnames(z_mXu) <- paste0('z_',rep(c('Smear', 'Mgit', 'Xpert', 'observe'),each=2), rep(c('[1]', '[2]'), 4))

z_plot_logit_m = (td.misc::mcmc_intervals_multi(list(Original = z, `Incomplete Xpert` = z_mX, `Incomplete Xpert + Wide prior` = z_mXu),
                              regex_pars = '^z',
                              point_est = c('median'),
                              point_size=2,
                              # multi_point_est = TRUE,
                              # point_size = 2.5,
                              prob_outer = .95) +
                    ylab('') +
                    geom_vline(aes(xintercept=0), color=grey(.5), alpha=.5) +
                    scale_color_discrete(type=RColorBrewer::brewer.pal(name='Set1', n=3)[c(2,1,3)]) +
                    scale_x_continuous(name='',
                                      breaks = qlogis(c(.00001,.0001,.001,.01, .05, .25,  .5, .75, .95, .99, .999, .9999)),
                                      labels = c('1e-4','1e-3','1e-2','.01', '.05', '.25', '.5', '.75', '.95', '.99', '.999', '.9999')) +
                    theme_bw()) |>
  td.misc::change_ylabs(
    c(
      "Response rate/TBM",
      "Response rate/non-TBM",
      "TPR<sub>Xpert</sub>",
      "FPR<sub>Xpert</sup>",
      "TPR<sub>MGIT</sub>",
      "FPR<sub>MGIT</sub>",
      "TPR<sub>ZN-Smear</sub>",
      "FPR<sub>ZN-Smear</sub>"
    ) |> rev(),
    top_down = T
  ) + theme(axis.text.y = ggtext::element_markdown(family=NULL), legend.position = 'bottom')

# z_plot_m_data_m1 = bayesplot::mcmc_intervals_data(z,
#                                                   regex_pars = '^z',
#                                                   transformations = plogis,
#                                                   prob_outer=.95,
#                                                   multi_point_est = TRUE,
#                                                   point_est = 'mean') %>%
#   mutate(across(c(ll:hh), ~ qlogis(.x))) %>% .$m

# z_plot_m_data_m2 = bayesplot::mcmc_interval
                                                  # regex_pars = '^z',
                                                  # transformations = plogis,
                                                  # prob_outer=.95,
                                                  # multi_point_est = TRUE,
                                                  # point_est = 'mean')%>%
  # mutate(across(c(ll:hh), ~ qlogis(.x))) %>% .$m

# z_plot_logit_m$data$m = c(z_plot_m_data_m1, z_plot_m_data_m2)

saveRDS(list(a_plot = a_plot_m, b_plot = b_plot_m, z_plot = z_plot_logit_m), 'export/m3m_plot.RDS')
