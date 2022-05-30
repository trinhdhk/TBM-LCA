library(ggplot2)
a_m = rstan::extract(m3mX$outputs, pars=c('a0', 'a'))
# a$a[,1:ncol(Xd)] = a$a[,1:ncol(Xd)] * 2
a_m$a[,18] = -a_m$a[,18]
a_m = cbind(a_m$a0, a_m$a) 
colnames(a_m) = c('a0', paste0('a[',1:(ncol(a)-1),']'))
labs = c(
  'Intercept (centred)', 
  'HIV +',
  'TB-suggested symptoms',
  'Focal neurological deficit +',
  'Cranial nerve palsy +',
  'Past noticed TB contact',
  'Glasgow coma score',
  'Pulmonary TB/X-ray',
  'Miliary TB/X-Ray',
  '*log<sub>2</sub>* (Days from onset)',
  '*log<sub>2</sub>* (Paired blood glucose)',
  '*log<sub>2</sub>* (CSF Glucose)',
  '*log<sub>2</sub>* (CSF Protein)',
  '*log<sub>2</sub>* (CSF Lactate)',
  '*log<sub>10</sub>* (CSF Lymphocyte)',
  '*log<sub>10</sub>* (CSF WBC)',
  '*log<sub>10</sub>* (CSF WBC)<sup>2</sup>',
  'CSF Eosinophil > 0',
  '*log<sub>10</sub>* (CSF Eosinophil)',
  '*log<sub>10</sub>* (CSF RBC)',
  'Cryptococcus +',
  'Gram stain +'
) 

cutpoints = c(-20,-16,-8,-4, -2,-1, log10(.5), log10(.9) ,0, log10(10/9),log10(2),1,2,4)
lab.cutpoints = c(paste0('1e', c(-20,-16,-8,-4)), .01, .1,.5, .9,1, 1.1, 2, 10 ,100, '1e4')

a_plot_m <-
  (td.misc::mcmc_intervals_multi(
    fits = list(Original = a, Missing = a_m),
    pars = c('a0',
             paste0('a[',c(1:5,18,6:7,11:13,15,16,14,17,21,10,19,20,8,9),']')),
    multi_point_est = TRUE,
    point_est = c('mean', 'median'),
    transformations = function(x) ifelse(x==0, 0, sign(x) * sqrt(abs(x))),
    prob_outer = .95) +
     geom_vline(aes(xintercept=0), color=grey(.5), alpha=.5) +
     scale_color_discrete(type=RColorBrewer::brewer.pal(name='Set1', n=3)[c(2,1,3)]) +
     scale_x_continuous(name='TBM odd ratio', 
                        breaks = (function(x) ifelse(x==0, 0, sign(x) * sqrt(abs(x))))(log(10^cutpoints)),
                        labels=lab.cutpoints) +
     ylab('')+
     theme_bw() +
     theme(axis.text.y = ggtext::element_markdown(face='bold'), legend.position = 'bottom')) |>
     # guides(color=guide_legend(title="Scale"))) %>%
     td.misc::change_ylabs(
       labs = labs,
       top_down = T
     ) 

a_plot_m_data_m1 = bayesplot::mcmc_intervals_data(a,
                                                  pars = c('a0',
                                                           paste0('a[',c(1:5,18,6:7,11:13,15,16,14,17,21,10,19,20,8,9),']')),
                                                  transformations = plogis, 
                                                  prob_outer=.95, point_est='mean') %>%
  mutate(across(c(ll:hh), ~ sign(qlogis(.x)) * sqrt(abs(qlogis(.x))))) %>% .$m

a_plot_m_data_m2 = bayesplot::mcmc_intervals_data(a_m,
                                                  pars = c('a0',
                                                           paste0('a[',c(1:5,18,6:7,11:13,15,16,14,17,21,10,19,20,8,9),']')),
                                                  transformations = plogis, 
                                                  prob_outer=.95, point_est='mean') %>%
  mutate(across(c(ll:hh), ~ sign(qlogis(.x)) * sqrt(abs(qlogis(.x))))) %>% .$m

a_plot_m$data$m = c(a_plot_m_data_m1, a_plot_m_data_m2)

b_m = rstan::extract(m3m$outputs, pars=c('b_HIV', 'b'))
b_m$b[,6]=-b_m$b[,6]
b_m = cbind(b_m$b_HIV, b_m$b)
colnames(b_m) = c('b_HIV', paste0('b[',1:(ncol(b)-1), ']'))

b_plot_m = td.misc::mcmc_intervals_multi(list(Original = b, Missing = b_m),
                                       pars = c('b_HIV',  paste0('b[',c(1,6,4,5,2,3), ']')),
                                       point_est = 'mean',
                                       prob_outer = .95) |>
  td.misc::change_ylabs(
    'HIV +',
    'GCS',
    '*log<sub>2</sub>* (CSF Glucose)',
    '*log<sub>2</sub>* (CSF Protein)',
    '*log<sub>2</sub>* (CSF Lactate)',
    '*log<sub>10</sub>* (CSF Lymphocyte)',
    '*log<sub>10</sub>* (CSF WBC)',
    top_down = T
  ) +
  geom_vline(aes(xintercept=0), color=grey(.5), alpha=.5) +
  scale_color_discrete(type=RColorBrewer::brewer.pal(name='Set1', n=3)[c(2,1,3)]) + 
  guides(color=guide_legend(title="Model")) +
  ylab('') +
  scale_x_continuous(name='Standardised bacillary burden', breaks=c(-3,-2,-1,0,1,2,3))+
  theme_bw() + theme(axis.text.y = ggtext::element_markdown(family=NULL, face='bold'), legend.position = 'bottom')

z_m = rstan::extract(m3m$outputs, pars=c('z_Smear', 'z_Mgit', 'z_Xpert', 'z_obs'))
z_m = do.call(cbind, z_m)
colnames(z_m) = paste0('z_',rep(c('Smear', 'Mgit', 'Xpert', 'observe'),each=2), rep(c('[1]', '[2]'), 4))

z_plot_logit_m = (td.misc::mcmc_intervals_multi(list(Original = z, Missing = z_m), 
                              regex_pars = '^z',
                              point_est = 'mean',
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
      "TPR<sub>Smear</sub>",
      "FPR<sub>Smear</sub>"
    ) |> rev(),
    top_down = T
  ) + theme(axis.text.y = ggtext::element_markdown(family=NULL, face='bold'), legend.position = 'bottom')

z_plot_m_data_m1 = bayesplot::mcmc_intervals_data(z,
                                                  regex_pars = '^z',
                                                  transformations = plogis, 
                                                  prob_outer=.95,  
                                                  multi_point_est = TRUE,
                                                  point_est = c('mean', 'median'),) %>%
  mutate(across(c(ll:hh), ~ qlogis(.x))) %>% .$m

z_plot_m_data_m2 = bayesplot::mcmc_intervals_data(z_m,
                                                  regex_pars = '^z',
                                                  transformations = plogis, 
                                                  prob_outer=.95, 
                                                  multi_point_est = TRUE,
                                                  point_est = c('mean', 'median'),)%>%
  mutate(across(c(ll:hh), ~ qlogis(.x))) %>% .$m

z_plot_logit_m$data$m = c(z_plot_m_data_m1, z_plot_m_data_m2)

saveRDS(list(a_plot = a_plot_m, b_plot = b_plot_m, z_plot = z_plot_logit_m), 'export/m3m_plot.RDS')
