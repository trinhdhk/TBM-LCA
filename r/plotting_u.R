z_plot_u =td.misc::mcmc_intervals_multi(
  fits = list(Selected = m3$outputs,
              "Incomplete Xpert" = m3mx$outputs,
              "Incomplete Xpert + Wide prior" = m3mxu$outputs),
  regex_pars = '^z_(Smear|Xpert|Mgit|obs)\\[',
  multi_point_est = TRUE,
  point_est = c('mean', 'median'),
  point_size=2,
  prob_outer = .95,
  transformations = {
    a <- list()
    within(a, {
      `z_Xpert[1]` <- `z_Mgit[1]` <- `z_Smear[1]` <- function(i) -i
      `z_Xpert[2]` <- `z_Mgit[2]` <- `z_Smear[2]` <- I
    
    })
  }
) |>
  td.misc::change_ylabs(
    'Spc<sub>ZN-Smear</sub>',
    'Sen<sub>ZN-Smear</sub>',
    'Spc<sub>MGIT</sub>',
    'Sen<sub>MGIT</sub>',
    'Spc<sub>Xpert</sub>',
    'Sen<sub>Xpert</sub>',
    'Xpert observed / non-TBM',
    'Xpert observed / TBM'
  ) + 
  ggplot2::theme_bw() + 
  ggplot2::scale_x_continuous(
    breaks = qlogis(c(.02, .1, .25, .5, .75, .9, .99, .999)),
    labels = c(.02, .1, .25, .5, .75, .9, .99, .999)
  ) + 
  ggplot2::theme(axis.text.y = ggtext::element_markdown(face='bold'), 
                 axis.title = ggplot2::element_blank(),
                 legend.position = 'bottom')

z_plot_dat = td.misc::mcmc_intervals_multi(
  fits = list(Selected = m3$outputs,
              "Incomplete Xpert" = m3mx$outputs,
              "Incomplete Xpert + Wide prior" = m3mxu$outputs),
  regex_pars = '^z_(Smear|Xpert|Mgit|obs)\\[',
  multi_point_est = TRUE,
  point_est = c('mean', 'median'),
  prob_outer = .95,
  transformations = {
    a <- list()
    within(a, {
      `z_Xpert[1]` <- `z_Mgit[1]` <- `z_Smear[1]` <- function(i) plogis(-i)
      `z_Xpert[2]` <- `z_Mgit[2]` <- `z_Smear[2]` <- plogis
      
    })
  }
)$data

z_plot_u$data$m[! z_plot_u$data$parameter %in% c("z_obs[1]", "z_obs[2]")] = qlogis(z_plot_dat$m[! z_plot_dat$parameter %in% c("z_obs[1]", "z_obs[2]")])
saveRDS(ggplotGrob(z_plot_u), 'export/z_plot_u.RDS')
