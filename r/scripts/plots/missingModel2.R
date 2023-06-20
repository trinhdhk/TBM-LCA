#Script to generate plots of coefficient for selected and models with missing Xpert
# Author: Trinh Dong
# Email: trinhdhk@oucru.org
z = rstan::extract(m3$outputs, pars=c('z_Smear', 'z_Mgit', 'z_Xpert')) |> 
  do.call(what=cbind,) |>
  `colnames<-`(paste0(rep(c('z_Smear', 'z_Mgit', 'z_Xpert'), each=2), "[", 1:2, "]"))
z_mX = rstan::extract(m3m$outputs, pars=c('z_Smear', 'z_Mgit', 'z_Xpert', 'z_obs')) |>
  do.call(what=cbind,) |>
  `colnames<-`(paste0(rep(c('z_Smear', 'z_Mgit', 'z_Xpert', 'z_obs'), each=2), "[", 1:2, "]"))
z_mXu= rstan::extract(m3mu$outputs, pars=c('z_Smear', 'z_Mgit', 'z_Xpert', 'z_obs')) |>
  do.call(what=cbind,) |>
  `colnames<-`(paste0(rep(c('z_Smear', 'z_Mgit', 'z_Xpert', 'z_obs'), each=2), "[", 1:2, "]"))
z_plot_u =td.misc::mcmc_intervals_multi(
  fits = list(Selected = z,
              "Incomplete Xpert" = z_mX,
              "Incomplete Xpert + Wide prior" = z_mXu),
  regex_pars = '^z_(Smear|Xpert|Mgit|obs)\\[',
  multi_point_est = FALSE,
  point_est =  c('median'),
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
  scale_color_discrete(type=RColorBrewer::brewer.pal(name='Set1', n=3)[c(2,1,3)]) +
  ggplot2::theme(axis.text.y = ggtext::element_markdown(), 
                 axis.title = ggplot2::element_blank(),
                 legend.position = 'bottom')

saveRDS((z_plot_u), 'export/z_plot_u.RDS')
