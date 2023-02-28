ppc = readRDS('outputs/ppc_check.RDS')
orig = readRDS('outputs/m3_t00_b345678_q7_r1_k1_3012.RDS')
ppc_merged = rstan::sflist2stanfit(ppc)

# We will visually check for a[1], a[2], a[3], a[11], a[18]
extract_param_data <- function(par){
  ppc = ppc
  if (!length(names(ppc))) names(ppc) = seq_along(ppc)
  param = list()
  param$orig = rstan::extract(orig$outputs, par)[[par]]
  param$rep = sapply(ppc, \(x) rstan::extract(x, par)[[par]], simplify = FALSE)
  param_data = data.frame(val = do.call(c, param$rep), 
                          id =  rep(names(param$rep), each=length(param$rep[[1]]))) |>
    rbind(data.frame(val = param$orig, id = 0))
  class(param_data) = c('param_data', 'data.frame')
  param_data
}

plot.param_data = function(x, grob=TRUE, ..., plot=TRUE){
  require(ggplot2)
  x_mean = dplyr::with_groups(dplyr::tibble(x), id, dplyr::summarise, med= median(val))
  plt = ggplot(data=x) + 
    geom_density(aes(x=val, color=id), adjust=1, size=.7) + 
    gghighlight::gghighlight(id == 0, use_group_by=FALSE, unhighlighted_params = list(size=.5, color = grey(.8))) +
    geom_segment(aes(y = -.05, yend = -.05, x = quantile(med, .025), xend = quantile(med, .975)), data=x_mean[x_mean$id!=0,], color='black')+
    geom_point(aes(y = -0.05, x=med), data=x_mean[x_mean$id==0,], color='red', size=1) +
    # ggdist::geom_pointinterval(aes(x = med, xmin=.lower, xmax=.upper), data=x_mean) +
    theme_classic() + 
    scale_color_brewer(type='qual', palette='Set1')+
    xlab('Log odds ratio') + 
    ylab('Density')+
    theme(legend.position='none')
  if (!grob) return(plt)
  if (plot) plot(plt)
  ggplotGrob(plt)
}

plt = list()

a1 = extract_param_data('a[1]')
plt$a1 = plot(a1)

a2 = extract_param_data('a[2]')
plt$a2 = plot(a2)

a3 = extract_param_data('a[3]')
plt$a3 = plot(a3)

a11 = extract_param_data('a[11]')
plt$a11 = plot(a11)

a18 = extract_param_data('a[18]')
plt$a18 = plot(a18)

a0 = extract_param_data('a0')
plt$a0 = plot(a0)

saveRDS(plt, 'export/ppc_checkresults.RDS')
