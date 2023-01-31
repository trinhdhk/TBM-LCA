library(data.table)
library(ggplot2)
model <- readRDS('outputs/m3_t00_b345678_q7_r1_k1_2912.RDS')$outputs
theta <- rstan::extract(model, 'theta')$theta
jsonlite::write_json(theta, 'export/m3_theta.json', )
mean_theta <- apply(theta, 2, median)
plt <- ggplot(mapping=aes(x=mean_theta, fill = data_19EI[, csf_smear+csf_xpert+csf_mgit>0])) + 
  geom_histogram() + ggsci::scale_fill_jco()
   #ill='#B5B991FF' 
plt_data <- ggplot_build(plt)$data[[1]]
bin_theta <- apply(theta, 1, 
    \(x) cut(x, breaks = c(0, unique(plt_data$xmax))) |> table() |> c())
bin_theta_confirmed <- apply(theta[,data_19EI[,(csf_smear + csf_mgit + csf_xpert > 0)]], 1, 
    \(x) cut(x, breaks = c(0, unique(plt_data$xmax))) |> table() |> c())
plt_data$count = c(
  apply(bin_theta-bin_theta_confirmed, 1,median), apply(bin_theta_confirmed,1,median)
)

upper_ci <- \(x) quantile(x, .975)
lower_ci <- \(x) quantile(x, .025)
plt2 = ggplot() + 
  stat_identity(aes(x=x, y=count, fill=rep(c( "#C7C7C7", '#FC8D62'), each=length(x)/2)), data=plt_data, geom='bar')+
  scale_fill_identity()+
  geom_errorbar(aes(x = unique(plt_data$x),
    ymin = apply(bin_theta,1,lower_ci),
    ymax = apply(bin_theta,1,upper_ci)),
    width = .01, inherit.aes = F) +
    theme_bw() + 
    labs(x='TBM probability', y='Count')

