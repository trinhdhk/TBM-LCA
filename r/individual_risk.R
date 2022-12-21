library(data.table)
library(ggplot2)
model <- readRDS('outputs/m3_t00_b345678_q7_r1_k1.RDS')$outputs
theta <- rstan::extract(model, 'theta')$theta
jsonlite::write_json(theta, 'export/m3_theta.json', )
mean_theta <- apply(theta, 2, mean)
plt <- ggplot(mapping=aes(x=mean_theta)) + geom_histogram(fill='#B5B991FF') 
plt_data <- ggplot_build(plt)$data[[1]]
bin_theta <- apply(theta, 1, 
    \(x) cut(x, breaks = c(0, plt_data$xmax)) |> table() |> c())

upper_ci <- \(x) quantile(x, .975)
lower_ci <- \(x) quantile(x, .025)
plt + geom_errorbar(aes(x = plt_data$x,
    ymin = apply(bin_theta,1,lower_ci),
    ymax = apply(bin_theta,1,upper_ci)),
    width = .01) +
    theme_bw() + 
    labs(x='TBM probability', y='Count')
    # scale_y_continuous(trans = 'log2')
