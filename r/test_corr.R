
library(data.table)
library(ggplot2)
load('data/cleaned/data_input.Rdata')
setDT(data_19EI)
data_19EI[, 
          culture_time := sapply(
            TimeToPositive, 
            function(time){
              if (is.na(time)) return(NA)
              time_str = strsplit(time, ';|:|,')[[1]]
              if (length(time_str)==1) return(as.integer(time_str) * 24)
              as.integer(time_str[[1]])*24 + as.integer(time_str[[2]])
            })]

# model = readRDS('outputs/m3_t00_b345678_q7_r1_k1_3012.RDS')$outputs
model = m3$outputs
X = rstan::extract(model, 'X')$X
b = rstan::extract(model, c('b', 'b_HIV', 'b_RE'))
RE = rstan::extract(model, 'RE')$RE
burden = sapply(
  seq_len(dim(X)[1]),
  \(i){
    hiv = X[i,,1]
    bd = X[i,,13:18]
    re = RE[i,]
    b_HIV = b$b_HIV[i]
    b_bd = b$b[i,,drop=FALSE]
    (hiv * b_HIV + bd %*% t(b_bd) + re)
  }) 

ztheta = rstan::extract(model, 'z_theta')$z_theta

mean_burden = apply(burden,1,mean)
mean_ztheta = apply(ztheta,2,mean)

# bind together
bd_dt = 
  data.table(
    burden = mean_burden,
    xpert_level = data_19EI[,ifelse(is.na(XpertLevel), 'Undetected', XpertLevel)] |> factor(levels = c('Undetected', 'Trace','Very low','Low','Medium'),exclude = NULL, ordered=T),
    culture_time = data_19EI[,ifelse(is.na(culture_time),56,culture_time)],
    growth_unit = data_19EI$GrowthUnit,
    xpert = data_19EI$csf_xpert,
    mgit = data_19EI$csf_mgit,
    tbm_dx = data_19EI[,tbm_dx|csf_smear|csf_mgit|csf_xpert],
    vol = data_19EI$Volume,
    hiv = data_19EI$hiv_stat
  )
bd_dt[, log_culture_time := (log2(culture_time+1))][,  gph := log(growth_unit*culture_time)]
# check correpondence with Xpert Level
xpert_level = 
  ggstatsplot::ggbetweenstats(bd_dt[xpert==T], x=xpert_level, y=burden, plot.type='box', pairwise.comparison=F, type='non-parametric', pairwise.display='none',
                              results.subtitle = FALSE,
                              xlab = "Xpert semi-quantification level", ylab = "Standardised Mycobacillary burden")
# check correspondence with culture time
cultime = list(
  pos=ggstatsplot::ggscatterstats(bd_dt[mgit==TRUE & growth_unit>0 & hiv], x=culture_time, y=burden, type = 'non-parametric', smooth.line.args = list(method='loess',span=2, se=F),
                                  results.subtitle = FALSE,
                                  xlab = "Time to Culture positive (hours)", ylab = "Standardised Mycobacillary burden") + ggtitle('HIV +'),
  neg=ggstatsplot::ggscatterstats(bd_dt[mgit==TRUE & growth_unit>0 & !hiv], x=culture_time, y=burden, type = 'non-parametric', smooth.line.args = list(method='loess',span=2, se=F),
                                  results.subtitle = FALSE,
                                  xlab = "Time to Culture positive (hours)", ylab = "Standardised Mycobacillary burden") + ggtitle('HIV -'),
  overall=ggstatsplot::ggscatterstats(bd_dt[mgit==TRUE & growth_unit>0], x=culture_time, y=burden, type = 'non-parametric', smooth.line.args = list(method='loess',span=2, se=F),
                                      results.subtitle = FALSE,
                                      xlab = "Time to Culture positive (hours)", ylab = "Standardised Mycobacillary burden")
)
cultimefit = lm(burden~hiv*log_culture_time,data=bd_dt[mgit==T])

# with current score
score_dt = data.table(
  score = data_19EI$crude_total_score,
  score2 = dplyr::case_when(
    # data_19EI$crude_total_score == 0 ~ 7,
    data_19EI$crude_total_score < 6 ~ data_19EI$crude_total_score + 5,
    data_19EI$crude_total_score == 6 ~ 11,
    T ~ data_19EI$crude_total_score
  ),
  ztheta = mean_ztheta,
  tbm_dx = fcase(data_19EI[,csf_smear|csf_mgit|csf_xpert], "Confirmed TBM", 
                 data_19EI$tbm_dx, "Clinical TBM",
                 default = "Not TBM") |> tidyr::replace_na('Not TBM')
)

# score_dt[, logscore2 := qlogis((score2-6)/8)]
area = expand.grid(score = seq(0, 14, .01), y=seq(-25, 20,.1)) |> setDT()
# score_dt[, def := fcase(score < 6, "Unknown", score < 9, "Possible TBM", default = "Probable TBM")]
area[, def := fcase(score < 6, "Unknown", score < 9, "Possible TBM", default = "Probable TBM") |> factor(levels = c('Unknown', 'Possible TBM', 'Probable TBM'), ordered=TRUE)]
defscore = ggstatsplot::ggscatterstats(score_dt, x=score, y=ztheta, type = 'non-parametric', smooth.line.args = list(method='loess',span=.8,se=F),
                                       point.args = list(alpha=0),
                                      xlab = "Total score", ylab = "Predicted TBM risk (%)",
                                      results.subtitle = FALSE,
                                      ggplot.component = scale_x_continuous(breaks = c(0, 3, 6, 9, 12, 14))) +
  geom_jitter(aes(color = tbm_dx), size=1, na.rm=TRUE, width=.3, height=1)+
  scale_color_brewer(type='qual', palette=7, name='Dx at discharge', guide = guide_legend(nrow=2, title.position = 'top')) +
  geom_raster(aes(x=score, y=y, fill = def), data=area, alpha=.3) + 
  scale_fill_brewer(palette=7, name='Case definition', guide = guide_legend(nrow=2, title.position = 'top')) +
  scale_y_continuous(breaks = qlogis(c(1e-10, 1e-4, 0.01, .1, .5, .9, .99, .9999, .99999999)), labels = c(1e-8, 1e-2, 1, 10, 50, 90, 99, 99.99, 100)) +
  theme(legend.position='bottom')

defscore2 = ggstatsplot::ggscatterstats(score_dt, x=score2, y=ztheta, type = 'non-parametric', smooth.line.args = list(method='loess',span=.8,se=F),
                                       point.args = list(alpha=0),
                                       xlab = "Total score", ylab = "Predicted TBM risk (%)",
                                       results.subtitle = FALSE,
                                       ggplot.component = scale_x_continuous(breaks = c(6, 9, 12, 14))) +
  geom_jitter(aes(color = tbm_dx), size=1, na.rm=TRUE, width=.3, height=1)+
  scale_color_brewer(type='qual', palette=7, name='Dx at discharge', guide = guide_legend(nrow=2, title.position = 'top')) +
  geom_raster(aes(x=score, y=y, fill = def), data=area[score>=6], alpha=.3) + 
  scale_fill_brewer(palette=7, name='Case definition', guide = guide_legend(nrow=2, title.position = 'top')) +
  scale_y_continuous(breaks = qlogis(c(1e-10, 1e-4, 0.01, .1, .5, .9, .99, .9999, .99999999)), labels = c(1e-8, 1e-2, 1, 10, 50, 90, 99, 99.99, 100)) +
  theme(legend.position='bottom')

saveRDS(list(
  defscore = defscore |> ggplotGrob(),
  defscore2 = defscore2 |> ggplotGrob(),
  defscorecombined = defscore+(defscore2+scale_fill_brewer(palette=7, name='Case definition', guide = guide_none())+theme(axis.title.y = element_blank()))+patchwork::plot_layout(guide='collect')&theme(legend.position = 'bottom'),
  xpertlv = xpert_level |> ggplotGrob(),
  cultime = cultime |> sapply(ggplotGrob, simplify = F),
  cultimefit = cultimefit
), file = 'export/test_corr.RDS')
