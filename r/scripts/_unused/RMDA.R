rmda_dt = data.table(
  Dx = m3$recipe$data_19EI[, as.numeric(tbm_dx | csf_smear | csf_mgit | csf_xpert)],
  prob = m3$p$theta$mean,
  smear = m3$recipe$data_19EI[, csf_smear],
  mgit = m3$recipe$data_19EI[, csf_mgit],
  xpert = m3$recipe$data_19EI[, csf_xpert],
  possible = m3$recipe$data_19EI[, crude_total_score>=6],
  probable = m3$recipe$data_19EI[, crude_total_score>=9]
)
rmda_dt[, test := smear | mgit | xpert]
des.curve=list()
des.curve$model = rmda::decision_curve(Dx ~ prob, data = rmda_dt, fitted.risk = TRUE, policy='opt-in')
des.curve$test = rmda::decision_curve(Dx ~ test, data = rmda_dt, policy='opt-in', fitted.risk = T)
# des.curve$mgit = rmda::decision_curve(Dx ~ mgit, data = rmda_dt, policy='opt-in', fitted.risk = T)
# des.curve$xpert = rmda::decision_curve(Dx ~ xpert, data = rmda_dt, policy='opt-in', fitted.risk = T)
des.curve$possible = rmda::decision_curve(Dx ~ possible, data = rmda_dt, policy='opt-in', fitted.risk = T)
des.curve$probable = rmda::decision_curve(Dx ~ probable, data = rmda_dt, policy='opt-in', fitted.risk = T)
rmda::plot_decision_curve(des.curve, 
                          curve.names = c('Model', 'Test', 'Possible TBM', 'Probable TBM'),
                          confidence.intervals = 'none', standardize=FALSE, legend.position = 'bottomright', xlim=c(0, 1.33))
saveRDS(des.curve, 'export/des_curve.RDS')
