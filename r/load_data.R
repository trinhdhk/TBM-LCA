# Load data
# Author: trinhdhk
# Version 0.1.2005

# 1. Load data
rm(list=ls())
cleaned_data_dir <- 'data/cleaned'
data.files <- list.files(cleaned_data_dir, '.RDS')
# data.names <- gsub('.RDS$', '', data.files)
for (f in data.files)
  assign(gsub('.RDS$', '', f), readRDS(file.path(cleaned_data_dir, f)))
# data_23TB <- readRDS(file.path(cleaned_data_dir, 'data_23TB.RDS'))
# data_05TB <- readRDS(file.path(cleaned_data_dir, 'data_05TB.RDS'))
data.names <- gsub('.RDS$', '', data.files)
cols <- Reduce(intersect, lapply(data.names, function(x) names(get(x))))

maindt <- dplyr::bind_rows(
  lapply(data.names, function(x){
    x <- get(x, envir=globalenv())
    x[, cols, with=FALSE]
  }))
# maindt <- rbind(data_23TB,
#                 data_05TB[,names(data_23TB),with = FALSE])[,-1]
maindt[, lym_glu_ratio := ((csf_lym_pct+1)/glucose_ratio)/100]
maindt[, `:=`(
  clin_symptoms = as.numeric(clin_symptoms),
  clin_nerve_palsy = as.numeric(clin_nerve_palsy),
  clin_motor_palsy = as.numeric(clin_motor_palsy),
  csf_clear = as.numeric(csf_clear),
  xray_miliary_tb = as.numeric(xray_miliary_tb),
  xray_pul_tb = as.numeric(xray_pul_tb)
)]

# 2. Build input

maindtcomplete <- na.omit(maindt, cols = c('age', 'sex', 'bmi', 'csf_smear', 'csf_mgit', 'csf_xpert', 'lym_glu_ratio', #'img_score'
                                           'hiv_stat','clin_illness_day', 'clin_symptoms', 'clin_gcs', 
                                           'clin_nerve_palsy', 'clin_motor_palsy', 'csf_clear', 'csf_protein',
                                           'xray_miliary_tb', 'xray_pul_tb'))
X <- maindtcomplete %$% as.matrix(cbind(hiv_stat, age, sex, bmi, clin_illness_day, #remove contact_tb due to many "unknown"s
                                        clin_nerve_palsy, clin_motor_palsy, glucose_ratio, csf_protein,
                                        xray_pul_tb))

model_input <- maindtcomplete %$% 
  list(
    N = nrow(maindtcomplete), nX = ncol(X),
    Y_Smear = csf_smear,
    Y_Mgit = csf_mgit,
    Y_Xpert = csf_xpert,
    Y_Img = xray_miliary_tb,
    Y_lympho = csf_lympho,
    # Y_lnLymGlu = log(lym_glu_ratio),
    X = X
  )
