# Load data
# Author: trinhdhk
# Version 0.1.2005

# 1. Load data
cleaned_data_dir <- 'data/cleaned'
data.files <- list.files(cleaned_data_dir, '.RDS')
# data.names <- gsub('.RDS$', '', data.files)
# for (f in data.files)
#   assign(gsub('.RDS$', '', f), readRDS(file.path(cleaned_data_dir, f)))
data_23TB <- readRDS(file.path(cleaned_data_dir, 'data_23TB.RDS'))
data_05TB <- readRDS(file.path(cleaned_data_dir, 'data_05TB.RDS'))
# data.names <- gsub('.RDS$', '', data.files)
# cols <- Reduce(intersect, lapply(data.names, function(x) names(get(x))))

# data.names <- data.names[data.names!='data_23TB']
# maindt <- dplyr::bind_rows(
  # lapply(data.names, function(x){
  #   x <- get(x, envir=globalenv())
  #   x[, cols, with=FALSE]
  # }))
# maindt <- data_05TB
maindt <- rbind(data_23TB,
                data_05TB[,names(data_23TB),with = FALSE])[,-1]
maindt[, lym_glu_ratio := ((csf_lympho+1)/glucose_ratio)]
maindt<-maindt[hiv_stat==0, `:=`(
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
X <- maindtcomplete %$% as.matrix(cbind(#hiv_stat, 
                                        log2(age), 
                                        sex, 
                                        log2(bmi), 
                                        clin_illness_day, #remove contact_tb due to many "unknown"s
                                        clin_nerve_palsy,
                                        clin_motor_palsy,
                                        sqrt(csf_lympho),
                                        sqrt(glucose_ratio),
                                        sqrt(csf_protein)
                                        # xray_pul_tb
                                        ))

X_full <- maindt %$% as.matrix(cbind(#hiv_stat, 
                                     age, 
                                     sex, 
                                     bmi, 
                                     clin_illness_day, #remove contact_tb due to many "unknown"s
                                     clin_nerve_palsy,
                                     clin_motor_palsy,
                                     csf_lympho,
                                     glucose_ratio,
                                     csf_protein,
                                     xray_pul_tb))

model_input_cont <- maindtcomplete %$% 
  list(
    N = nrow(maindtcomplete), nX = ncol(X),
    Y_Smear = as.integer(csf_smear),
    Y_Mgit = as.integer(csf_mgit),
    Y_Xpert = as.integer(csf_xpert),
    Y_Img = as.integer(xray_miliary_tb),
    # Y_lympho = csf_lympho,
    Y_lnLymGlu = lym_glu_ratio,
    X = X
  )

model_input_disc <- maindtcomplete %$%
  list(
    N = nrow(maindtcomplete), nX = ncol(X),
    Y_Smear = as.integer(csf_smear),
    Y_Mgit = as.integer(csf_mgit),
    Y_Xpert = as.integer(csf_xpert),
    Y_Img = as.integer(xray_miliary_tb),
    # Y_CSF = (csf_wbc>10&csf_wbc<500) & (csf_lym_pct>.5) & (glucose_ratio<.5|csf_glucose<2.2) & (csf_protein>1),
    X = X
  )

model_input_fiml <- maindt %$%
  list(
    N = nrow(maindt), nX = ncol(X), 
    # N_Smear = length(na.omit(csf_smear)), N_Mgit = length(na.omit(csf_mgit)),
    # N_Xpert = length(na.omit(csf_xpert)), N_Img = length(na.omit(img_score)),
    # N_Xray = length(na.omit(xray_miliary_tb)), N_CSF = length(na.omit(xray_miliary_tb)),
    Y_Smear = as.integer(tidyr::replace_na(csf_smear, 0)), obs_Smear = (!is.na(csf_smear)), 
    Y_Mgit = as.integer(tidyr::replace_na(csf_mgit, 0)), obs_Mgit = (!is.na(csf_mgit)),
    Y_Xpert = as.integer(tidyr::replace_na(csf_xpert, 0)), obs_Xpert = (!is.na(csf_xpert)),
    Y_Img = as.integer(tidyr::replace_na(img_score, 0)), obs_Img = (!is.na(img_score)),
    Y_Xray = as.integer(tidyr::replace_na(xray_miliary_tb, 0)), obs_Xray = (!is.na(xray_miliary_tb)),
    Y_CSF = tidyr::replace_na(csf_lym_pct>.5 & (glucose_ratio<.5|csf_glucose<2.2) & csf_protein>1, 0), obs_CSF = (!is.na(csf_lym_pct+glucose_ratio+csf_glucose+csf_protein)),
    X = tidyr::replace_na(X_full, 0), obs_X = !apply(X_full, c(1,2), is.na)
  )
