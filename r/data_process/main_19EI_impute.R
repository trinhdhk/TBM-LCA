library(data.table)
library(magrittr)

data_19EI <- readRDS('data/cleaned/data_19EI.RDS')
data_19EI <- na.omit(data_19EI, cols = c('csf_smear', 'csf_mgit', 'csf_xpert'))
data_19EI <- data_19EI[, c('hiv_stat', 'age', 'csf_smear', 'csf_mgit', 'csf_xpert',
                           'clin_illness_day', 'clin_gcs', 
                           'clin_nerve_palsy',
                           'csf_clear', 'csf_protein', 'csf_lympho', 'glucose_ratio', 'csf_lactate',
                           'xray_miliary_tb', 'xray_pul_tb', 'ISDIABETE', 'BLDGLU', 'ISNSWEAT', 'ISCOUGH', 'ISWEIGHT',
                           'HEMIPLEGIA', 'PARAPLEGIA', 'TETRAPLEGIA'), with=F]

data_19EI[, `:=`(log_illness_day = log2(clin_illness_day), clin_illness_day=NULL,
                 log_lympho = log2(csf_lympho+1),csf_lympho=NULL,
                 log_protein = log2(csf_protein),csf_protein=NULL,
                 log_lactate = log2(csf_lactate),csf_lactate=NULL)]
set.seed(2906)

pred.mat <- mice::make.predictorMatrix(data_19EI)
pred.mat[3:5,] <- 0
pred.mat[14:19,] <- 0
pred.mat[14:16,c(14:16, 21:23)] <- 1
pred.mat[14:16,1] <- 1
pred.mat[17:19,] <- 0
pred.mat[17:19,c(1,12:13,17:19, 21:23)] <- 1
pred.mat[17:19,1] <- 1
pred.mat[, 3:5] <- 1

data_19EI_imp <- mice::parlmice(data_19EI, n.core = 18, maxit=100, m=126, n.imp.core = 7, predictorMatrix = pred.mat, defaultMethod = c("pmm", "pmm", "polyreg", "polr"))
data_19EI_complete <- lapply(seq_len(data_19EI_imp$m), mice::complete, data = data_19EI_imp)
data_19EI_complete <- lapply(data_19EI_complete, function(dt)
  dplyr::mutate(dt,clin_symptoms = ISNSWEAT|ISCOUGH|ISWEIGHT, clin_motor_palsy = HEMIPLEGIA|PARAPLEGIA|TETRAPLEGIA))
saveRDS(data_19EI_complete, 'data/impute_19EI.RDS')
