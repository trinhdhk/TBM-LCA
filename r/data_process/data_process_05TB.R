# Data processing for 05TB data
# Author: trinhdhk
# Ver 0.1.2004.2

library(data.table)
rm(list=ls())
diag  <- haven::read_sas("data/raw/05TB/SASvad/tbmdiagnosis.sas7bdat")
assay <- haven::read_sas("data/raw/05TB/SASvad/bl_inv.sas7bdat")
hist <- haven::read_sas("data/raw/05TB/SASvad/bl_demohist.sas7bdat")
neuro <- haven::read_sas("data/raw/05TB/SASvad/bl_neurotbm.sas7bdat")
# bl <- readxl::read_excel(file.path('data', 'raw', clinical_data_file), sheet='SCREX')

setDT(diag)
setDT(assay)
setDT(hist)
setDT(neuro)
# setDT(bl)

hist[, `:=`(clin_symptoms=(weightloss=='yes')|(nightsweats=='yes'),
            clin_contact_tb=contact_tb1y=='yes',
            hiv_stat = hiv_strat == 'positive')]
neuro[, `:=`(clin_nerve_palsy = cnp=='yes',
            clin_motor_palsy = hemiplegia == 'yes' | paraplegia == 'yes' | quadriplegia == 'yes'), by=ParNo]

assay[, `:=`(
  xray_miliary_tb=(xrayresult=='abnormal miliary TB'), 
  xray_pul_tb=(xrayresult=='abnormal consistent with TB'),
  csf_lympho=CsfBC_bl*CsfLympho_bl/100,
  csf_clear=csf_appearance_screen=='clear'
), by=ParNo]

diag[, `:=`(
  img_hydro=cere_1,
  img_basal=cere_2,
  img_infarct=cere_4,
  img_tuber=cere_3,
  img_precon=cere_5,
  img_score=cere_score
)]

maindt <- plyr::join_all(
  list(
    diag[,.(ParNo, 
            img_hydro, img_basal, img_infarct, img_tuber, img_precon, img_score,
            csf_smear = csf_afb_pos==1, csf_mgit = csf_mtb_pos==1, csf_xpert = csf_Xpert_pos==1,
            cere_score, crude_total_score, diagnostic_score)],
    hist[,.(ParNo, age, sex=as.numeric(sex=='male'), bmi, hiv_stat, clin_illness_day=illness_duration, clin_symptoms, clin_contact_tb)],
    neuro[,.(ParNo, clin_gcs=gcs_bl, clin_nerve_palsy, clin_motor_palsy)],
    assay[,.(ParNo,csf_clear,csf_lympho,csf_lym_pct=CsfLympho_bl,csf_wbc=CsfBC_bl,csf_protein=CsfProtein_bl,
             csf_lactate=CsfLactate_bl,glucose_ratio=CsfGlucose_bl/BlGlucose_bl,csf_glucose=CsfGlucose_bl,
             xray_miliary_tb,xray_pul_tb)]
  ),
  by = 'ParNo',
  type = 'left'
)[diagnostic_score!='confirmed other diagnosis'][, diag := factor(
  dplyr::case_when(
    is.na(cere_score) & crude_total_score >= 10 | crude_total_score >= 12 ~ 'Probable',
    !is.na(crude_total_score) ~ 'Possible'
  ),
  levels = c('Possible', 'Probable'),
  ordered = TRUE
)]

setnames(maindt, 
         c('cere_score', 'crude_total_score', 'diagnostic_score'), 
         c('csf_score', 'total_score', 'diag'))

saveRDS(maindt, file='data/cleaned/data_05TB.RDS')
