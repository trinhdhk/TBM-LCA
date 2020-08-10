# Data processing for 23TB data
# Author: trinhdhk
# Ver 0.1.2004

library(data.table)
diag  <- haven::read_sas("data/raw/SASvad/tbmdiagnosis.sas7bdat")
assay <- haven::read_sas("data/raw/SASvad/bl_inv.sas7bdat")
hist <- haven::read_sas("data/raw/SASvad/bl_demohist.sas7bdat")
neuro <- haven::read_sas("data/raw/SASvad/bl_neurotbm.sas7bdat")

data_files <- list.files('data/raw', pattern = '^\\w')
clinical_data_file <-  grep('DATA final', data_files, value = TRUE)
scr <- readRDS('data/raw/raw_23TB.RDS')$tables$SCR
bl <- readRDS('data/raw/raw_23TB.RDS')$tables$BAS

setDT(diag)
setDT(assay)
setDT(hist)
setDT(neuro)
setDT(scr)
setDT(bl)

hist[, `:=`(clin_symptoms=(WeightLoss=='Yes')%in%T|(NightSweats=='Yes')%in%T,
            clin_contact_tb=IsContact=='Yes')]

neuro[, clin_nerve_palsy := CranialPalsy=='Yes']

bl[, clin_motor_palsy := (Paraplegia=='Yes')%in%T|(Hemiplegia=='Yes')%in%T|(Tetraplegia=='Yes')%in%T]

assay[, `:=`(
  xray_miliary_tb=(MiliaryTB==1)%in%T, 
  xray_pul_tb=(TypicalPulTB==1)%in%T|(SuspectedPulTB==1)%in%T,
  csf_lympho=whitecells_bl*lym_bl/100, csf_protein=protein_bl,
  csf_clear=Appearance_bl=='Clear',
  glucose_ratio=csfglu_ratio_bl,
  img_hydro=hydrocephalus_bl,
  img_basal=Meningeal_Enhancement_bl,
  img_infarct=infarcts_bl,
  img_tuber=tuberculomas_bl,
  img_precon=communicating_bl, #???
  img_score=sum(hydrocephalus_bl, Meningeal_Enhancement_bl, infarcts_bl, tuberculomas_bl, communicating_bl,
                na.rm=TRUE)
), by=PatientCode]

maindt <- plyr::join_all(
  list(
    diag[,.(PatientCode, 
            csf_smear = csf_afb_pos==1, csf_mgit = csf_mtb_pos==1, csf_xpert = csf_Xpert_pos==1,
            cere_score, crude_total_score, diagnostic_score)],
    hist[,.(PatientCode,age, sex=as.numeric(Gender=='Male'), bmi, clin_illness_day=illness_duration, clin_symptoms, clin_contact_tb)],
    neuro[,.(PatientCode, clin_gcs=gcs_bl, clin_nerve_palsy)],
    bl[,.(PatientCode, clin_motor_palsy)],
    assay[,.(PatientCode,csf_clear,csf_lympho,csf_lym_pct=lym_bl,csf_wbc=whitecells_bl,csf_protein=protein_bl,
          csf_lactate=lactate_bl,glucose_ratio,csf_glucose=csfglu_bl,
          xray_miliary_tb,xray_pul_tb, img_hydro, img_basal, img_infarct, img_tuber, img_precon, img_score)],
    scr[,.(PatientCode, hiv_stat=as.numeric(HIV == 'Yes'))]
  ),
  by = 'PatientCode',
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
         c('PatientCode', 'cere_score', 'crude_total_score', 'diagnostic_score'), 
         c('ParNo', 'csf_score', 'total_score', 'diag'))

saveRDS(maindt, file='data/cleaned/data_23TB.RDS')
