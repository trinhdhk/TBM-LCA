# Data process code for 19EI
# Author: Trinhdhk
# Version: 2101.1

library(data.table); library(dplyr)
data.folder <- 'data/raw/19EI'
data.file <- file.path(data.folder, "1-3-2021-_19EI_V1_Data.xls")
sheets <- readxl::excel_sheets(data.file)
data <- sapply(sheets, function(sheet) setDT(readxl::read_excel(data.file, sheet=sheet)))
myco <- readRDS("data/cleaned//mycol19.RDS")

joindt <- plyr::join_all(
  list(
    data$INEX[, .(USUBJID, INITIAL, YOB, SEX, ENROLL)],
    data$BASE[, .(USUBJID, ISHIV, ISCHRONIC, ISDIABETE, ISFEVER, ISCONTUBER, ILLNESSDAY, 
             ISREDUCED, ISWEIGHT, ISNSWEAT, ISCOUGH, NONE, ISHEADACHE, ISLOCALSEIZURE, ISGENNERALSEIZURE, 
             ISPSYCHOSIS, ISLANGCHANGE, ISMOVEMENT,
             HEMIPLEGIA, PARAPLEGIA, TETRAPLEGIA, 
             GENCONVUL, LOCALCONVUL, 
             GLASCOW, GCSE, GCSM, GCSV)],
    data$BASE_LP[, .(USUBJID, LUMBARDATE, APP_A, GRAM, INDIAINK, CRYTO,
                     BACSPE, OTH, JEVIGM, DENGUEIMG, DENGUEPCR, JEVCSF, VIRO_OTH, NMDAR, #other diag
               CSF_VOL= fcase(#!is.na(as.numeric(TB_VOL)) & (TB_VOL != "" %in% T), as.numeric(TB_VOL),
                              !is.na(as.numeric(CSFVOLUME)) & (CSFVOLUME != "" %in% T), as.numeric(CSFVOLUME),
                              !is.na(as.numeric(CSFVOL)) & (CSFVOL != "" %in% T), as.numeric(CSFVOL),
                              default=NA),
               REDCELL, WHITECELL, NEUPER, LYMPER, EOSPER,
               PROTEIN, CSFGLU, BLDGLU, CSFLAC,
               #ZN,
               TBNAAT, NAATSPEC, 
               MYCORESULT, GRAM, BACCUL)],
    data$IMAG[, .(USUBJID, PTB, MTB, CTRESULT)],
    data$ADDINFO[, .(USUBJID, HIV, DENIGM, DENNS1, ADD_JEVIGM = JEVIGM, HEMO, WHITE, NEUTRO, LYMP, EOSI, PLATE)],
    data$OUTC[, .(USUBJID, DISDIA, ISTBTREAT, TBTREATDATE=TBDATE)],
    
    myco
  ), 
  by='USUBJID',
  type='left'
)

joindt <- mutate(joindt,
                 across(c(ISCHRONIC, ISFEVER,
                          ISWEIGHT, ISNSWEAT, ISCOUGH, ISCONTUBER, HEMIPLEGIA, PARAPLEGIA, TETRAPLEGIA, ISREDUCED),
                        function(x) {fcase(x == 'C49488', TRUE,
                                           x == 'N', FALSE
                        )}))
# joindt <- joindt[-(1:35),] # remove since they are A MESS!!!!!

# joindt <- joindt[RANDO=='YES']
joindt[,`:=`(
  age = fifelse(is.na(ENROLL),
                lubridate::year(lubridate::ymd(LUMBARDATE)),
                lubridate::year(lubridate::ymd(ENROLL))) - YOB,
  sex = as.numeric(SEX == 'M'),
  hiv_stat = fcase(
    HIV %in% c('POS', 'UNCERTAIN') , TRUE,
    HIV %in% 'NEG' , FALSE,
    (HIV %in% c('', 'NOT DONE', 'UNKNOWN') | is.na(HIV)) & ISHIV == 'C49488' , TRUE,
    (HIV %in% c('', 'NOT DONE', 'UNKNOWN') | is.na(HIV)) & ISHIV == 'N' , FALSE
  ),
  clin_headache = (ISHEADACHE == 'C49488'),
  clin_psychosis = (ISPSYCHOSIS == 'C49488'),
  clin_illness_day = as.numeric(ILLNESSDAY),
  clin_symptoms = ISWEIGHT | ISNSWEAT | ISCOUGH,
  clin_contact_tb = ISCONTUBER,
  clin_motor_palsy = HEMIPLEGIA | PARAPLEGIA | TETRAPLEGIA,
  clin_nerve_palsy = !NONE,
  clin_gcs = fifelse(is.na(GLASCOW), fifelse(!ISREDUCED,15, GCSV+GCSE+GCSM), as.numeric(GLASCOW)),
  csf_clear = APP_A == 'CLEAR',
  csf_vol = CSF_VOL,
  csf_rbc = as.numeric(REDCELL),
  csf_wbc = as.numeric(WHITECELL),
  csf_lym_pct = fifelse(is.na(LYMPER), WHITECELL/WHITECELL, LYMPER/100),
  csf_lympho = fifelse(is.na(LYMPER), WHITECELL, LYMPER/100*WHITECELL),
  csf_eos = fifelse(is.na(EOSPER), 0, EOSPER/100*WHITECELL),
  csf_protein = PROTEIN,
  csf_glucose = CSFGLU,
  csf_lactate = CSFLAC,
  glucose_ratio = CSFGLU/BLDGLU,
  csf_gram = (GRAM=="POS")%in%TRUE,
  csf_ink = (INDIAINK=="POS")%in% TRUE,
  csf_crypto = (CRYTO=="POS")%in%TRUE,
  # img_hydro = MRIRESULT=='Hydrocephalus',
  # img_basal = MRIRESULT == 'Basal meningeal enhancement',
  # img_tuber = MRIRESULT == 'Tuberculoma',
  # img_infarct = MRIRESULT == 'Infarction',
  # img_precon = MRIRESULT == 'Basal hyperintensity',
  xray_miliary_tb = MTB,
  xray_pul_tb = PTB,
  #csf_smear = fcase(ZN == 'POS', TRUE, 
  #                  ZN == 'NEG', FALSE,
  #                  default = NA),
  csf_smear = ZN,
  csf_xpert = XpertResult,
  #csf_xpert = fcase(TBNAAT == 'POS' | XpertResult, TRUE, 
  #                  TBNAAT == 'NEG' | !XpertResult, FALSE,
  #                  default = NA),
  csf_mgit = MycoResult,
  other_dx = !is.na(BACSPE) | stringr::str_detect(tolower(OTH), "[dương,+]")| JEVIGM == "POS" | ADD_JEVIGM == "POS" |  DENGUEIMG == "POS" | DENGUEPCR == "POS" | JEVCSF == "POS" | 
    stringr::str_detect(tolower(VIRO_OTH), "[dương,+]") | NMDAR == "POS",
  other_dis_dx_conf = !DISDIA %in% c("DIA1", "DIA10", "DIA16", "DIA2", "DIA4"),
  other_dis_dx = !DISDIA %in% c("DIA1", "DIA10"),
  tbm_dx = ifelse(is.na(DISDIA), NA, (DISDIA %in% c("DIA1", "DIA10"))) | (ISTBTREAT %in% 'C49488' & (as.Date(TBTREATDATE) > as.Date(LUMBARDATE)) %in% TRUE)
)][, `:=`(
  clin_score = pmin(4*clin_illness_day+2*as.numeric(clin_symptoms)+(2*clin_contact_tb)+1*clin_motor_palsy+1*clin_nerve_palsy+1*(clin_gcs<15),6),
  csf_score = pmin(4, csf_clear + (csf_wbc >= 10 & csf_wbc <= 500) + (csf_lym_pct > .5) + (csf_protein > 1) + (glucose_ratio < .5 | csf_glucose < 2.2)),
  # img_score = img_hydro+2*img_basal+2*img_tuber+img_infarct+2*img_precon,
  img_score = 0,
  tube_score = 2*xray_pul_tb+4*xray_miliary_tb
)][, crude_total_score := rowSums(.SD[, .(pmin(clin_score,6), pmin(csf_score,4), pmin(img_score,6), pmin(tube_score,4))], na.rm=T)]

# joindt <- joindt[USUBJID!='003-335'&!is.na(RANDO)]

summarise_all(joindt, ~is.na(.x) %>% sum) %>% t -> ms

maindt <- joindt[!is.na(WHITECELL), csf_rbc := ifelse(is.na(csf_rbc), 0, csf_rbc)]
saveRDS(maindt, file='data/cleaned/data_19EI_3.RDS')

