# Data process code for 19EI
# Author: Trinhdhk
# Version: 2101.1

library(data.table); library(dplyr)
data.folder <- 'data/raw/19EI'
inex <- fread(file.path(data.folder,'19EI_INEX.csv'))
base <- fread(file.path(data.folder,'19EI_BASE.csv'))
baselp <- fread(file.path(data.folder,'19EI_BASE_LP.csv'))
imag <- fread(file.path(data.folder,'19EI_IMAG.csv'))
addinfo <- fread(file.path(data.folder,'19EI_ADDINFO.csv'))

joindt <- plyr::join_all(
  list(
    inex[, .(USUBJID, YOB, SEX, ENROLL)],
    base[, .(USUBJID, ISHIV, ISCHRONIC, ISDIABETE, ISFEVER, ISCONTUBER, ILLNESSDAY, 
             ISREDUCED, ISWEIGHT, ISNSWEAT, ISCOUGH, NONE, ISHEADACHE, ISLOCALSEIZURE, ISGENNERALSEIZURE, 
             ISPSYCHOSIS, ISLANGCHANGE, ISMOVEMENT,
             HEMIPLEGIA, PARAPLEGIA, TETRAPLEGIA, 
             GENCONVUL, LOCALCONVUL, 
             GLASCOW, GCSE, GCSM, GCSV)],
    baselp[, .(USUBJID, LUMBARDATE, APP_A, 
               CSF_VOL= fcase(!is.na(as.numeric(TB_VOL)) & (TB_VOL != "" %in% T), as.numeric(TB_VOL),
                              !is.na(as.numeric(CSFVOLUME)) & (CSFVOLUME != "" %in% T), as.numeric(CSFVOLUME),
                              !is.na(as.numeric(CSFVOL)) & (CSFVOL != "" %in% T), as.numeric(CSFVOL),
                              default=NA),
               REDCELL, WHITECELL, NEUPER, LYMPER, EOSPER,
               PROTEIN, CSFGLU, BLDGLU, CSFLAC, ZN, TBNAAT, RANDO, NAATSPEC, 
               MYCORESULT, LUNGPOS, REP_NAAT, REP_RES, GRAM, BACCUL)],
    imag[, .(USUBJID, PTB, MTB, MRIRESULT, CTRESULT)],
    addinfo[, .(USUBJID, HIV)]
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
joindt <- joindt[-(1:35),] # remove since they are A MESS!!!!!

# joindt <- joindt[RANDO=='YES']
joindt[,`:=`(
  age = if_else(is.na(ENROLL),
                lubridate::year(lubridate::mdy_hm(joindt$LUMBARDATE)),
                lubridate::year(lubridate::dmy(joindt$ENROLL))) - YOB,
  sex = as.numeric(SEX == 'M'),
  hiv_stat = dplyr::case_when(
    HIV %in% c('POS', 'UNCERTAIN') ~ TRUE,
    HIV %in% 'NEG' ~ FALSE,
    (HIV %in% c('', 'NOT DONE', 'UNKNOWN') | is.na(HIV)) & ISHIV == 'C49488' ~ TRUE,
    (HIV %in% c('', 'NOT DONE', 'UNKNOWN') | is.na(HIV)) & ISHIV == 'N' ~ FALSE,
  ),
  clin_illness_day = as.numeric(ILLNESSDAY),
  clin_symptoms = ISWEIGHT | ISNSWEAT | ISCOUGH,
  clin_contact_tb = ISCONTUBER,
  clin_motor_palsy = HEMIPLEGIA | PARAPLEGIA | TETRAPLEGIA,
  clin_nerve_palsy = !NONE,
  clin_gcs = ifelse(is.na(GLASCOW), ifelse(!ISREDUCED,15,NA), as.numeric(GLASCOW)),
  csf_clear = APP_A == 'CLEAR',
  csf_vol = CSF_VOL,
  csf_rbc = as.numeric(REDCELL),
  csf_wbc = as.numeric(WHITECELL),
  csf_lym_pct = ifelse(is.na(LYMPER), WHITECELL/WHITECELL, LYMPER/100),
  csf_lympho = ifelse(is.na(LYMPER), WHITECELL, LYMPER/100*WHITECELL),
  csf_protein = PROTEIN,
  csf_glucose = CSFGLU,
  csf_lactate = CSFLAC,
  glucose_ratio = CSFGLU/BLDGLU,
  img_hydro = MRIRESULT=='Hydrocephalus',
  img_basal = MRIRESULT == 'Basal meningeal enhancement',
  img_tuber = MRIRESULT == 'Tuberculoma',
  img_infarct = MRIRESULT == 'Infarction',
  img_precon = MRIRESULT == 'Basal hyperintensity',
  xray_miliary_tb = MTB,
  xray_pul_tb = PTB,
  csf_smear = ZN == 'POS',
  csf_xpert = fcase(TBNAAT == 'POS', TRUE, 
                    TBNAAT == 'NEG', FALSE,
                    default = NA),
  csf_mgit = MYCORESULT == 'POS'
)][, `:=`(
  clin_score = pmin(4*clin_illness_day+2*as.numeric(clin_symptoms)+(2*clin_contact_tb)+1*clin_motor_palsy+1*clin_nerve_palsy+1*(clin_gcs<15),6),
  csf_score = pmin(4, csf_clear + (csf_wbc >= 10 & csf_wbc <= 500) + (csf_lym_pct > .5) + (csf_protein > 1) + (glucose_ratio < .5 | csf_glucose < 2.2)),
  img_score = img_hydro+2*img_basal+2*img_tuber+img_infarct+2*img_precon,
  tube_score = 2*xray_pul_tb+4*xray_miliary_tb
)][, crude_total_score := rowSums(.SD[, .(clin_score, csf_score, img_score, tube_score)], na.rm=T)]

joindt <- joindt[USUBJID!='003-335'&!is.na(RANDO)]

summarise_all(joindt, ~is.na(.x) %>% sum) %>% t -> ms

maindt <- joindt[!is.na(WHITECELL), csf_rbc := ifelse(is.na(csf_rbc), 0, csf_rbc)]
saveRDS(maindt, file='data/cleaned/data_19EI_2.RDS')

