# Data process code for 19EI
# Author: Trinhdhk
# Version: 2005.1

library(data.table)
data.folder <- 'data/raw/19EI'
inex <- fread(file.path(data.folder,'19EI_INEX.csv'))
base <- fread(file.path(data.folder,'19EI_BASE.csv'))
baselp <- fread(file.path(data.folder,'19EI_BASE_LP.csv'))
imag <- fread(file.path(data.folder,'19EI_IMAG.csv'))
addinfo <- fread(file.path(data.folder,'19EI_ADDINFO.csv'))

joindt <- plyr::join_all(
    list(
        inex[, .(USUBJID, YOB, SEX, ENROLL)],
        base[, .(USUBJID, ISCONTUBER, ILLNESSDAY, 
            ISREDUCED, ISWEIGHT, ISNSWEAT, ISCOUGH, NONE,
            HEMIPLEGIA, PARAPLEGIA, TETRAPLEGIA, 
            GENCONVUL, LOCALCONVUL, 
            GLASCOW)],
        baselp[, .(USUBJID, APP_A, WHITECELL, LYMPER,
            PROTEIN, CSFGLU, BLDGLU, ZN, TBNAAT, RANDO, NAATSPEC, 
            MYCORESULT, LUNGPOS, REP_NAAT, REP_RES, GRAM, BACCUL)],
        imag[, .(USUBJID, PTB, MTB, MRIRESULT, CTRESULT)],
        addinfo[, .(USUBJID, HIV)]
    ), 
    by='USUBJID',
    type='left'
)

joindt <- joindt[-1:-35][RANDO=='YES']
maindt <- joindt[, .(
    age = 2019 - YOB,
    hiv_stat = HIV == 'POS',
    clin_illness_day = as.numeric(ILLNESSDAY),
    clin_symptoms = ISWEIGHT == 'C49488' | ISNSWEAT == 'C49488' | ISCOUGH == 'C49488',
    clin_contact_tb = ISCONTUBER == 'C49488',
    clin_motor_palsy = HEMIPLEGIA == 'C49488' | PARAPLEGIA == 'C49488' | TETRAPLEGIA == 'C49488',
    clin_nerve_palsy = !NONE,
    clin_gcs = ifelse(is.na(GLASCOW), ifelse(ISREDUCED!='C49488',15,NA), as.numeric(GLASCOW)),
    csf_clear = APP_A == 'CLEAR',
    csf_wbc = as.numeric(WHITECELL),
    csf_lym_pct = as.numeric(LYMPER)/100,
    csf_lympho = LYMPER/100*WHITECELL,
    csf_protein = PROTEIN,
    csf_glucose = CSFGLU,
    glucose_ratio = CSFGLU/BLDGLU,
    img_hydro = MRIRESULT=='Hydrocephalus',
    img_basal = MRIRESULT == 'Basal meningeal enhancement',
    img_tuber = MRIRESULT == 'Tuberculoma',
    img_infarct = MRIRESULT == 'Infarction',
    img_precon = MRIRESULT == 'Basal hyperintensity',
    xray_miliary_tb = MTB,
    xray_pul_tb = PTB,
    csf_smear = ZN == 'POS',
    csf_xpert = TBNAAT == 'POS',
    csf_mgit = MYCORESULT == 'POS'
)][, `:=`(
    clin_score = pmin(4*clin_illness_day+2*as.numeric(clin_symptoms)+(2*clin_contact_tb)+1*clin_motor_palsy+1*clin_nerve_palsy+1*(clin_gcs<15),6),
    csf_score = pmin(4, csf_clear + (csf_wbc >= 10 & csf_wbc <= 500) + (csf_lym_pct > .5) + (csf_protein > 1) + (glucose_ratio < .5 | csf_glucose < 2.2)),
    img_score = img_hydro+2*img_basal+2*img_tuber+img_infarct+2*img_precon,
    tube_score = 2*xray_pul_tb+4*xray_miliary_tb
)][, crude_total_score := rowSums(.SD[, .(clin_score, csf_score, img_score, tube_score)], na.rm=T)]

saveRDS(maindt, file='data/cleaned/data_19EI.RDS')

