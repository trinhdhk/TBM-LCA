# Data process for 27TB
# Author: trinhdhk
# Version: 2005.1

rm(list=ls())
library(data.table)
source('r/fn_def.R')
data.folder <- 'data/raw/2627TB/'
data.file <- grep('_27TB', dir(data.folder), value=TRUE)[1]
.dt <- import.Excel(file.path(data.folder,data.file))
for (i in .dt$tableNames) assign(tolower(i), .dt$tables[[i]])
rm(.dt)

setDT(base)
setDT(scr)
setDT(stdr_stdrug)
setDT(basecsf)
setDT(xraye)
setDT(baselaboth)

# Demographic info
.base <- scr[ENROLLED=='Y', .(USUBJID, DATERANDOM, DATEICF, BL_DATE_BASE=as.Date(DATERANDOM))]
.datedrug <- stdr_stdrug[STDR_STDRUG_SEQ==1,.(USUBJID, BL_DATE_DR=as.Date(DATESTART))]
.demo <- base[,.(USUBJID, BIRTHYR, SEX, DATEADM, BL_DATE_DEMO=as.Date(DATEADM))]
.tmp <- purrr::reduce(list(.base,.datedrug,.demo), merge, all.x=TRUE,by='USUBJID')
.tmp[, BL_DATE:=dplyr::case_when(!is.na(BL_DATE_DR)~BL_DATE_DR,!is.na(BL_DATE_DEMO)~BL_DATE_DEMO,TRUE~BL_DATE_BASE)]
.tmp[, age:=lubridate::year(BL_DATE)-BIRTHYR]
demoinfo <- .tmp[,.(USUBJID, age, sex=SEX=='M')]

# CLinical score
clin <- base[, .(
    ISDIABETES = DIABETES == "Y",
    clin_illness_day = max(HEADACHEDAY, IRRITABILITYDAY, VOMITDAY, FEVERDAY, NECKSTIFFDAY, SEIZURESDAY, NEURODAY, CONSCIOUSDAY, LETHARGYDAY, 0, na.rm=TRUE),
    clin_symptoms = COUGH=='Y'|WEIGHTLOSS=='Y'|NIGHTSWEATS=='Y',
    clin_contact_tb = LIVELUNGTB=='Y' %in% TRUE,
    ISCOUGH = fcase(COUGH == "Y", TRUE, COUGH == "N", FALSE, default = NA),
    ISNSWEAT = fcase(NIGHTSWEATS == "Y", TRUE, NIGHTSWEATS == "N", FALSE, default = NA),
    ISWEIGHT = fcase(WEIGHTLOSS == "Y", TRUE, WEIGHTLOSS == "N", FALSE, default = NA),
    HEMIPLEGIA = fcase(HEMIPLEGIA == "Y", TRUE, 
                       HEMIPLEGIA == "N", FALSE,
                       default = NA),
    PARAPLEGIA = fcase(PARAPLEGIA == "Y", TRUE, 
                       PARAPLEGIA == "N", FALSE,
                       default = NA),
    TETRAPLEGIA = fcase(TETRAPLEGIA == "Y", TRUE, 
                        TETRAPLEGIA == "N", FALSE,
                        default = NA),
    clin_nerve_palsy = CNP=='Y',
    GCSE, GCSM, GCSV,
    clin_gcs = dplyr::case_when(
        !is.na(GCS) ~ as.numeric(GCS),
        ALTEREDCONSCIOUS != 'Y' ~ 15
    )
), by=USUBJID][,`:=`(clin_symptoms = ISCOUGH|ISNSWEAT|ISWEIGHT,
                     clin_motor_palsy = HEMIPLEGIA|PARAPLEGIA|TETRAPLEGIA)][, 
    clin_score := min(6,sum(4*(clin_illness_day>5), 2*clin_symptoms, clin_contact_tb, clin_motor_palsy, (clin_gcs<15), na.rm=TRUE)),
    by=USUBJID]

# CSF score
csf <- basecsf[, .(
    csf_clear = APPEARANCE == '1',
    csf_wbc = CSFWBC,
    csf_lym_pct = CSFLYMLE/100,
    csf_lympho = CSFLYMLE*CSFWBC/100,
    csf_protein = PROTEIN,
    csf_glucose = CSFGLUC,
    csf_lactate = LACTATE,
    bld_glucose = PAIREDGLUC), by=USUBJID][,`:=`(
        glucose_ratio = csf_glucose/bld_glucose,
        csf_lympho = csf_lym_pct*csf_wbc)][,
    csf_score := min(sum(csf_clear, csf_wbc>=10&csf_wbc<=500, csf_lym_pct>.5, csf_protein>1, glucose_ratio<.5|csf_glucose<2.2, na.rm=TRUE),4), 
    by=USUBJID]

# Evidence elsewhere
.tube <- xraye[VISITNUM=='BASE', .(
    xray_miliary_tb = XRAYRESULT=='2',
    xray_pul_tb = XRAYRESULT=='3'
), by=USUBJID]
.other <- baselaboth[,other_test:=ZNSMEAR=='1'|MYCORESULT=='1'][!is.na(other_test), .(tube3=max(other_test)), by=USUBJID]
tube <- merge(.tube, .other, all.x=TRUE, by='USUBJID')
tube[, tube_score := min(4,sum(4*xray_miliary_tb,2*xray_pul_tb,4*tube3, na.rm=TRUE)),by=USUBJID]

# Confirmation tests
confirmed <- basecsf[, .(
    csf_smear = any(ZNSMEAR=='1'),
    csf_xpert = any(TBNAAT=='1'),
    csf_mgit = any(MYCORESULT=='1')
), by=USUBJID]

maindt <- purrr::reduce(list(demoinfo, clin, csf, tube, confirmed), merge, all.x=TRUE, by='USUBJID')
maindt$hiv_stat <- 0

saveRDS(maindt, 'data/cleaned/data_27TB.RDS')
