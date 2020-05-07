# Data process for 26TB
# Author: trinhdhk
# Version: 2005.1

rm(list=ls())
library(data.table)
data.folder <- 'data/raw/2627TB/'
data.file <- grep('_26TB', dir(data.folder), value=TRUE)[1]
.dt <- import.Excel(file.path(data.folder,data.file))
for (i in .dt$tableNames) assign(tolower(i), .dt$tables[[i]])
rm(.dt)

setDT(base)
setDT(scr)
setDT(stdr_stdrug)
setDT(basecsf)
setDT(basecsf_jkt)
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
    clin_illness_day = max(HEADACHEDAY, IRRITABILITYDAY, VOMITDAY, FEVERDAY, NECKSTIFFDAY, SEIZURESDAY, NEURODAY, CONSCIOUSDAY, LETHARGYDAY, 0, na.rm=TRUE),
    clin_symptoms = COUGH=='Y'|WEIGHTLOSS=='Y'|NIGHTSWEATS=='Y',
    clin_contact_tb = LIVELUNGTB=='Y',
    clin_motor_palsy = HEMIPLEGIA=='Y'|PARAPLEGIA=='Y'|TETRAPLEGIA=='Y',
    clin_nerve_palys = CNP=='Y',
    clin_gcs = dplyr::case_when(
        !is.na(GCS) ~ as.numeric(GCS),
        ALTEREDCONSCIOUS != 'Y' ~ 15
    )
), by=USUBJID][, 
    clin_score := min(6,sum(4*(clin_illness_day>5), 2*clin_symptoms, clin_contact_tb, clin_motor_palsy, (clin_gcs<15), na.rm=TRUE)),
    by=USUBJID]

# CSF score
.basecsf <- dplyr::bind_rows(basecsf, basecsf_jkt)
# .basecsf <- merge(.basecsf, base[,.(USUBJID, SITEID)], by='USUBJID', all.x=TRUE)
csf <- .basecsf[, .(
    csf_clear = APPEARANCE == '1',
    csf_wbc = CSFWBC,
    csf_lym_pct = CSFLYMLE/100,
    csf_protein = PROTEIN*ifelse(SITEID=='044',.01,1),
    csf_glucose = CSFGLUC*ifelse(SITEID=='044',.0555,1),
    bld_glucose = PAIREDGLUC*ifelse(SITEID=='044',.0555,1)), by=USUBJID][, glucose_ratio:=csf_glucose/bld_glucose][,
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
confirmed <- .basecsf[, .(
    csf_smear = any(ZNSMEAR=='1'),
    csf_mgit = any(TBNAAT=='1'),
    csf_xpert = any(MYCORESULT=='1')
), by=USUBJID]

maindt <- purrr::reduce(list(demoinfo, clin, csf, tube, confirmed), merge, all.x=TRUE, by='USUBJID')
maindt$hiv_stat <- 1

saveRDS(maindt, 'data/cleaned/data_26TB.RDS')