# Data processing for 05TB data
# Author: trinhdhk
# Ver 0.1.2004

library(data.table)
library(dplyr)
rm(list=ls())
source('r/fn_def.r')
data_files <- list.files('data/raw', pattern = '^\\w')

# Get data names
clinical_data_file <-  grep('DATA final', data_files, value = TRUE)
lab_data_file <- grep('lab', data_files, value = TRUE)
lab_data_2 <- grep('CSFResults05TB', data_files, value=TRUE)

# Load clinical data
full_sheet <- readxl::excel_sheets(file.path('data','raw', clinical_data_file))
for (sheet in full_sheet) assign(tolower(sheet), as.data.table(readxl::read_excel(file.path('data', 'raw', clinical_data_file), sheet=sheet)))

#Load lab data
lab_data <- readRDS(file.path('data', 'raw', lab_data_file))$tables
for (tbl in names(lab_data)) assign(paste('lab', tolower(tbl), sep='_'), as.data.table(lab_data[[tbl]]))
rm(lab_data)

lab_data2 <- as.data.table(readxl::read_excel(file.path('data', 'raw', lab_data_2)))
lab_data2[, ParNo:=substring(StudyCode, 5)]

# Get clinical & csf info from screening
maindt <-
  screx[,
        list(
          ParNo,
          clin_illness_day = DayIlless,
          clin_symptoms = Cough2W == 1 | WeightLoss == 1 | NightSweat == 1,
          clin_contact_tb = ContactTB == 1,
          clin_motor_palsy = Hemi == 1 | Para == 1 | Quadri == 1 ,
          clin_nerve_palsy = !CNPNone,
          scr_gcs = Gcs,
          csf_clear = CsfApp == 1
        )
        ][, Hospital := ifelse(strtrim(ParNo,1) %in% c('1','2'), 'PNT', 'HTD')] %>% 
  merge(gcs[WeekMonth == 1, list(ParNo, enr_gcs = GCS)], all.x = TRUE)

#update csf_clear from new lab_data2
new_csf <- lab_data2[, .(new_csf_clear = any(CSF_Appearance==1%in%TRUE)), by=ParNo]
maindt <- merge(maindt, new_csf, all.x=TRUE)
maindt[, csf_clear:= csf_clear %in% TRUE | new_csf_clear %in% TRUE, by = ParNo]
setkey(maindt, ParNo)

#get HIV status
maindt <- merge(maindt, hivdr[, list(ParNo, HIVStat)], all.x=TRUE)
setnames(maindt, 'HIVStat', 'hiv_stat')

# Get GCS
maindt[, clin_gcs := ifelse(is.na(enr_gcs), scr_gcs, enr_gcs)][, c('scr_gcs', 'enr_gcs'):=NULL]

# Blood count nearest the date of CSF take
hae_merge <- merge(hae, 
                   csf_hae[, head(.SD, 1), .SDcols = c('DateCsfH', 'CsfNeutro', 'CsfLympho', 'CsfBC', 'CsfRC'), by=ParNo],
                   all.y=TRUE)
hae_merge[, deltaD:=abs(difftime(DateHema, DateCsfH, units='days'))]
hae_merge <- unique(hae_merge[,
                              c('WBC', 'Neutro', 'Lympho') := 
                                ifelse(is.na(deltaD), rep(NA, 3),
                                       list(get_nearest(WBC, deltaD), get_nearest(Neutro, deltaD), get_nearest(Lympho, deltaD))),
                              by=ParNo][, list(ParNo, WBC, Neutro, Lympho, CsfBC, CsfRC, CsfNeutro, CsfLympho)])

maindt <- merge(maindt,
                hae_merge[, list(ParNo, 
                                 bld_wbc = WBC, bld_neu_pct = ifelse(is.na(Neutro), 100-Lympho, Neutro), bld_lym_pct = Lympho,
                                 csf_wbc = CsfBC, csf_neu_pct =  ifelse(is.na(CsfNeutro), 100-CsfLympho, CsfNeutro),
                                 csf_lym_pct = CsfLympho, csf_rbc = CsfRC)][, head(.SD,1), by=ParNo],
                all.x=TRUE)
maindt[, csf_lympho := csf_lym_pct*csf_wbc/100]

# Biomarkers nearest to the date of CSF

bio_merge <- merge(bio,
                   csf_bio[, head(.SD, 1), .SDcols = c('DateCsfB', 'CsfGlucose', 'CsfProtein', 'CsfLactate'), by=ParNo],
                   all.y=TRUE)
bio_merge[,deltaD:=abs(difftime(DateBio, DateCsfB, units='days'))]
bio_merge <- unique(bio_merge[,BlGlucose:=ifelse(is.na(deltaD), NA, get_nearest(BlGlucose, deltaD)),by=ParNo][, list(ParNo, BlGlucose, CsfGlucose, CsfProtein, CsfLactate)])

maindt <- merge(maindt, 
                bio_merge[, list(ParNo, bld_glucose = BlGlucose, csf_glucose = CsfGlucose,
                                 csf_protein = CsfProtein, csf_lactate = CsfLactate)][, head(.SD,1), by=ParNo],
                all.x=TRUE)
maindt[, glucose_ratio := csf_glucose/bld_glucose]

#imaging

ct_res <- 
  addinv_ct[, lapply(.SD, any),
            .SDcols = c('CtAbnHydro', 'CtAbnBasal', 'CtAbnTuber', 'CtAbnInfarct', 'CtAbnPreCon'),
            by = ParNo]

mri_res <- 
  addinv_mri[, lapply(.SD, any), 
             .SDcols = c('MriAbnHydro', 'MriAbnBasal', 'MriAbnTuber', 'MriAbnInfarct', 'MriAbnPreCon'),
             by = ParNo]

image_res <- merge(ct_res, mri_res, all=TRUE)[, lapply(.SD, function(x) ifelse(is.na(x), FALSE, x))]
image_res <- image_res[,.(
  img_hydro = CtAbnHydro | MriAbnHydro,
  img_basal = CtAbnBasal | MriAbnBasal,
  img_tuber = CtAbnTuber | MriAbnTuber,
  img_infarct = CtAbnInfarct | MriAbnInfarct,
  img_precon = CtAbnPreCon | MriAbnPreCon
), by= ParNo]

image_res[, img_score := img_hydro + 2*img_basal + 2*img_tuber + img_infarct + 2*img_precon]

# Add imaging result to main data
maindt <- merge(maindt, image_res, all.x=TRUE)

#evidence of TB elsewhere
chest_xray <- screx[, .(xray=XrayResult,
                        xray_miliary_tb=XrayResult==2,
                        xray_pul_tb=XrayResult==3), by=ParNo]
ultrasound <- addinv_us[, .(ultrasound=any(EviCns == 1)%in%TRUE), by=ParNo]
afb_othersource <- 
  plyr::join_all(
    list(
      ## Clinical database
      cult_blood[, .(blood = any(BlCulture == 1) %in% TRUE | any(BlZnStain == 1) %in% TRUE), by = ParNo],
      cult_urine[, .(urine = any(UriCulture == 1) %in% TRUE| any(UrZnStain == 1) %in% TRUE), by = ParNo],
      cult_lymph[, .(lymph = any(LnaCulture == 1) %in% TRUE| any(LnaZnStain == 1) %in% TRUE), by = ParNo],
      cult_others[, .(other = any(OthersCulture == 1) %in% TRUE | any(OthersZnStain == 1) %in% TRUE), by = ParNo],
      
      ## Lab database
      lab_sputum[, .(sputum = 
                       any(SpuZNstain == 1) %in% TRUE | any(SpuZNstain2 == 1) %in% TRUE | any(SpuZNstain3 == 1) %in% TRUE |
                       any(SpuZNstainCheck == 1) %in% TRUE | any(SpuZNstain2Check == 1) %in% TRUE | 
                       any(SpuZNstain3Check == 1) %in% TRUE | 
                       any(SpuMGITCult == 1) %in% TRUE | any(SpuGeneXpert %in% 1:4) %in% TRUE), by=ParNo],
      lab_gastricfluid[, .(gastric = 
                             any(GFZNstain == 1) %in% TRUE | any(GFZNstain2 == 1) %in% TRUE | any(GFZNstain3 == 1) %in% TRUE |
                             any(GFZNstainCheck == 1) %in% TRUE | any(GFZNstain2Check == 1) %in% TRUE | 
                             any(GFZNstain3Check == 1) %in% TRUE | 
                             any(GFMGITCult == 1) %in% TRUE | any(GFGeneXpert %in% 1:4) %in% TRUE), by=ParNo],
      lab_othersamples[, .(lab_other = any(LJCult == 1) %in% TRUE), by=ParNo]
    ),
    by = 'ParNo', type = 'full'
  )
tube <- plyr::join_all(list(chest_xray, ultrasound, afb_othersource), 
                       type='full', by='ParNo')

maindt <- merge(maindt, tube, all.x=TRUE, by='ParNo')

# Confirmed other diags
other_diag <- outc[, .(confirm_other = REASONOTH %in% 1:4), by=ParNo]
maindt <- merge(maindt, other_diag, all.x=TRUE)


# Indicators
csf_lab_bio <- lab_csf[, .(csf_lab_positive = 
                             any(csfZNstain == 1) %in% TRUE | any(csfZNstainCheck == 1) %in% TRUE |
                             any(csfMGITCult == 1) %in% TRUE | any(csfGeneXpert  %in% 1:4) %in% TRUE),
                       by=ParNo]
csf_lab_bio_details <-
  lab_csf[, .(csf_smear = any(csfZNstain == 1) %in% TRUE | any(csfZNstainCheck == 1) %in% TRUE,
              csf_mgit = any(csfMGITCult == 1) %in% TRUE,
              csf_xpert = any(csfGeneXpert  %in% 1:4) %in% TRUE), by=ParNo]
maindt <- merge(maindt, csf_lab_bio, all.x=TRUE, by='ParNo')
maindt <- merge(maindt, csf_lab_bio_details, all.x=TRUE, by='ParNo')

# Get consensus diag
maindt[,
       c('clin_score', 'csf_score','tube_score') := 
         list(4*(clin_illness_day>5) + 2*clin_symptoms + 2*clin_contact_tb + 
                clin_motor_palsy + clin_nerve_palsy + (clin_gcs<15),
              
              csf_clear %in% TRUE + (csf_wbc>=100 & csf_wbc <= 500) %in% TRUE + (csf_lym_pct>50) %in% TRUE + (csf_protein>1) %in% TRUE + (glucose_ratio<50|csf_glucose<2.2) %in% TRUE,
              
              case_when(xray == 3 ~ 2, xray == 2 ~ 4, TRUE ~ 0) +
                2*(ultrasound%in%TRUE) + 4*((blood | urine | lymph | other | sputum | gastric | lab_other) %in% TRUE)
         ), by=ParNo]
maindt[, csf_clear := tidyr::replace_na(csf_clear, FALSE)][, total_score := sum(min(clin_score, 6), min(csf_score, 4), min(tube_score, 4)), by = ParNo][, total_score := sum(total_score, ifelse(is.na(img_score), 0, min(img_score, 6))), by = ParNo]

maindt[, diag := 
         case_when(
           total_score + (2 * is.na(img_score)) >= 12 ~ 'Probable',
           total_score >= 6 ~ 'Possible',
           !is.na(total_score) ~ 'Unlikely'
         )][,diag_abs := case_when(
           csf_lab_positive %in% TRUE ~'Definite',
           confirm_other %in% TRUE ~ 'Definite_Not',
           TRUE ~ diag)][,diag_abs:=factor(diag_abs, levels = c('Definite_Not', 'Unlikely', 'Possible', 'Probable', 'Definite'))][,diag:=factor(diag, levels=c('Definite_Not', 'Unlikely', 'Possible', 'Probable'))]


#Cleanup
maindt <- maindt[diag_abs != 'Definite_Not'
                 ][,
                   c('lymph', 'other', 'new_csf_clear', 'confirm_other' , 'xray', 'urine'):=NULL]


saveRDS(maindt, file='data/cleaned/data_05TB.RDS')
