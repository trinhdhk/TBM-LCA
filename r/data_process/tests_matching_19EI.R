unmatched_xpert = joindt[(((!is.na(TBNAAT) & TBNAAT != "NOT DONE") & is.na(XpertResult)) | 
                           ((TBNAAT == "NOT DONE") & !is.na(XpertResult))), 
                         .(USUBJID, ResultFrom19EIDat = TBNAAT, ResultfromTBGroup = XpertResult)]
unmatched_xpert = merge(data$INEX[,.(USUBJID, INITIAL, HOSP_NB, SEX)], unmatched_xpert, all.y = TRUE)
write.csv(unmatched_xpert, "data/unmatched_xpert.csv", row.names = F)
write.csv(filter(m2, !((abs(as.Date(LUMBARDATE)-DateSample) < 2)%in%T) & ZN1 != "NOT DONE" & !(((ZN1 == "POS") == ZN) %in% T)), "mismatchresdate.csv", na = "")
write.csv(filter(m2, !((abs(as.Date(LUMBARDATE)-DateSample) < 2)%in%T) & ZN1 != "NOT DONE"), "mismatchdate.csv", na = "")
write.csv(filter(m2, ((abs(as.Date(LUMBARDATE)-DateSample) < 2)%in%T) & ZN1 != "NOT DONE" & !(((ZN1 == "POS") == ZN) %in% T)), "mismatchres.csv", na = "")