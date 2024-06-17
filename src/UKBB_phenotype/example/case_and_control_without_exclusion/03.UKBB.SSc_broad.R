source("/data1/sanghyeon/Projects/mr_drug_repurposing/src/PRS/PRS_calculation/phenotype_and_covariate/v2/fnc.query_date.R")

icd10 <- c("M34", "L94")
icd9 <- c("710.1", "517.2")
# opcs4 <- c("K40", "K41", "K42", "K44", "K45", "K46", "K47.1", "K49", "K50", "K75")
self_report <- "1384"

start <- Sys.time()
df.icd10 <- query_date(code_type="ICD10", code=icd10, simplify=TRUE)
cat(nrow(df.icd10 %>% filter(!is.na(date_earliest.ICD10))), "\n")
cat("Elapsed time: ", Sys.time() - start, "\n")

start <- Sys.time()
df.icd9 <- query_date(code_type="ICD9", code=icd9, simplify=TRUE)
cat(nrow(df.icd9 %>% filter(!is.na(date_earliest.ICD9))), "\n")
cat("Elapsed time: ", Sys.time() - start, "\n")

# df.opcs4 <- query_date(code_type="OPCS4", code=opcs4, simplify=TRUE)

start <- Sys.time()
df.self_report <- query_date(code_type="Self-report", code=self_report, self_report_instance="all", simplify=TRUE)
cat(nrow(df.self_report %>% filter(!is.na(date_earliest.selfreport_all))), "\n")
cat("Elapsed time: ", Sys.time() - start, "\n")

cat("Merging...")
start <- Sys.time()
df.final <- query_date_multiple_source(queried_code_type=c("ICD10", "ICD9", "Self-report"),
                                    ICD10 = df.icd10, ICD9 = df.icd9,
                                    Self_report = df.self_report)
cat("Elapsed time: ", Sys.time() - start, "\n")

df.SSc_broad <- df.final %>%
    mutate(SSc_broad = ifelse(is.na(date_earliest), 0, 1),
        FID = f.eid, IID = f.eid) %>%
    dplyr::select(FID, IID, date_earliest, SSc_broad)

cat(nrow(df.SSc_broad %>% filter(SSc_broad == 1)), "\n")

write.table(df.SSc_broad,
            "UKBB.SSc_broad.txt",
            quote=F, row.names=F, sep=" "
            )
