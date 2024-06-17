source("/data1/sanghyeon/Projects/mr_drug_repurposing/src/PRS/PRS_calculation/phenotype_and_covariate/v2/fnc.query_date.R")

icd10 <- c("Q874", "Q796", "Q61.1", "Q61.2", "Q61.3")


start <- Sys.time()
df.icd10 <- query_date(code_type="ICD10", code=icd10, simplify=TRUE)
cat(nrow(df.icd10 %>% filter(!is.na(date_earliest.ICD10))), "\n")
cat("Elapsed time: ", Sys.time() - start, "\n")

df.IA_control_exclude <- df.icd10 %>%
    mutate(IA_control_exclude = ifelse(is.na(date_earliest.ICD10), 0, 1),
        FID = f.eid, IID = f.eid) %>%
    dplyr::select(FID, IID, date_earliest.ICD10, IA_control_exclude)

cat(nrow(df.IA_control_exclude %>% filter(IA_control_exclude == 1)), "\n")

write.table(df.IA_control_exclude,
            "UKBB.IA_IA_control_exclude.txt",
            quote=F, row.names=F, sep=" "
            )
