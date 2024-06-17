source("/data1/sanghyeon/Projects/mr_drug_repurposing/src/PRS/PRS_calculation/phenotype_and_covariate/v2/fnc.query_date.R")

icd10 <- c("I67.1", "I60")


start <- Sys.time()
df.icd10 <- query_date(code_type="ICD10", code=icd10, simplify=TRUE)
cat(nrow(df.icd10 %>% filter(!is.na(date_earliest.ICD10))), "\n")
cat("Elapsed time: ", Sys.time() - start, "\n")

df.IA <- df.icd10 %>%
    mutate(IA = ifelse(is.na(date_earliest.ICD10), 0, 1),
        FID = f.eid, IID = f.eid) %>%
    dplyr::select(FID, IID, date_earliest.ICD10, IA)

cat(nrow(df.IA %>% filter(IA == 1)), "\n")

write.table(df.IA,
            "UKBB.IA_case.txt",
            quote=F, row.names=F, sep=" "
            )
