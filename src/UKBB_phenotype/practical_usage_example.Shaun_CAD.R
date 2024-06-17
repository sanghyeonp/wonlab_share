source("fnc.query_date.R")

icd10 <- c("I21", "I22", "I23", "I24", "I25.2", "I46")
icd9 <- c("410", "411", "412", "414.2", "427.5")
opcs4 <- c("K40", "K41", "K42", "K44", "K45", "K46", "K47.1", "K49", "K50", "K75")
self_report <- "1075"

df.icd10 <- query_date(code_type="ICD10", code=icd10, simplify=TRUE)
df.icd9 <- query_date(code_type="ICD9", code=icd9, simplify=TRUE)
df.opcs4 <- query_date(code_type="OPCS4", code=opcs4, simplify=TRUE)
df.self_report <- query_date(code_type="Self-report", code=self_report, self_report_instance=0, simplify=TRUE)


df.final <- query_date_multiple_source(queried_code_type=c("ICD10", "ICD9", "OPCS4", "Self-report"),
                                    ICD10 = df.icd10, ICD9 = df.icd9,
                                    OPCS4 = df.opcs4, Self_report = df.self_report)

