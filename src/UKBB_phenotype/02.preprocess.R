# Sanghyeon Park
# 2024.06.06
# Credit: Shaun Chen

library(dplyr)
library(tidyr)

######################################
### Preprocess ICD10 code and date ###
######################################
# print("Preprocess ICD10 code and date")
# start_time <- Sys.time()
# df.icd10 <- readRDS("UKBB.41270.diagnoses_main_ICD10.rds")
# df.icd10_date <- readRDS("UKBB.41280.diagnoses_main_ICD10_date.rds")

# df.icd10.collapse <- as.data.frame(df.icd10 %>%
#         rowwise() %>%
#         mutate(ICD10_main = paste(c_across(-f.eid), collapse = ";")) %>%
#         ungroup() %>%
#         dplyr::select(f.eid, ICD10_main))

# idate_columns <- sapply(df.icd10_date, inherits, "IDate")
# df.icd10_date.collapse <- as.data.frame(df.icd10_date %>%
#     mutate(across(all_of(names(which(idate_columns))), as.character)) %>%
#     mutate_all(~replace_na(., "")) %>%
#     rowwise() %>%
#     mutate(ICD10_date = paste(c_across(-f.eid), collapse = ";")) %>%
#     ungroup() %>%
#     dplyr::select(f.eid, ICD10_date))

# df.icd10_code_date <- merge(df.icd10.collapse, df.icd10_date.collapse, by="f.eid", all=T)

# saveRDS(df.icd10_code_date, "UKBB.41270_41280.ICD10_code_and_date_merged.rds")

# end_time <- Sys.time()
# elapsed_time <- end_time - start_time
# print(elapsed_time)

### This is how to query a single ICD-10 code date
# ICD10_code <- "M161"
# df.icd10_date.quried <- as.data.frame(
#     df.icd10_code_date %>%
#     rowwise %>%
#     mutate(index = ifelse(grepl(ICD10_code, ICD10_main),
#                           which(strsplit(ICD10_main, ";")[[1]] == ICD10_code), NA),
#            date = ifelse(is.na(index), NA, strsplit(ICD10_date, ";")[[1]][index])) %>%
#     dplyr::select(f.eid, date)
#     )
### 

#####################################
### Preprocess ICD9 code and date ###
#####################################
# print("Preprocess ICD9 code and date")
# start_time <- Sys.time()
# df.icd9 <- readRDS("UKBB.41271.diagnoses_main_ICD9.rds")
# df.icd9_date <- readRDS("UKBB.41281.diagnoses_main_ICD9_date.rds")

# integer_columns <- sapply(df.icd9, inherits, "integer")
# df.icd9 <- df.icd9 %>%
#     mutate(across(all_of(names(which(integer_columns))), as.character)) %>%
#     mutate_all(~replace_na(., ""))

# df.icd9.collapse <- as.data.frame(df.icd9 %>%
#                        rowwise() %>%
#                        mutate(ICD9_main = paste(c_across(-f.eid), collapse = ";")) %>%
#                        ungroup() %>%
#                        dplyr::select(f.eid, ICD9_main))

# idate_columns <- sapply(df.icd9_date, inherits, "IDate")
# df.icd9_date.collapse <- as.data.frame(df.icd9_date %>%
#                             mutate(across(all_of(names(which(idate_columns))), as.character)) %>%
#                             mutate_all(~replace_na(., "")) %>%
#                             rowwise() %>%
#                             mutate(ICD9_date = paste(c_across(-f.eid), collapse = ";")) %>%
#                             ungroup() %>%
#                             dplyr::select(f.eid, ICD9_date))

# df.icd9_code_date <- merge(df.icd9.collapse, df.icd9_date.collapse, by="f.eid", all=T)

# saveRDS(df.icd9_code_date, "UKBB.41271_41281.ICD9_code_and_date_merged.rds")

# end_time <- Sys.time()
# elapsed_time <- end_time - start_time
# print(elapsed_time)

### This is how to query a single ICD-9 code date
# ICD9_code <- "3540"
# df.icd9_date.quried <- as.data.frame(
#     df.icd9_code_date %>%
#     rowwise %>%
#     mutate(index = ifelse(grepl(ICD9_code, ICD9_main),
#                           which(strsplit(ICD9_main, ";")[[1]] == ICD9_code), NA),
#            date = ifelse(is.na(index), NA, strsplit(ICD9_date, ";")[[1]][index])) %>%
#     dplyr::select(f.eid, date)
#     )
### 

######################################
### Preprocess OPCS4 code and date ###
######################################
# print("Preprocess OPCS4 code and date")
# start_time <- Sys.time()

# df.opcs4 <- readRDS("UKBB.41272.OPCS4.rds")
# df.opcs4_date <- readRDS("UKBB.41282.OPCS4_date.rds")

# df.opcs4.collapse <- as.data.frame(df.opcs4 %>%
#                        rowwise() %>%
#                        mutate(OPCS4 = paste(c_across(-f.eid), collapse = ";")) %>%
#                        ungroup() %>%
#                        dplyr::select(f.eid, OPCS4))

# idate_columns <- sapply(df.opcs4_date, inherits, "IDate")
# df.opcs4_date.collapse <- as.data.frame(df.opcs4_date %>%
#                             mutate(across(all_of(names(which(idate_columns))), as.character)) %>%
#                             mutate_all(~replace_na(., "")) %>%
#                             rowwise() %>%
#                             mutate(OPCS4_date = paste(c_across(-f.eid), collapse = ";")) %>%
#                             ungroup() %>%
#                             dplyr::select(f.eid, OPCS4_date))

# df.opcs4_code_date <- merge(df.opcs4.collapse, df.opcs4_date.collapse, by="f.eid", all=T)

# saveRDS(df.opcs4_code_date, "UKBB.41272_41282.OPCS4_code_and_date_merged.rds")

# end_time <- Sys.time()
# elapsed_time <- end_time - start_time
# print(elapsed_time)

### This is how to query a single OPCS4 code date
# OPCS4_code <- "E492"
# df.opcs4_date.quried <- as.data.frame(
#     df.opcs4_code_date %>%
#     rowwise %>%
#     mutate(index = ifelse(grepl(OPCS4_code, OPCS4), 
#                           which(strsplit(OPCS4, ";")[[1]] == OPCS4_code), NA),
#            date = ifelse(is.na(index), NA, strsplit(OPCS4_date, ";")[[1]][index])) %>%
#     dplyr::select(f.eid, date)
#     )
### 

##########################################
### Preprocess Self-report code and age/year ###
##########################################
print("Preprocess self-report code and age/year")
start_time <- Sys.time()

df.self_report <- readRDS("UKBB.20002.Self_report.rds")
df.self_report_age_year <- readRDS("UKBB.87.Self_report_age_year.rds")
df.year_of_birth <- readRDS("UKBB.34.year_of_birth.rds")

integer_columns <- sapply(df.self_report, inherits, "integer")
logical_columns <- sapply(df.self_report, inherits, "logical")
df.self_report <- df.self_report %>%
    mutate(across(all_of(names(which(logical_columns))), as.character)) %>%
    mutate(across(all_of(names(which(integer_columns))), as.character)) %>%
    mutate_all(~replace_na(., ""))

integer_columns <- sapply(df.self_report_age_year, inherits, "integer")
logical_columns <- sapply(df.self_report_age_year, inherits, "logical")
df.self_report_age_year <- df.self_report_age_year %>%
    mutate(across(all_of(names(which(logical_columns))), as.character)) %>%
    mutate(across(all_of(names(which(integer_columns))), as.character)) %>%
    mutate_all(~replace_na(., ""))

integer_columns <- sapply(df.year_of_birth, inherits, "integer")
df.year_of_birth <- df.year_of_birth %>%
    mutate(across(all_of(names(which(integer_columns))), as.character)) %>%
    mutate_all(~replace_na(., ""))

process_element <- function(element, year_of_birth) {
    # Overall, last day of the year is used to avoid overriding the date obtained from ICD codes and/or OPCS4
    # He/she had disease but don't know when
    if (element == "-1" | element == "-3") return ("3000-12-31")
    # If age of diagnosis is given, year of birth + age and last day of the year
    if (nchar(element) %in% c(1, 2)) return (paste(as.numeric(year_of_birth) + as.numeric(element), 12, 31, sep="-"))
    # If year of diagnosis is given, year of diagnosis and last day of the year
    if (nchar(element) == 4) return (paste(element, 12, 31, sep="-"))
    return (element)
}


## All instances
df.self_report.collapse <- as.data.frame(df.self_report %>%
                                         rowwise() %>%
                                         mutate(selfreport_all = paste(c_across(-f.eid), collapse = ";")) %>%
                                         ungroup() %>%
                                         dplyr::select(f.eid, selfreport_all))
df.self_report_age_year.collapse <- as.data.frame(df.self_report_age_year %>%
                                                  rowwise() %>%
                                                  mutate(selfreport_age_year_all = paste(c_across(-f.eid), collapse = ";")) %>%
                                                  ungroup() %>%
                                                  dplyr::select(f.eid, selfreport_age_year_all))

df.self_report_code_date <- merge(df.self_report.collapse, df.self_report_age_year.collapse, by="f.eid", all=T)
df.self_report_code_date <- merge(df.self_report_code_date, df.year_of_birth, by="f.eid", all.x=T)

df.self_report_date.collapse <- as.data.frame(df.self_report_code_date %>%
    rowwise() %>%
    mutate(selfreport_date_all = paste(sapply(strsplit(selfreport_age_year_all, ";")[[1]], 
                                               process_element, 
                                               year_of_birth=f.34.0.0), collapse=";")) %>%
    dplyr::select(f.eid, selfreport_all, selfreport_age_year_all, selfreport_date_all))

## First instance (Initial visit)
columns_of_instance0 <- c("f.eid", names(df.self_report)[grepl("f.20002.0.*", names(df.self_report))])
df.self_report.instance0 <- df.self_report %>% select(all_of(columns_of_instance0))
df.self_report.instance0.collapse <- as.data.frame(df.self_report.instance0 %>%
                                 rowwise() %>%
                                 mutate(selfreport_inst0 = paste(c_across(-f.eid), collapse = ";")) %>%
                                 ungroup() %>%
                                 dplyr::select(f.eid, selfreport_inst0))

columns_of_instance0 <- c("f.eid", names(df.self_report_age_year)[grepl("f.87.0.*", names(df.self_report_age_year))])
df.self_report_age_year.instance0 <- df.self_report_age_year %>% select(all_of(columns_of_instance0))
df.self_report_age_year.instance0.collapse <- as.data.frame(df.self_report_age_year.instance0 %>%
                                                    rowwise() %>%
                                                    mutate(selfreport_age_year_inst0 = paste(c_across(-f.eid), collapse = ";")) %>%
                                                    ungroup() %>%
                                                    dplyr::select(f.eid, selfreport_age_year_inst0))

df.self_report_code_date.instance0 <- merge(df.self_report.instance0.collapse, df.self_report_age_year.instance0.collapse, by="f.eid", all=T)
df.self_report_code_date.instance0 <- merge(df.self_report_code_date.instance0, df.year_of_birth, by="f.eid", all.x=T)

df.self_report_date.instance0.collapse <- as.data.frame(df.self_report_code_date.instance0 %>%
                                      rowwise() %>%
                                      mutate(selfreport_date_inst0 = paste(sapply(strsplit(selfreport_age_year_inst0, ";")[[1]], 
                                                                                process_element, 
                                                                                year_of_birth=f.34.0.0), collapse=";")) %>%
                                      dplyr::select(f.eid, selfreport_inst0, selfreport_age_year_inst0, selfreport_date_inst0))

## Second instance
columns_of_instance1 <- c("f.eid", names(df.self_report)[grepl("f.20002.1.*", names(df.self_report))])
df.self_report.instance1 <- df.self_report %>% select(all_of(columns_of_instance1))
df.self_report.instance1.collapse <- as.data.frame(df.self_report.instance1 %>%
                                                     rowwise() %>%
                                                     mutate(selfreport_inst1 = paste(c_across(-f.eid), collapse = ";")) %>%
                                                     ungroup() %>%
                                                     dplyr::select(f.eid, selfreport_inst1))

columns_of_instance1 <- c("f.eid", names(df.self_report_age_year)[grepl("f.87.1.*", names(df.self_report_age_year))])
df.self_report_age_year.instance1 <- df.self_report_age_year %>% select(all_of(columns_of_instance1))
df.self_report_age_year.instance1.collapse <- as.data.frame(df.self_report_age_year.instance1 %>%
                                                              rowwise() %>%
                                                              mutate(selfreport_age_year_inst1 = paste(c_across(-f.eid), collapse = ";")) %>%
                                                              ungroup() %>%
                                                              dplyr::select(f.eid, selfreport_age_year_inst1))

df.self_report_code_date.instance1 <- merge(df.self_report.instance1.collapse, df.self_report_age_year.instance1.collapse, by="f.eid", all=T)
df.self_report_code_date.instance1 <- merge(df.self_report_code_date.instance1, df.year_of_birth, by="f.eid", all.x=T)

df.self_report_date.instance1.collapse <- as.data.frame(
            df.self_report_code_date.instance1 %>%
                rowwise() %>%
                mutate(selfreport_date_inst1 = paste(sapply(strsplit(selfreport_age_year_inst1, ";")[[1]], 
                                                            process_element, 
                                                            year_of_birth=f.34.0.0), collapse=";")) %>%
                dplyr::select(f.eid, selfreport_inst1, selfreport_age_year_inst1, selfreport_date_inst1)
            )

### Third instance
columns_of_instance2 <- c("f.eid", names(df.self_report)[grepl("f.20002.2.*", names(df.self_report))])
df.self_report.instance2 <- df.self_report %>% select(all_of(columns_of_instance2))
df.self_report.instance2.collapse <- as.data.frame(df.self_report.instance2 %>%
                                                     rowwise() %>%
                                                     mutate(selfreport_inst2 = paste(c_across(-f.eid), collapse = ";")) %>%
                                                     ungroup() %>%
                                                     dplyr::select(f.eid, selfreport_inst2))

columns_of_instance2 <- c("f.eid", names(df.self_report_age_year)[grepl("f.87.2.*", names(df.self_report_age_year))])
df.self_report_age_year.instance2 <- df.self_report_age_year %>% select(all_of(columns_of_instance2))
df.self_report_age_year.instance2.collapse <- as.data.frame(df.self_report_age_year.instance2 %>%
                                                              rowwise() %>%
                                                              mutate(selfreport_age_year_inst2 = paste(c_across(-f.eid), collapse = ";")) %>%
                                                              ungroup() %>%
                                                              dplyr::select(f.eid, selfreport_age_year_inst2))

df.self_report_code_date.instance2 <- merge(df.self_report.instance2.collapse, df.self_report_age_year.instance2.collapse, by="f.eid", all=T)
df.self_report_code_date.instance2 <- merge(df.self_report_code_date.instance2, df.year_of_birth, by="f.eid", all.x=T)

df.self_report_date.instance2.collapse <- as.data.frame(
    df.self_report_code_date.instance2 %>%
        rowwise() %>%
        mutate(selfreport_date_inst2 = paste(sapply(strsplit(selfreport_age_year_inst2, ";")[[1]], 
                                                    process_element, 
                                                    year_of_birth=f.34.0.0), collapse=";")) %>%
        dplyr::select(f.eid, selfreport_inst2, selfreport_age_year_inst2, selfreport_date_inst2)
)

### Fourth instance
columns_of_instance3 <- c("f.eid", names(df.self_report)[grepl("f.20002.3.*", names(df.self_report))])
df.self_report.instance3 <- df.self_report %>% select(all_of(columns_of_instance3))
df.self_report.instance3.collapse <- as.data.frame(df.self_report.instance3 %>%
                                                       rowwise() %>%
                                                       mutate(selfreport_inst3 = paste(c_across(-f.eid), collapse = ";")) %>%
                                                       ungroup() %>%
                                                       dplyr::select(f.eid, selfreport_inst3))

columns_of_instance3 <- c("f.eid", names(df.self_report_age_year)[grepl("f.87.2.*", names(df.self_report_age_year))])
df.self_report_age_year.instance3 <- df.self_report_age_year %>% select(all_of(columns_of_instance3))
df.self_report_age_year.instance3.collapse <- as.data.frame(df.self_report_age_year.instance3 %>%
                                                                rowwise() %>%
                                                                mutate(selfreport_age_year_inst3 = paste(c_across(-f.eid), collapse = ";")) %>%
                                                                ungroup() %>%
                                                                dplyr::select(f.eid, selfreport_age_year_inst3))

df.self_report_code_date.instance3 <- merge(df.self_report.instance3.collapse, df.self_report_age_year.instance3.collapse, by="f.eid", all=T)
df.self_report_code_date.instance3 <- merge(df.self_report_code_date.instance3, df.year_of_birth, by="f.eid", all.x=T)

df.self_report_date.instance3.collapse <- as.data.frame(
    df.self_report_code_date.instance3 %>%
        rowwise() %>%
        mutate(selfreport_date_inst3 = paste(sapply(strsplit(selfreport_age_year_inst3, ";")[[1]], 
                                                    process_element, 
                                                    year_of_birth=f.34.0.0), collapse=";")) %>%
        dplyr::select(f.eid, selfreport_inst3, selfreport_age_year_inst3, selfreport_date_inst3)
)

###
df.self_report_code_date.all_instances <- merge(df.self_report_date.collapse,
                                                df.self_report_date.instance0.collapse, by="f.eid", all=T)
df.self_report_code_date.all_instances <- merge(df.self_report_code_date.all_instances, 
                                                df.self_report_date.instance1.collapse, by="f.eid", all=T)
df.self_report_code_date.all_instances <- merge(df.self_report_code_date.all_instances, 
                                                df.self_report_date.instance2.collapse, by="f.eid", all=T)
df.self_report_code_date.all_instances <- merge(df.self_report_code_date.all_instances, 
                                                df.self_report_date.instance3.collapse, by="f.eid", all=T)

saveRDS(df.self_report_code_date.all_instances, "UKBB.20002_87.Selfreport_code_and_age_year_merged.rds")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

### This is how to query a single Interview code age/date from first instance
# Interview_code <- "1065"
# df.self_report_date.quried <- as.data.frame(
#     df.self_report_code_date.all_instances %>%
#     dplyr::select(f.eid, Interview, Interview_age_year) %>%
#     rowwise %>%
#     mutate(index = ifelse(grepl(Interview_code, Interview),
#                           which(strsplit(Interview, ";")[[1]] == Interview_code), NA),
#            date = ifelse(is.na(index), NA, strsplit(Interview_age_year, ";")[[1]][index])) %>%
#     dplyr::select(f.eid, date)
#     )
### 