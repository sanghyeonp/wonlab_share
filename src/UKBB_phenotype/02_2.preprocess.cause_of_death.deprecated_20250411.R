library(dplyr)
library(tidyr)

df.assessment_date <- readRDS("UKBB.53.assessment_center_date.rds")


#####
df.cause_of_death <- readRDS("UKBB.40001.cause_of_death.rds")


######################################
### Preprocess cause of death and date ###
######################################
print("Preprocess cause of death and date")
start_time <- Sys.time()
df.cause_of_death <- readRDS("UKBB.40001.cause_of_death.rds")
df.death_date <- readRDS("UKBB.40000.date_of_death.rds")

df.cause_of_death.instance_collapse <- as.data.frame(df.cause_of_death %>%
        rowwise() %>%
        mutate(Cause_of_death = paste(c_across(-f.eid), collapse = ";")) %>%
        ungroup() %>%
        dplyr::select(f.eid, Cause_of_death))

idate_columns <- sapply(df.death_date, inherits, "IDate")
df.death_date.instance_collapse <- as.data.frame(df.death_date %>%
    mutate(across(all_of(names(which(idate_columns))), as.character)) %>%
    mutate_all(~replace_na(., "")) %>%
    rowwise() %>%
    mutate(Cause_of_death_date = paste(c_across(-f.eid), collapse = ";")) %>%
    ungroup() %>%
    dplyr::select(f.eid, Cause_of_death_date))

df.cause_of_death_and_date <- merge(df.cause_of_death.instance_collapse, 
                                    df.death_date.instance_collapse, by="f.eid", all=T)

saveRDS(df.cause_of_death_and_date, "UKBB.40000_40001.Cause_of_death_and_date_merged.rds")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

### This is how to query a single ICD-10 code date
ICD10_code <- "C844"
df.cause_of_death.quried <- as.data.frame(
    df.cause_of_death_and_date %>%
    rowwise %>%
    mutate(index = ifelse(grepl(ICD10_code, Cause_of_death),
                          which(strsplit(Cause_of_death, ";")[[1]] == ICD10_code), NA),
           date = ifelse(is.na(index), NA, strsplit(Cause_of_death_date, ";")[[1]][index])) %>%
    dplyr::select(f.eid, date)
    )

x <- df.cause_of_death_and_date %>%
dplyr::rowwise() %>%
    dplyr::mutate(
        # Index of ICD9 code of interest
        index = ifelse(grepl(search_pattern, Cause_of_death), 
                       paste(which(grepl(exact_pattern, strsplit(Cause_of_death, ";")[[1]])), collapse=";"), NA), 
        # Query ICD9 code of interest using the index (double check if correct code is queried)
        code_query = ifelse(!is.na(index), 
                            paste(strsplit(Cause_of_death, ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
        # Query date of diagnosis of interest using the index
        date_query = ifelse(!is.na(index), 
                            paste(strsplit(Cause_of_death_date, ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
        # Obtain the earliest date of diagnosis
        date_earliest = ifelse(is.na(date_query),
                               NA, ifelse(!grepl(";", date_query), 
                                          date_query,
                                          as.character(min(as.Date(strsplit(date_query, ";")[[1]], format="%Y-%m-%d"), na.rm=T))))) %>%
    # Drop collapsed data columns (dplyr::select(-ICD9_main, -ICD9_date))
    dplyr::select(f.eid, index, code_query, date_query, date_earliest)

### 

