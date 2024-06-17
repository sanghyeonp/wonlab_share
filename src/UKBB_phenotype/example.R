source("fnc.query_date.single_code.R")

icd10 <- c("I21", "I22", "I23", "I24", "I25.2", "I46")
icd9 <- c("410", "411", "412", "414.2", "427.5")

df <- query_date_ICD10(code=icd10[1])


####
reformat_icd9 <- function(code){
    all_subclass <- !grepl("\\.", code)
    code <- gsub("\\.", "", code)
    if(all_subclass) code <- paste0(code, "\\d*") else code <- paste0(code, "$")
    code <- paste0("^", code)
    return (code)
}

reformat_icd9 <- function(code){
    all_subclass <- nchar(code) == 3
    code <- gsub("\\.", "", code)
    if(all_subclass) code <- paste0(code, "\\d*") else code <- paste0(code, "$")
    code <- paste0("^", code)
    return (code)
}

df.opcs4 <- readRDS("UKBB.41272_41282.OPCS4_code_and_date_merged.rds")

opcs4 <- c("K40", "K41", "K42", "K44")

search_pattern <- paste(sapply(opcs4, reformat_icd9), collapse = "|")
search_pattern

df.icd91 <- as.data.frame(df.icd9 %>%
                              rowwise() %>%
                              mutate(index = ifelse(grepl(search_pattern, ICD9_main), 
                                                  paste(which(grepl(search_pattern, strsplit(ICD9_main, ";")[[1]])), collapse=";"), NA), #
                                     ICD9_query = ifelse(!is.na(index), 
                                                  paste(strsplit(ICD9_main, ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
                                     date_query = ifelse(!is.na(index), 
                                                   paste(strsplit(ICD9_date, ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
                                     date_earliest = ifelse(is.na(date_query),
                                                       NA, ifelse(!grepl(";", date_query), 
                                                                  date_query,
                                                                  as.character(min(as.Date(strsplit(date_query, ";")[[1]], format="%Y-%m-%d"))))))
                        )

df.icd101 <- df.icd10

df.icd101 <- df.icd101 %>%
    rowwise() %>%
    mutate(min_date = ifelse(is.na(date),
                             NA, ifelse(!grepl(";", date), 
                                        date,
                                        as.character(min(as.Date(strsplit(date, ";")[[1]], format="%Y-%m-%d"))))))

row <- df.icd10[17563, ]
date1 <- row$date

min(as.Date(strsplit(date1, ";")[[1]], format="%Y-%m-%d"))

y <- as.Date(date1, format="%Y-%m-%d")
typeof(y)

library(lubridate)
x <- ymd(date1)
typeof(x)

index1 <- row$index

strsplit(date1, ";")[[1]][as.integer(strsplit(index1, ";")[[1]])]
