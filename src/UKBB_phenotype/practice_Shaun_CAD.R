source("fnc.query_date.R")

icd10 <- c("I21", "I22", "I23", "I24", "I25.2", "I46")
icd9 <- c("410", "411", "412", "414.2", "427.5")
opcs4 <- c("K40", "K41", "K42", "K44", "K45", "K46", "K47.1", "K49", "K50", "K75")
self_report <- "1075"

start <- Sys.time()
df.icd10 <- query_date(code_type="ICD10", code=icd10, simplify=TRUE)
print(paste0("Elapsed: ", Sys.time() - start))

start <- Sys.time()
df.icd9 <- query_date(code_type="ICD9", code=icd9, simplify=TRUE)
print(paste0("Elapsed: ", Sys.time() - start))

start <- Sys.time()
df.opcs4 <- query_date(code_type="OPCS4", code=opcs4, simplify=TRUE)
print(paste0("Elapsed: ", Sys.time() - start))

start <- Sys.time()
df.self_report <- query_date(code_type="Self-report", code=self_report, self_report_instance=0, simplify=TRUE)
print(paste0("Elapsed: ", Sys.time() - start))

#####

df.final <- query_date_multiple_source(queried_code_type=c("ICD10", "ICD9", "OPCS4", "Self-report"),
                                       df_query.ICD10 = df.icd10, df_query.ICD9 = df.icd9,
                                       df_query.OPCS4 = df.opcs4, df_query.self_report = df.self_report)

query_date_multiple_source <- function(queried_code_type, df_query.ICD10=NA, df_query.ICD9=NA, 
                                        df_query.OPCS4=NA, df_query.self_report=NA, simplify=FALSE){
    # queried_code_type: 
    n_quried <- length(queried_code_type)
    for (idx in 1:n_quried){
        if (idx == 1){
            if (queried_code_type[idx] == "ICD10"){
                df <- df_query.ICD10
            } else if (queried_code_type[idx] =="ICD9"){
                df <- df_query.ICD9
            } else if (queried_code_type[idx] == "OPCS4"){
                df <- df_query.OPCS4
            } else if (queried_code_type[idx] == "Self-report"){
                df <- df_query.self_report
            }
        } else{
            if (queried_code_type[idx] == "ICD10"){
                df <- merge(df, df_query.ICD10, by="f.eid", all=T)
            } else if (queried_code_type[idx] =="ICD9"){
                df <- merge(df, df_query.ICD9, by="f.eid", all=T)
            } else if (queried_code_type[idx] == "OPCS4"){
                df <- merge(df, df_query.OPCS4, by="f.eid", all=T)
            } else if (queried_code_type[idx] == "Self-report"){
                df <- merge(df, df_query.self_report, by="f.eid", all=T)
            }
        }
    }
    
    df <- as.data.frame(df %>%
             rowwise() %>%
             mutate(date_earliest.collapse = paste(c_across(-f.eid), collapse = ";"),
                    date_earliest = ifelse(grepl(paste0("^", paste(rep(";", n_quried-1), collapse=""), "$"), date_earliest.collapse),
                                           NA, 
                                           as.character(min(as.Date(strsplit(date_earliest.collapse, ";")[[1]], format="%Y-%m-%d"), na.rm=T)))) %>%
             dplyr::select(-date_earliest.collapse)
            )
    
    if (simplify) df <- df %>% dplyr::select(f.eid, date_earliest)
    
    return (df)
}

df.opcs41 <- NA

df2 <- df.icd10 %>%
    left_join(df.icd9, by="f.eid") %>%
    left_join(df.opcs41, by="f.eid") %>%
    left_join(df.self_report, by="f.eid") %>%
    mutate_all(~replace_na(., ""))

col_name <- colnames(df)[2:length(colnames(df))]

df1 <- as.data.frame(df %>%
    rowwise() %>%
    mutate(date_earliest.collapse = paste(c_across(-f.eid), collapse = ";"),
           date_earliest = ifelse(grepl(paste0("^", paste(rep(";", length(col_name)-1), collapse=""), "$"), date_earliest.collapse),
                                  NA, 
                                  as.character(min(as.Date(strsplit(date_earliest.collapse, ";")[[1]], format="%Y-%m-%d"), na.rm=T)))) %>%
    dplyr::select(-date_earliest.collapse)
)
    
#####
row <- df[c(3, 5,8), ]
row

row <- as.data.frame(row %>%
    rowwise() %>%
    mutate(date_earliest.collapse = paste(c_across(-f.eid), collapse = ";")))

row1 <- row %>%
    rowwise() %>%
    mutate(date_earliest = ifelse(grepl(paste0("^", paste(rep(";", length(col_name)-1), collapse=""), "$"), date_earliest.collapse),
                          NA, 
                          as.character(min(as.Date(strsplit(date_earliest.collapse, ";")[[1]], format="%Y-%m-%d"), na.rm=T))))

min(as.Date(strsplit(row[2,]$date_earliest.collapse, ";")[[1]], format="%Y-%m-%d"), na.rm=T)

paste0("^", paste(rep(";", length(col_name)), collapse=""), "$")
row[1,]$date_earliest.collapse
