### Updates
# 2025.04.12:
#    - '1065;1223;;;;;;;;;;;1143;;;;;;;;;;;' 에서 ^1223을 찾을 수는 없음.

######
library(dplyr)

query_date <- function(code_type, code, instance=0, simplify=TRUE){
    # code_type: either "ICD10", "ICD9", "OPCS4", "Self-report" or "Cause of death"
    # code: single or a vector character of ICD-10, ICD-9, or OPCS4 or Self-report codes
    # instance: either all, 0, 1, 2, 3
    # simplify: if TRUE, return f.eid and only the date of the earliest diagnosis,

    # Elapsed time: average 1 minutes per code type
    if (code_type == "ICD10"){
        df.query <- query_date.ICD10(code=code, simplify=simplify)
    } else if (code_type == "ICD9"){
        df.query <- query_date.ICD9(code=code, simplify=simplify)
    } else if (code_type == "OPCS4"){
        df.query <- query_date.OPCS4(code=code, simplify=simplify)
    } else if (code_type == "Self-report") {
        df.query <- query_date.self_report_verbal(code=code, instance=instance, simplify=simplify)
    } else if (code_type == "Cause of death"){
        df.query <- query_date.cause_of_death(code=code, simplify=simplify)
    } else if (code_type == "Secondary cause of death"){
        df.query <- query_date.secondary_cause_of_death(code=code, simplify=simplify)
    } else if (code_type == "Surgical operation"){
        df.query <- query_date.surgical_operation(code=code, instance=instance, simplify=simplify)
    } else if (code_type == "Touchscreen heart problem"){
        df.query <- query_date.heart_problem(code=code, instance=instance, simplify=simplify)
    }
    else {
        stop("Invalid code_type. Please choose from 'ICD10', 'ICD9', 'OPCS4', 'Self-report', 'Cause of death', 'Secondary cause of death', 'Surgical operation' or 'Touchscreen heart problem'.")
    }
    return (df.query)
}


#########################################################
# reformat_into_search_pattern <- function(code, exact_match, self_report=FALSE){
#     # 
#     # Example
#     # Looking for all sub-classes of ICD10 L71: "L71*"
#     # Looking for specific sub-class (ICD10 L48.6): either "L48.6" or "L486"
#     if (exact_match){
#         if (!self_report){
#             # Check if looking for all subclassifications
#             # all_subclass <- grepl("*", code)
#             # Remove "." from the code
#             code <- gsub("\\.", "", code)
#             # Add "*" if looking for all subclassifications, otherwise add "$"
#             if(!all_subclass) code <- paste0(code, "$")
#             # Start with the given code itself
#             code <- paste0("^", code)
#         } else{
#             code <- paste0("^", code, "$")
#         }
#     } else{ # For general search among collapsed codes
#         code <- gsub("\\.", "", code)
#     }
#     return (code)
# }

####################################################
# Query the earilest date of diagnosis from ICD-10 #
####################################################
query_date.ICD10 <- function(code, simplify=TRUE){
    # code: single or a vector character of ICD-10 codes
    #   - Example) c("I21", "I24", "I25.2", "I46")
    # simplify: if TRUE, return f.eid and only the date of the earliest diagnosis, 
    #            otherwise return f.eid, index, quried ICD-10 code, quried date of diagnosis, and the earliest date of diagnosis

    # Read the data
    df.data <- readRDS("/data1/sanghyeon/wonlab_contribute/combined/src/UKBB_phenotype/UKBB.41270_41280.ICD10_code_and_date_merged.rds")

    # Make the search pattern of given ICD-10 codes
    code <- as.character(code)
    search_pattern <- paste(sapply(code, function(x){ifelse(grepl("\\.", x), paste0("(^|;)", gsub("\\.", "", x), "(;|$)"),
                                                        paste0("(^|;)", gsub("\\.", "", x), "[^;]*(;|$)"))}),
                        collapse="|")
    # search_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=FALSE, self_report=FALSE), collapse = "|")
    # exact_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=TRUE, self_report=FALSE), collapse = "|")
    cat(paste0("Specified ICD10: ", paste(code, collapse=", "), "\n", 
            "Search pattern: ", search_pattern, "\n"))

    # Query the ICD10 codes
    df.query <- as.data.frame(df.data %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            # Index of ICD10 code of interest
            index = ifelse(grepl(search_pattern, ICD10_main),
                            paste(which(grepl(search_pattern, strsplit(ICD10_main, ";")[[1]])), collapse=";"), NA), 
            # Query ICD10 code of interest using the index (double check if correct code is queried)
            code_query = ifelse(!is.na(index), 
                        paste(strsplit(ICD10_main, ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
            # Query date of diagnosis of interest using the index
            date_query = ifelse(!is.na(index), 
                        paste(strsplit(ICD10_date, ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
            # Obtain the earliest date of diagnosis
            date_earliest = ifelse(is.na(date_query),
                            NA, ifelse(!grepl(";", date_query), 
                                        date_query,
                                        as.character(min(as.Date(strsplit(date_query, ";")[[1]], format="%Y-%m-%d")))))) %>%
        # Drop collapsed data columns (dplyr::select(-ICD10_main, -ICD10_date))
        dplyr::select(f.eid, index, code_query, date_query, date_earliest) %>%
        dplyr::mutate(f.eid = as.character(f.eid), index = as.character(index), 
                    code_query = as.character(code_query), date_query = as.character(date_query), 
                    date_earliest = as.character(date_earliest)) %>%
        dplyr::rename(index.ICD10=index, code_query.ICD10=code_query, 
                    date_query.ICD10=date_query, date_earliest.ICD10=date_earliest)
    )

    if (simplify) df.query <- df.query %>% dplyr::select(f.eid, date_earliest.ICD10)
    
    return(df.query)
}

####################################################
# Query the earilest date of diagnosis from ICD-9 #
####################################################
query_date.ICD9 <- function(code, simplify=TRUE){
    # code: single or a vector character of ICD-9 codes
    #   - Example) c("410", "411", "414.2", "427.5")
    # simplify: if TRUE, return f.eid and only the date of the earliest diagnosis, 
    #            otherwise return f.eid, index, quried ICD-9 code, quried date of diagnosis, and the earliest date of diagnosis

    # Read the data
    df.data <- readRDS("/data1/sanghyeon/wonlab_contribute/combined/src/UKBB_phenotype/UKBB.41271_41281.ICD9_code_and_date_merged.rds")

    # Make the search pattern of given ICD-9 codes
    code <- as.character(code)
    search_pattern <- paste(sapply(code, function(x){ifelse(grepl("\\.", x), paste0("(^|;)", gsub("\\.", "", x), "(;|$)"),
                                                        paste0("(^|;)", gsub("\\.", "", x), "[^;]*(;|$)"))}),
                        collapse="|")
    # search_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=FALSE, self_report=FALSE), collapse = "|")
    # exact_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=TRUE, self_report=FALSE), collapse = "|")
    cat(paste0("Specified ICD9: ", paste(code, collapse=", "), "\n", 
            "Search pattern: ", search_pattern, "\n"))

    # Query the ICD9 codes
    df.query <- as.data.frame(df.data %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            # Index of ICD9_main code of interest
            index = ifelse(grepl(search_pattern, ICD9_main),
                        paste(which(grepl(search_pattern, strsplit(ICD9_main, ";")[[1]])), collapse=";"), NA), 
            # Query ICD9 code of interest using the index (double check if correct code is queried)
            code_query = ifelse(!is.na(index), 
                        paste(strsplit(ICD9_main, ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
            # Query date of diagnosis of interest using the index
            date_query = ifelse(!is.na(index), 
                        paste(strsplit(ICD9_date, ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
            # Obtain the earliest date of diagnosis
            date_earliest = ifelse(is.na(date_query),
                            NA, ifelse(!grepl(";", date_query), 
                                        date_query,
                                        as.character(min(as.Date(strsplit(date_query, ";")[[1]], format="%Y-%m-%d"), na.rm=T))))) %>%
        # Drop collapsed data columns (dplyr::select(-ICD9_main, -ICD9_date))
        dplyr::select(f.eid, index, code_query, date_query, date_earliest) %>%
        dplyr::mutate(f.eid = as.character(f.eid), index = as.character(index), 
            code_query = as.character(code_query), date_query = as.character(date_query), 
            date_earliest = as.character(date_earliest)) %>%
        dplyr::rename(index.ICD9=index, code_query.ICD9=code_query, 
                    date_query.ICD9=date_query, date_earliest.ICD9=date_earliest)
    )

    if (simplify) df.query <- df.query %>% dplyr::select(f.eid, date_earliest.ICD9)
    
    return(df.query)
}

####################################################
# Query the earilest date of diagnosis from OPCS4 #
####################################################
query_date.OPCS4 <- function(code, simplify=TRUE){
    # code: single or a vector character of OPCS4 codes
    #   - Example) c("K40", "K41", "K42", "K44.1", "K44.2")
    # simplify: if TRUE, return f.eid and only the date of the earliest diagnosis, 
    #            otherwise return f.eid, index, quried OPCS4 code, quried date of diagnosis, and the earliest date of diagnosis

    # Read the data
    df.data <- readRDS("/data1/sanghyeon/wonlab_contribute/combined/src/UKBB_phenotype/UKBB.41272_41282.OPCS4_code_and_date_merged.rds")

    # Make the search pattern of given OPCS4 codes
    code <- as.character(code)
    search_pattern <- paste(sapply(code, function(x){ifelse(grepl("\\.", x), paste0("(^|;)", gsub("\\.", "", x), "(;|$)"),
                                                        paste0("(^|;)", gsub("\\.", "", x), "[^;]*(;|$)"))}),
                        collapse="|")
    # search_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=FALSE, self_report=FALSE), collapse = "|")
    # exact_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=TRUE, self_report=FALSE), collapse = "|")
    cat(paste0("Specified OPCS4: ", paste(code, collapse=", "), "\n", 
            "Search pattern: ", search_pattern, "\n"))

    # Query the OPCS4 codes
    df.query <- as.data.frame(df.data %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            # Index of OPCS4 code of interest
            index = ifelse(grepl(search_pattern, OPCS4),
                        paste(which(grepl(search_pattern, strsplit(OPCS4, ";")[[1]])), collapse=";"), NA), 
            # Query OPCS4 code of interest using the index (double check if correct code is queried)
            code_query = ifelse(!is.na(index), 
                        paste(strsplit(OPCS4, ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
            # Query date of diagnosis of interest using the index
            date_query = ifelse(!is.na(index), 
                        paste(strsplit(OPCS4_date, ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
            # Obtain the earliest date of diagnosis
            date_earliest = ifelse(is.na(date_query),
                            NA, ifelse(!grepl(";", date_query), 
                                        date_query,
                                        as.character(min(as.Date(strsplit(date_query, ";")[[1]], format="%Y-%m-%d"), na.rm=T))))) %>%
        # Drop collapsed data columns (dplyr::select(-OPCS4, -OPCS4_date))
        dplyr::select(f.eid, index, code_query, date_query, date_earliest) %>%
        dplyr::mutate(f.eid = as.character(f.eid), index = as.character(index), 
                    code_query = as.character(code_query), date_query = as.character(date_query), 
                    date_earliest = as.character(date_earliest)) %>%
        dplyr::rename(index.OPCS4=index, code_query.OPCS4=code_query, 
                    date_query.OPCS4=date_query, date_earliest.OPCS4=date_earliest)
    )

    if (simplify) df.query <- df.query %>% dplyr::select(f.eid, date_earliest.OPCS4)
    
    return(df.query)
}

#########################################################
# Query the earilest date of diagnosis from self-report #
#########################################################
query_date.self_report_verbal <- function(code, instance=0, simplify=TRUE){
    # code: single or a vector character of self-report codes
    #   - Example) 1065
    # simplify: if TRUE, return f.eid and only the date of the earliest diagnosis, 
    #            otherwise return f.eid, index, quried self-report code, quried date of diagnosis, and the earliest date of diagnosis

    # Note: date with 3000-12-12 is considered as having disease but the exact date is unknown (according to coding 37 from field 87) 
    #    -3 represents "Preferred not to answer"
    #    -1 represents "Time uncertain/unknown"

    # Read the data
    df.data <- readRDS("/data1/sanghyeon/wonlab_contribute/combined/src/UKBB_phenotype/UKBB.20002_87.Selfreport_code_and_age_year_merged.rds")

    # Make the search pattern of given self-report codes
    code <- as.character(code)
    code <- gsub("\\.", "", code)
    search_pattern <- paste(paste0("(^|;)", code, "(;|$)"), collapse="|")
    # search_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=FALSE, self_report=FALSE), collapse = "|")
    # exact_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=TRUE, self_report=FALSE), collapse = "|")
    cat(paste0("Specified self-report codes: ", paste(code, collapse=", "), "\n", 
            "Search pattern: ", search_pattern, "\n"))
    
    if (instance == "all"){
        instance_to_retain <- c("f.eid", names(df.data)[grepl("all", names(df.data))])
        df.data <- df.data %>% select(all_of(instance_to_retain))
        suffix <- "selfreport_all"
    } else{
        instance_to_retain <- c("f.eid", names(df.data)[grepl(paste0("inst", instance), names(df.data))])
        df.data <- df.data %>% select(all_of(instance_to_retain))
        suffix <- paste0("selfreport_inst", instance)
    }

    code_col <- colnames(df.data)[2]; date_col <- colnames(df.data)[4]

    # Query the Self-report codes
    df.query <- as.data.frame(df.data %>%
                dplyr::rowwise() %>%
                dplyr::mutate(
                    # Index of self-report code of interest
                    index = ifelse(grepl(search_pattern, !!as.name(code_col)), 
                                    paste(which(grepl(search_pattern, strsplit(!!as.name(code_col), ";")[[1]])), collapse=";"), NA), 
                    # Query self-report code of interest using the index (double check if correct code is queried)
                    code_query = ifelse(!is.na(index), 
                                        paste(strsplit(!!as.name(code_col), ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
                    # Query date of diagnosis of interest using the index
                    date_query = ifelse(!is.na(index), 
                                        paste(strsplit(!!as.name(date_col), ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
                    # Obtain the earliest date of diagnosis
                    date_earliest = ifelse(is.na(date_query),
                                            NA, ifelse(!grepl(";", date_query), 
                                                        date_query,
                                                        as.character(min(as.Date(strsplit(date_query, ";")[[1]], format="%Y-%m-%d"), na.rm=T))))) %>%
                # Drop collapsed data columns (dplyr::select(-!!as.name(code_col), -!!as.name(date_col)))
                dplyr::select(f.eid, index, code_query, date_query, date_earliest) %>%
                dplyr::mutate(f.eid = as.character(f.eid), index = as.character(index), 
                    code_query = as.character(code_query), date_query = as.character(date_query), 
                    date_earliest = as.character(date_earliest)) %>%
                dplyr::rename(!!as.name(paste0("index.", suffix)) := index,
                            !!as.name(paste0("code_query.", suffix)) := code_query,
                            !!as.name(paste0("date_query.", suffix)) := date_query,
                            !!as.name(paste0("date_earliest.", suffix)) := date_earliest)
                )

    

    if (simplify) df.query <- df.query %>% dplyr::select(f.eid, !!as.name(paste0("date_earliest.", suffix)))

    cat("\n::Please note to check for missing date coded as 3000-12-31\n")

    return (df.query)
}


###########################################
# Query the death date for cause of death #
###########################################
query_date.cause_of_death <- function(code, simplify=TRUE){
    # code: single or a vector character of ICD-10 codes
    #   - Example) c("I21", "I24", "I25.2", "I46")
    # simplify: if TRUE, return f.eid and only the date of the earliest diagnosis, 
    #            otherwise return f.eid, index, quried ICD-10 code for death cause and quried death date

    # Read the data
    df.data <- readRDS("/data1/sanghyeon/wonlab_contribute/combined/src/UKBB_phenotype/UKBB.40000_40001.Cause_of_death_and_date_merged.rds")

    # Make the search pattern of given ICD-10 codes
    code <- as.character(code)
    search_pattern <- paste(sapply(code, function(x){ifelse(grepl("\\.", x), paste0("(^|;)", gsub("\\.", "", x), "(;|$)"),
                                                        paste0("(^|;)", gsub("\\.", "", x), "[^;]*(;|$)"))}),
                        collapse="|")
    # search_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=FALSE, self_report=FALSE), collapse = "|")
    # exact_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=TRUE, self_report=FALSE), collapse = "|")
    cat(paste0("Specified ICD10: ", paste(code, collapse=", "), "\n", 
            "Search pattern: ", search_pattern, "\n"))

    # Query the ICD10 codes
    df.query <- as.data.frame(df.data %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            # Index of ICD10 code of interest
            index = ifelse(grepl(search_pattern, Cause_of_death),
                            paste(which(grepl(search_pattern, strsplit(Cause_of_death, ";")[[1]])), collapse=";"), NA), 
            # Query ICD10 code of interest using the index (double check if correct code is queried)
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
        dplyr::select(f.eid, index, code_query, date_query, date_earliest) %>%
        dplyr::mutate(f.eid = as.character(f.eid), index = as.character(index), 
            code_query = as.character(code_query), date_query = as.character(date_query), 
            date_earliest = as.character(date_earliest)) %>%
        dplyr::rename(index.cause_of_death=index, code_query.cause_of_death=code_query, 
                    date_query.cause_of_death=date_query, date_earliest.cause_of_death=date_earliest)
    )

    if (simplify) df.query <- df.query %>% dplyr::select(f.eid, date_earliest.cause_of_death)
    
    return(df.query)
}

###########################################
# Query the death date for secondary cause of death #
###########################################
query_date.secondary_cause_of_death <- function(code, simplify=TRUE){
    # code: single or a vector character of ICD-10 codes
    #   - Example) c("I21", "I24", "I25.2", "I46")
    # simplify: if TRUE, return f.eid and only the date of the earliest diagnosis, 
    #            otherwise return f.eid, index, quried ICD-10 code for death cause and quried death date

    # Read the data
    df.data <- readRDS("/data1/sanghyeon/wonlab_contribute/combined/src/UKBB_phenotype/UKBB.40002_40001.Secondary_cause_of_death_and_date_merged.rds")

    # Make the search pattern of given ICD-10 codes
    code <- as.character(code)
    search_pattern <- paste(sapply(code, function(x){ifelse(grepl("\\.", x), paste0("(^|;|,)", gsub("\\.", "", x), "(,|;|$)"),
                                                            paste0("(^|;|,)", gsub("\\.", "", x), "[^;]*(,|;|$)"))}),
                            collapse="|")
    # search_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=FALSE, self_report=FALSE), collapse = "|")
    # exact_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=TRUE, self_report=FALSE), collapse = "|")
    cat(paste0("Specified ICD10: ", paste(code, collapse=", "), "\n", 
            "Search pattern: ", search_pattern, "\n"))

    # Query the ICD10 codes in first instance
    df.query <- df.data
    df.query$has_code <- sapply(df.query$Secondary_cause_of_death, function(x){grepl(search_pattern, x)})
    df.query$index_inst <- apply(df.query, 1, function(x){ifelse(x[["has_code"]], 
                                                        paste(which(grepl(search_pattern, strsplit(x[["Secondary_cause_of_death"]], ";")[[1]])), collapse=";"), 
                                                        NA)})
    df.query$index_code <- apply(df.query, 1, function(x){ifelse(x[["has_code"]], 
                                                            paste(which(grepl(search_pattern, strsplit(strsplit(x[["Secondary_cause_of_death"]], ";")[[1]][as.integer(x[["index_inst"]])], ",")[[1]])), collapse=","), 
                                                            NA)})
    df.query$code_query <- apply(df.query, 1 ,function(x){ifelse(x[["has_code"]],
                                                                paste(strsplit(x[["Secondary_cause_of_death"]], ";")[[1]][as.integer(strsplit(x[["index_inst"]], ";")[[1]])], collapse=";"),
                                                                NA)})
    df.query$date_query <- apply(df.query, 1 ,function(x){ifelse(x[["has_code"]],
                                                                paste(strsplit(x[["Secondary_cause_of_death_date"]], ";")[[1]][as.integer(strsplit(x[["index_inst"]], ";")[[1]])], collapse=";"),
                                                                NA)})
    df.query$date_earliest <- apply(df.query, 1, function(x){ifelse(x[["has_code"]], 
                                                                            as.character(min(as.Date(strsplit(x[["date_query"]], ";")[[1]], format="%Y-%m-%d"), na.rm=T)),
                                                                            NA)})
    df.query <- df.query %>%
        dplyr::select(f.eid, index_inst, index_code, code_query, date_query, date_earliest) %>%
        dplyr::mutate(f.eid = as.character(f.eid), index_inst = as.character(index_inst), index_code = as.character(index_code),
                    code_query = as.character(code_query), date_query = as.character(date_query), 
                    date_earliest = as.character(date_earliest)) %>%
        dplyr::rename(index_inst.secondary_cause_of_death=index_inst, index_code.secondary_cause_of_death=index_code,
                    code_query.secondary_cause_of_death=code_query, date_query.secondary_cause_of_death=date_query, date_earliest.secondary_cause_of_death=date_earliest)

    # Query the ICD10 codes in thrid instance
    if (simplify) df.query <- df.query %>% dplyr::select(f.eid, date_earliest.secondary_cause_of_death)
    
    return(df.query)
}


#########################################################
# Query the earilest date of surgical operation #
#########################################################
query_date.surgical_operation <- function(code, instance=0, simplify=TRUE){
    # code: single or a vector character of self-report codes
    #   - Example) 1065
    # simplify: if TRUE, return f.eid and only the date of the earliest diagnosis, 
    #            otherwise return f.eid, index, quried self-report code, quried date of diagnosis, and the earliest date of diagnosis

    # Note: date with 3000-12-12 is considered as having disease but the exact date is unknown (according to coding 37 from field 87) 
    #    -3 represents "Preferred not to answer"
    #    -1 represents "Time uncertain/unknown"

    # Read the data
    df.data <- readRDS("/data1/sanghyeon/wonlab_contribute/combined/src/UKBB_phenotype/UKBB.20004_92.Surgical_operation_and_age_year_merged.rds")

    # Make the search pattern of given self-report codes
    code <- as.character(code)
    code <- gsub("\\.", "", code)
    search_pattern <- paste(paste0("(^|;)", code, "(;|$)"), collapse="|")
    # search_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=FALSE, self_report=FALSE), collapse = "|")
    # exact_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=TRUE, self_report=FALSE), collapse = "|")
    cat(paste0("Specified surgical operation codes: ", paste(code, collapse=", "), "\n", 
            "Search pattern: ", search_pattern, "\n"))
    
    instance_to_retain <- c("f.eid", names(df.data)[grepl(paste0("inst", instance), names(df.data))])
    df.data <- df.data %>% select(all_of(instance_to_retain))
    suffix <- paste0("surgical_op_inst", instance)
    
    code_col <- colnames(df.data)[2]; date_col <- colnames(df.data)[3]

    # Query the Self-report codes
    df.query <- as.data.frame(df.data %>%
                dplyr::rowwise() %>%
                dplyr::mutate(
                    # Index of self-report code of interest
                    index = ifelse(sum(grepl(search_pattern, !!as.name(code_col))) > 0, 
                                    paste(which(grepl(search_pattern, strsplit(!!as.name(code_col), ";")[[1]])), collapse=";"), NA), 
                    # Query self-report code of interest using the index (double check if correct code is queried)
                    code_query = ifelse(!is.na(index), 
                                        paste(strsplit(!!as.name(code_col), ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
                    # Query date of diagnosis of interest using the index
                    date_query = ifelse(!is.na(index), 
                                        paste(strsplit(!!as.name(date_col), ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
                    # Obtain the earliest date of diagnosis
                    date_earliest = ifelse(is.na(date_query),
                                            NA, ifelse(!grepl(";", date_query), 
                                                        date_query,
                                                        as.character(min(as.Date(strsplit(date_query, ";")[[1]], format="%Y-%m-%d"), na.rm=T))))) %>%
                # Drop collapsed data columns (dplyr::select(-!!as.name(code_col), -!!as.name(date_col)))
                dplyr::select(f.eid, index, code_query, date_query, date_earliest) %>%
                dplyr::mutate(f.eid = as.character(f.eid), index = as.character(index), 
                    code_query = as.character(code_query), date_query = as.character(date_query), 
                    date_earliest = as.character(date_earliest)) %>%
                dplyr::rename(!!as.name(paste0("index.", suffix)) := index,
                            !!as.name(paste0("code_query.", suffix)) := code_query,
                            !!as.name(paste0("date_query.", suffix)) := date_query,
                            !!as.name(paste0("date_earliest.", suffix)) := date_earliest)
                )

    if (simplify) df.query <- df.query %>% dplyr::select(f.eid, !!as.name(paste0("date_earliest.", suffix)))

    cat("\n::Please note to check for missing date coded as 3000-12-31\n")

    return (df.query)
}

query_date.heart_problem <- function(code, instance=0, simplify=TRUE){
    # Read the data
    df.data <- readRDS("/data1/sanghyeon/wonlab_contribute/combined/src/UKBB_phenotype/UKBB.6150_3894_3627.heart_problem_date_merged.rds")

    # Make the search pattern of given self-report codes
    code <- as.character(code)
    code <- gsub("\\.", "", code)
    search_pattern <- paste(paste0("(^|;)", code, "(;|$)"), collapse="|")
    # search_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=FALSE, self_report=FALSE), collapse = "|")
    # exact_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=TRUE, self_report=FALSE), collapse = "|")
    cat(paste0("Specified heart problem codes: ", paste(code, collapse=", "), "\n", 
            "Search pattern: ", search_pattern, "\n"))

    instance_to_retain <- c("f.eid", names(df.data)[grepl(paste0("inst", instance), names(df.data))])
    df.data <- df.data %>% select(all_of(instance_to_retain))
    if (code == "1"){
        col_to_retain <- c("f.eid", paste0("touchscreen.inst", instance), paste0("date.heart_attack.inst", instance))
        df.data <- df.data %>% select(all_of(col_to_retain))
    } else if (code == "2"){
        col_to_retain <- c("f.eid", paste0("touchscreen.inst", instance), paste0("date.angina.inst", instance))
        df.data <- df.data %>% select(all_of(col_to_retain))
    }

    suffix <- paste0("heart_problem_code", code, "_inst", instance)
    code_col <- colnames(df.data)[2]; date_col <- colnames(df.data)[3]

    # Query the Self-report codes
    df.query <- as.data.frame(df.data %>%
                dplyr::rowwise() %>%
                dplyr::mutate(
                    # Index of self-report code of interest
                    index = ifelse(grepl(search_pattern, !!as.name(code_col)), 
                                    paste(which(grepl(search_pattern, strsplit(!!as.name(code_col), ";")[[1]])), collapse=";"), NA), 
                    # Query self-report code of interest using the index (double check if correct code is queried)
                    code_query = ifelse(!is.na(index), 
                                        paste(strsplit(!!as.name(code_col), ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
                    # Query date of diagnosis of interest using the index
                    date_query = ifelse(!is.na(index), 
                                        paste(strsplit(!!as.name(date_col), ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
                    # Obtain the earliest date of diagnosis
                    date_earliest = ifelse(is.na(date_query),
                                            NA, ifelse(!grepl(";", date_query), 
                                                        date_query,
                                                        as.character(min(as.Date(strsplit(date_query, ";")[[1]], format="%Y-%m-%d"), na.rm=T))))) %>%
                # Drop collapsed data columns (dplyr::select(-!!as.name(code_col), -!!as.name(date_col)))
                dplyr::select(f.eid, index, code_query, date_query, date_earliest) %>%
                dplyr::mutate(f.eid = as.character(f.eid), index = as.character(index), 
                    code_query = as.character(code_query), date_query = as.character(date_query), 
                    date_earliest = as.character(date_earliest)) %>%
                dplyr::rename(!!as.name(paste0("index.", suffix)) := index,
                            !!as.name(paste0("code_query.", suffix)) := code_query,
                            !!as.name(paste0("date_query.", suffix)) := date_query,
                            !!as.name(paste0("date_earliest.", suffix)) := date_earliest)
                )

    if (simplify) df.query <- df.query %>% dplyr::select(f.eid, !!as.name(paste0("date_earliest.", suffix)))

    cat("\n::Please note to check for missing date coded as 3000-12-31\n")

    return (df.query)
}


#########################################################
parse_pheno <- function(code_list, pheno_col, ...){
    # code_list: a list of code names queried separated by semicolon. 
    #       For example, if ICD10, ICD9, Self-report were used, then x(code_list="ICD10;ICD9;Self-report")
    #       Possible options: "ICD10", "ICD9", "Self-report", "Cause of death", 
    # pheno_col: Name of the phenotype of interest
    #       Output dataframe will have <pheno_col> and <pheno_col>.Date
    # dataframe_list: Specify dataframe object in the same order as code listed in code_list.
    #       For example (continued from above example), x(code_list="ICD10;ICD9;Self-report", df.icd10, df.icd9, df.self)
    # Overall example: phenotype of interest is T2D and queried ICD10 and Self-report where the queried dataframe is held by df.icd10 and df.sr, respectively
    #       parse_pheno("ICD10;ICD9", "T2D", df.icd10, df.sr)
    
    df_list <- list(...); code_list <- strsplit(code_list, ";")[[1]]

    code_rename <- c("ICD10"="ICD10", "ICD9"="ICD9",
                     "Self-report"="selfreport", 
                     "Cause of death"="cause_of_death",
                     "Secondary cause of death"="secondary_cause_of_death")
    
    df.merged <- data.frame(); selfreport_col_suffix <- "NA"
    for (i in 1:length(code_list)){
        df.tmp <- df_list[[i]]; code <- as.vector(code_rename[code_list[i]])
        df.tmp <- df.tmp %>%
            dplyr::select(matches("^f\\.eid$|^code_query\\.|^date_earliest\\."))
        if (code == "ICD10"){
            df.tmp <- df.tmp %>% mutate(pheno.ICD10=ifelse(is.na(code_query.ICD10) & is.na(date_earliest.ICD10), 0, 
                                                           ifelse(!is.na(code_query.ICD10) & is.na(date_earliest.ICD10), NA, 1)))
        } else if (code == "ICD9"){
            df.tmp <- df.tmp %>% mutate(pheno.ICD9=ifelse(is.na(code_query.ICD9) & is.na(date_earliest.ICD9), 0, 
                                                          ifelse(!is.na(code_query.ICD9) & is.na(date_earliest.ICD9), NA, 1)))
        } else if (code == "selfreport"){
            selfreport_col_suffix <- strsplit(colnames(df.tmp)[2], "\\.")[[1]][2]
            col_code <- colnames(df.tmp)[2]; col_date <- colnames(df.tmp)[3]
            df.tmp <- df.tmp %>% mutate(!!as.name(paste0("pheno.", selfreport_col_suffix)):=ifelse(is.na(!!as.name(col_code)) & is.na(!!as.name(col_date)), 0, 
                                                                                                   ifelse(!is.na(!!as.name(col_code)) & (is.na(!!as.name(col_date)) | !!as.name(col_date) == "3000-12-31"), NA, 1)))
        } else if (code == "cause_of_death"){
            df.tmp <- df.tmp %>% mutate(pheno.cause_of_death=ifelse(is.na(code_query.cause_of_death) & is.na(date_earliest.cause_of_death), 0, 
                                                                    ifelse(!is.na(code_query.cause_of_death) & is.na(date_earliest.cause_of_death), NA, 1)))
        } else if (code == "secondary_cause_of_death"){
            df.tmp <- df.tmp %>% mutate(pheno.secondary_cause_of_death=ifelse(is.na(code_query.secondary_cause_of_death) & is.na(date_earliest.secondary_cause_of_death), 0, 
                                                                    ifelse(!is.na(code_query.secondary_cause_of_death) & is.na(date_earliest.secondary_cause_of_death), NA, 1)))
        } else{
            stop("Unknown code type")
        }
        if (i == 1){
            df.merged <- df.tmp
        } else{
            df.merged <- df.merged %>% left_join(df.tmp, by="f.eid")
        }
    }
    
    df.pheno <- df.merged %>%
        dplyr::select(f.eid, matches("^pheno\\.")) %>%
        mutate(pheno.final = ifelse(rowSums(select(., -1), na.rm = TRUE) > 0, 1, 0)) %>%
        dplyr::select(f.eid, pheno.final)
    
    df.merged <- df.merged %>%
        left_join(df.pheno, by="f.eid")
    
    if ("Self-report" %in% code_list){
        df.merged <- df.merged %>%
            rename(pheno.tmp = pheno.final) %>%
            mutate(pheno.final = ifelse(pheno.tmp == 0 & is.na(!!as.name(paste0("pheno.", selfreport_col_suffix))), NA, pheno.tmp)) %>%
            dplyr::select(-pheno.tmp)
    }
    
    df.pheno_date <- df.merged %>% 
        dplyr::select(f.eid, pheno.final, matches("^date_earliest\\.")) %>%
        mutate(across(3:last_col(), ~as.Date(.x, format = "%Y-%m-%d"))) %>%
        mutate(pheno.date = do.call(pmin, c(select(., 3:last_col()), na.rm = TRUE)),
               pheno.date = as.character(pheno.date),
               pheno.date = ifelse(is.na(pheno.final), NA, pheno.date)) %>%
        dplyr::select(f.eid, pheno.date)
    
    df.merged <- df.merged %>%
        left_join(df.pheno_date, by="f.eid")
    
    if (pheno_col != "NA"){
        df.merged <- df.merged %>%
            rename(!!as.name(pheno_col):=pheno.final,
                   !!as.name(paste0(pheno_col, ".Date")):=pheno.date)
    }
    
    return(df.merged)
}


#########################################################

# query_date_multiple_source <- function(queried_code_type, ICD10=NA, ICD9=NA, 
#                                         OPCS4=NA, Self_report=NA, simplify=FALSE){
#     # queried_code_type: a string or vector of strings of ICD10, ICD9, OPCS4, or Self-report
#     # Note: date with 3000-12-12 is considered as having disease but the exact date is unknown (according to coding 37 from field 87) from Self-report
#     #       - These should be considered as case. But if date is required, it should be considered as missing. 

#     n_quried <- length(queried_code_type)
#     for (idx in 1:n_quried){
#         if (queried_code_type[idx] == "ICD10"){
#             date_earliest_col <- colnames(ICD10)[grep("date_earliest*", colnames(ICD10))]
#             ICD10[[date_earliest_col]] <- as.Date(ICD10[[date_earliest_col]], format="%Y-%m-%d")
#         } else if (queried_code_type[idx] =="ICD9"){
#             date_earliest_col <- colnames(ICD9)[grep("date_earliest*", colnames(ICD9))]
#             ICD9[[date_earliest_col]] <- as.Date(ICD9[[date_earliest_col]], format="%Y-%m-%d")   
#         } else if (queried_code_type[idx] == "OPCS4"){
#             date_earliest_col <- colnames(OPCS4)[grep("date_earliest*", colnames(OPCS4))]
#             OPCS4[[date_earliest_col]] <- as.Date(OPCS4[[date_earliest_col]], format="%Y-%m-%d")   
#         } else if (queried_code_type[idx] == "Self-report"){
#             date_earliest_col <- colnames(Self_report)[grep("date_earliest*", colnames(Self_report))]
#             Self_report[[date_earliest_col]] <- as.Date(Self_report[[date_earliest_col]], format="%Y-%m-%d")        
#         }
#         if (idx == 1){
#             if (queried_code_type[idx] == "ICD10"){
#                 df <- ICD10
#             } else if (queried_code_type[idx] =="ICD9"){
#                 df <- ICD9
#             } else if (queried_code_type[idx] == "OPCS4"){
#                 df <- OPCS4
#             } else if (queried_code_type[idx] == "Self-report"){
#                 df <- Self_report
#             }
#         } else{
#             if (queried_code_type[idx] == "ICD10"){
#                 df <- merge(df, ICD10, by="f.eid", all=T)
#             } else if (queried_code_type[idx] =="ICD9"){
#                 df <- merge(df, ICD9, by="f.eid", all=T)
#             } else if (queried_code_type[idx] == "OPCS4"){
#                 df <- merge(df, OPCS4, by="f.eid", all=T)
#             } else if (queried_code_type[idx] == "Self-report"){
#                 df <- merge(df, Self_report, by="f.eid", all=T)
#             }
#         }
#     }

#     col_name <- colnames(df)[2:length(colnames(df))]
#     df$date_earliest <- apply(df[, col_name], 1, function(x) min(x, na.rm = TRUE))

#     date_columns <- sapply(df, inherits, "Date")
#     df <- df %>%
#         mutate(across(all_of(names(which(date_columns))), as.character))

#     if (simplify) df <- df %>% dplyr::select(f.eid, date_earliest)
    
#     return (df)
# }
