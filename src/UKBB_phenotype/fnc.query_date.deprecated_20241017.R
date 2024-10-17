library(dplyr)

query_date <- function(code_type, code, self_report_instance=0, simplify=TRUE){
    # code_type: either "ICD10", "ICD9", "OPCS4", "Self-report" or "Cause of death"
    # code: single or a vector character of ICD-10, ICD-9, or OPCS4 or Self-report codes
    # self_report_instance: either all, 0, 1, 2, 3
    # simplify: if TRUE, return f.eid and only the date of the earliest diagnosis,

    # Elapsed time: average 1 minutes per code type
    if (code_type == "ICD10"){
        df.query <- query_date_ICD10(code=code, simplify=simplify)
    } else if (code_type == "ICD9"){
        df.query <- query_date_ICD9(code=code, simplify=simplify)
    } else if (code_type == "OPCS4"){
        df.query <- query_date_OPCS4(code=code, simplify=simplify)
    } else if (code_type == "Self-report") {
        df.query <- query_date_Self_report(code=code, self_report_instance=self_report_instance, simplify=simplify)
    } else if (code_type == "Cause of death"){
        df.query <- query_date_cause_of_death(code=code, simplify=simplify)
    } else {
        stop("Invalid code_type. Please choose from 'ICD10', 'ICD9', 'OPCS4', 'Self-report' or 'Cause of death'.")
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
query_date_ICD10 <- function(code, simplify=TRUE){
    # code: single or a vector character of ICD-10 codes
    #   - Example) c("I21", "I24", "I25.2", "I46")
    # simplify: if TRUE, return f.eid and only the date of the earliest diagnosis, 
    #            otherwise return f.eid, index, quried ICD-10 code, quried date of diagnosis, and the earliest date of diagnosis

    # Read the data
    df.data <- readRDS("/data1/sanghyeon/wonlab_contribute/combined/src/UKBB_phenotype/UKBB.41270_41280.ICD10_code_and_date_merged.rds")

    # Make the search pattern of given ICD-10 codes
    code <- as.character(code)
    search_pattern <- paste(paste0("^", gsub("\\.", "", code), collapse="|"))
    # search_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=FALSE, self_report=FALSE), collapse = "|")
    # exact_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=TRUE, self_report=FALSE), collapse = "|")
    cat(paste0("Specified ICD10: ", paste(code, collapse=", "), "\n", 
            "Search pattern: ", search_pattern, "\n"))

    # Query the ICD10 codes
    df.query <- as.data.frame(df.data %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            # Index of ICD10 code of interest
            index = ifelse(sum(grepl(search_pattern, strsplit(ICD10_main, ";")[[1]])) > 0, 
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
query_date_ICD9 <- function(code, simplify=TRUE){
    # code: single or a vector character of ICD-9 codes
    #   - Example) c("410", "411", "414.2", "427.5")
    # simplify: if TRUE, return f.eid and only the date of the earliest diagnosis, 
    #            otherwise return f.eid, index, quried ICD-9 code, quried date of diagnosis, and the earliest date of diagnosis

    # Read the data
    df.data <- readRDS("/data1/sanghyeon/wonlab_contribute/combined/src/UKBB_phenotype/UKBB.41271_41281.ICD9_code_and_date_merged.rds")

    # Make the search pattern of given ICD-9 codes
    code <- as.character(code)
    search_pattern <- paste(paste0("^", gsub("\\.", "", code), collapse="|"))
    # search_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=FALSE, self_report=FALSE), collapse = "|")
    # exact_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=TRUE, self_report=FALSE), collapse = "|")
    cat(paste0("Specified ICD9: ", paste(code, collapse=", "), "\n", 
            "Search pattern: ", search_pattern, "\n"))

    # Query the ICD9 codes
    df.query <- as.data.frame(df.data %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            # Index of ICD9_main code of interest
            index = ifelse(sum(grepl(search_pattern, strsplit(ICD9_main, ";")[[1]])) > 0, 
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
query_date_OPCS4 <- function(code, simplify=TRUE){
    # code: single or a vector character of OPCS4 codes
    #   - Example) c("K40", "K41", "K42", "K44.1", "K44.2")
    # simplify: if TRUE, return f.eid and only the date of the earliest diagnosis, 
    #            otherwise return f.eid, index, quried OPCS4 code, quried date of diagnosis, and the earliest date of diagnosis

    # Read the data
    df.data <- readRDS("/data1/sanghyeon/wonlab_contribute/combined/src/UKBB_phenotype/UKBB.41272_41282.OPCS4_code_and_date_merged.rds")

    # Make the search pattern of given OPCS4 codes
    code <- as.character(code)
    search_pattern <- paste(paste0("^", gsub("\\.", "", code), collapse="|"))
    # search_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=FALSE, self_report=FALSE), collapse = "|")
    # exact_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=TRUE, self_report=FALSE), collapse = "|")
    cat(paste0("Specified OPCS4: ", paste(code, collapse=", "), "\n", 
            "Search pattern: ", search_pattern, "\n"))

    # Query the OPCS4 codes
    df.query <- as.data.frame(df.data %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            # Index of OPCS4 code of interest
            index = ifelse(sum(grepl(search_pattern, strsplit(OPCS4, ";")[[1]])) > 0, 
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
query_date_Self_report <- function(code, self_report_instance=0, simplify=TRUE){
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
    search_pattern <- paste(paste0("^", gsub("\\.", "", code), collapse="|"))
    # search_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=FALSE, self_report=FALSE), collapse = "|")
    # exact_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=TRUE, self_report=FALSE), collapse = "|")
    cat(paste0("Specified self-report codes: ", paste(code, collapse=", "), "\n", 
            "Search pattern: ", search_pattern, "\n"))
    
    if (self_report_instance == "all"){
        instance_to_retain <- c("f.eid", names(df.data)[grepl("all", names(df.data))])
        df.data <- df.data %>% select(all_of(instance_to_retain))
        suffix <- "selfreport_all"
    } else{
        instance_to_retain <- c("f.eid", names(df.data)[grepl(paste0("inst", self_report_instance), names(df.data))])
        df.data <- df.data %>% select(all_of(instance_to_retain))
        suffix <- paste0("selfreport_inst", self_report_instance)
    }

    code_col <- colnames(df.data)[2]; date_col <- colnames(df.data)[4]

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


###########################################
# Query the death date for cause of death #
###########################################
query_date_cause_of_death <- function(code, simplify=TRUE){
    # code: single or a vector character of ICD-10 codes
    #   - Example) c("I21", "I24", "I25.2", "I46")
    # simplify: if TRUE, return f.eid and only the date of the earliest diagnosis, 
    #            otherwise return f.eid, index, quried ICD-10 code for death cause and quried death date

    # Read the data
    df.data <- readRDS("/data1/sanghyeon/wonlab_contribute/combined/src/UKBB_phenotype/UKBB.40000_40001.Cause_of_death_and_date_merged.rds")

    # Make the search pattern of given ICD-10 codes
    code <- as.character(code)
    search_pattern <- paste(paste0("^", gsub("\\.", "", code), collapse="|"))
    # search_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=FALSE, self_report=FALSE), collapse = "|")
    # exact_pattern <- paste(sapply(code, reformat_into_search_pattern, exact_match=TRUE, self_report=FALSE), collapse = "|")
    cat(paste0("Specified ICD10: ", paste(code, collapse=", "), "\n", 
            "Search pattern: ", search_pattern, "\n"))

    # Query the ICD10 codes
    df.query <- as.data.frame(df.data %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
            # Index of ICD10 code of interest
            index = ifelse(sum(grepl(search_pattern, strsplit(Cause_of_death, ";")[[1]])) > 0, 
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


#########################################################

query_date_multiple_source <- function(queried_code_type, ICD10=NA, ICD9=NA, 
                                        OPCS4=NA, Self_report=NA, simplify=FALSE){
    # queried_code_type: a string or vector of strings of ICD10, ICD9, OPCS4, or Self-report
    # Note: date with 3000-12-12 is considered as having disease but the exact date is unknown (according to coding 37 from field 87) from Self-report
    #       - These should be considered as case. But if date is required, it should be considered as missing. 

    n_quried <- length(queried_code_type)
    for (idx in 1:n_quried){
        if (queried_code_type[idx] == "ICD10"){
            date_earliest_col <- colnames(ICD10)[grep("date_earliest*", colnames(ICD10))]
            ICD10[[date_earliest_col]] <- as.Date(ICD10[[date_earliest_col]], format="%Y-%m-%d")
        } else if (queried_code_type[idx] =="ICD9"){
            date_earliest_col <- colnames(ICD9)[grep("date_earliest*", colnames(ICD9))]
            ICD9[[date_earliest_col]] <- as.Date(ICD9[[date_earliest_col]], format="%Y-%m-%d")   
        } else if (queried_code_type[idx] == "OPCS4"){
            date_earliest_col <- colnames(OPCS4)[grep("date_earliest*", colnames(OPCS4))]
            OPCS4[[date_earliest_col]] <- as.Date(OPCS4[[date_earliest_col]], format="%Y-%m-%d")   
        } else if (queried_code_type[idx] == "Self-report"){
            date_earliest_col <- colnames(Self_report)[grep("date_earliest*", colnames(Self_report))]
            Self_report[[date_earliest_col]] <- as.Date(Self_report[[date_earliest_col]], format="%Y-%m-%d")        
        }
        if (idx == 1){
            if (queried_code_type[idx] == "ICD10"){
                df <- ICD10
            } else if (queried_code_type[idx] =="ICD9"){
                df <- ICD9
            } else if (queried_code_type[idx] == "OPCS4"){
                df <- OPCS4
            } else if (queried_code_type[idx] == "Self-report"){
                df <- Self_report
            }
        } else{
            if (queried_code_type[idx] == "ICD10"){
                df <- merge(df, ICD10, by="f.eid", all=T)
            } else if (queried_code_type[idx] =="ICD9"){
                df <- merge(df, ICD9, by="f.eid", all=T)
            } else if (queried_code_type[idx] == "OPCS4"){
                df <- merge(df, OPCS4, by="f.eid", all=T)
            } else if (queried_code_type[idx] == "Self-report"){
                df <- merge(df, Self_report, by="f.eid", all=T)
            }
        }
    }

    col_name <- colnames(df)[2:length(colnames(df))]
    df$date_earliest <- apply(df[, col_name], 1, function(x) min(x, na.rm = TRUE))

    date_columns <- sapply(df, inherits, "Date")
    df <- df %>%
        mutate(across(all_of(names(which(date_columns))), as.character))

    if (simplify) df <- df %>% dplyr::select(f.eid, date_earliest)
    
    return (df)
}
