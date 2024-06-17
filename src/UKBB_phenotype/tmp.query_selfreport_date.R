df.data <- readRDS("UKBB.20002_87.Selfreport_code_and_age_year_merged.rds")

code <- "1065"
self_report_instance <- 0

reformat_into_search_pattern <- function(code, exact_match, self_report=FALSE){
    if (exact_match){
        if (!self_report){
            # Check if looking for all subclassifications
            all_subclass <- !grepl("\\.", code)
            # Remove "." from the code
            code <- gsub("\\.", "", code)
            # Add "*" if looking for all subclassifications, otherwise add "$"
            if(all_subclass) code <- paste0(code, "\\d*") else code <- paste0(code, "$")
            # Start with the given code itself
            code <- paste0("^", code)
        } else{
            code <- paste0("^", code, "$")
        }
    } else{
        code <- gsub("\\.", "", code)
    }
    return (code)
}

self_report_instance <- as.character(self_report_instance)
search_pattern <- reformat_into_search_pattern(code, exact=F, self_report=T)
search_pattern
exact_pattern <- reformat_into_search_pattern(code, exact=T, self_report=T)
exact_pattern

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

# Query the ICD9 codes
df.query <- as.data.frame(df.data %>%
              dplyr::rowwise() %>%
              dplyr::mutate(
                  # Index of ICD9 code of interest
                  index = ifelse(grepl(search_pattern, !!as.name(code_col)), 
                                 paste(which(grepl(exact_pattern, strsplit(!!as.name(code_col), ";")[[1]])), collapse=";"), NA), 
                  # Query ICD9 code of interest using the index (double check if correct code is queried)
                  code_query = ifelse(!is.na(index), 
                                      paste(strsplit(!!as.name(code_col), ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
                  # Query date of diagnosis of interest using the index
                  date_query = ifelse(!is.na(index), 
                                      paste(strsplit(!!as.name(date_col), ";")[[1]][as.integer(strsplit(index, ";")[[1]])], collapse=";"), NA),
                  # Obtain the earliest date of diagnosis
                  date_earliest = ifelse(is.na(date_query),
                                         NA, ifelse(!grepl(";", date_query), 
                                                    date_query,
                                                    as.character(min(as.Date(strsplit(date_query, ";")[[1]], format="%Y-%m-%d")))))) %>%
              # Drop collapsed data columns (dplyr::select(-!!as.name(code_col), -!!as.name(date_col)))
              dplyr::select(f.eid, index, code_query, date_query, date_earliest) %>%
              dplyr::rename(!!as.name(paste0("index.", suffix)) := index,
                        !!as.name(paste0("code_query.", suffix)) := code_query,
                        !!as.name(paste0("date_query.", suffix)) := date_query,
                        !!as.name(paste0("date_earliest.", suffix)) := date_earliest)
              )

df.data[37015, ]
grepl(search_pattern, df.data[2, "selfreport_inst0"])
grepl(exact_pattern, strsplit(df.data[2, "selfreport_inst0"], ";")[[1]])
