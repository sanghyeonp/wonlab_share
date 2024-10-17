library(dplyr)
library(tidyr)

######################################
### Preprocess Surgical Operation code ###
######################################
df <- readRDS("UKBB.20004.surgical_operation.rds")
df[is.na(df)] <- ""

col_inst0 <- paste0("f.20004.0.", c(0:31))
df.instance0_collapse <- df %>% dplyr::select(f.eid, all_of(col_inst0))
df.instance0_collapse$surgical_operation.inst0 <- apply(df.instance0_collapse[, 2:33], 1, function(x){paste(x, collapse=";")})
df.instance0_collapse <- df.instance0_collapse %>% dplyr::select(f.eid, surgical_operation.inst0)

col_inst1 <- paste0("f.20004.1.", c(0:31))
df.instance1_collapse <- df %>% dplyr::select(f.eid, all_of(col_inst1))
df.instance1_collapse$surgical_operation.inst1 <- apply(df.instance1_collapse[, 2:33], 1, function(x){paste(x, collapse=";")})
df.instance1_collapse <- df.instance1_collapse %>% dplyr::select(f.eid, surgical_operation.inst1)

col_inst2 <- paste0("f.20004.2.", c(0:31))
df.instance2_collapse <- df %>% dplyr::select(f.eid, all_of(col_inst2))
df.instance2_collapse$surgical_operation.inst2 <- apply(df.instance2_collapse[, 2:33], 1, function(x){paste(x, collapse=";")})
df.instance2_collapse <- df.instance2_collapse %>% dplyr::select(f.eid, surgical_operation.inst2)

col_inst3 <- paste0("f.20004.3.", c(0:31))
df.instance3_collapse <- df %>% dplyr::select(f.eid, all_of(col_inst3))
df.instance3_collapse$surgical_operation.inst3 <- apply(df.instance3_collapse[, 2:33], 1, function(x){paste(x, collapse=";")})
df.instance3_collapse <- df.instance3_collapse %>% dplyr::select(f.eid, surgical_operation.inst3)


df.age <- readRDS("UKBB.92.surgical_operation_age.rds")
df.age[is.na(df.age)] <- ""
df.birth_year <- readRDS("UKBB.34.year_of_birth.rds")

parse_date <- function(collapsed_data, birth_year){
    d <- strsplit(collapsed_data, ";")[[1]]
    parsed_d <- c()
    for (v in d){
        if (nchar(v) == 0){
            parsed_d <- c(parsed_d, "")
        } else{
            if (v %in% c("-1", "3")){
                parsed_d <- c(parsed_d, "3000-12-31")
            } else {
                if (nchar(v) == 4){
                    parsed_d <- c(parsed_d, paste(v, 12, 31, sep="-"))
                } else {
                    parsed_d <- c(parsed_d, paste(as.integer(birth_year) + as.integer(v), 12, 31, sep="-"))
                }
            }  
        }
    }
    return(paste(parsed_d, collapse=";"))
}

col_inst0 <- paste0("f.92.0.", c(0:31))
df.age.instance0_collapse <- df.age %>% dplyr::select(f.eid, all_of(col_inst0))
df.age.instance0_collapse$age_date.inst0 <- apply(df.age.instance0_collapse[, 2:33], 1, function(x){paste(x, collapse=";")})
df.age.instance0_collapse <- df.age.instance0_collapse %>% dplyr::select(f.eid, age_date.inst0)
df.age.instance0_collapse <- df.age.instance0_collapse %>% left_join(df.birth_year, by="f.eid")
df.age.instance0_collapse$date.inst0 <- apply(df.age.instance0_collapse, 1, function(row){parse_date(row[2], row[3])})
df.age.instance0_collapse <- df.age.instance0_collapse %>% dplyr::select(f.eid, date.inst0)

col_inst1 <- paste0("f.92.1.", c(0:31))
df.age.instance1_collapse <- df.age %>% dplyr::select(f.eid, all_of(col_inst1))
df.age.instance1_collapse$age_date.inst1 <- apply(df.age.instance1_collapse[, 2:33], 1, function(x){paste(x, collapse=";")})
df.age.instance1_collapse <- df.age.instance1_collapse %>% dplyr::select(f.eid, age_date.inst1)
df.age.instance1_collapse <- df.age.instance1_collapse %>% left_join(df.birth_year, by="f.eid")
df.age.instance1_collapse$date.inst1 <- apply(df.age.instance1_collapse, 1, function(row){parse_date(row[2], row[3])})
df.age.instance1_collapse <- df.age.instance1_collapse %>% dplyr::select(f.eid, date.inst1)

col_inst2 <- paste0("f.92.2.", c(0:31))
df.age.instance2_collapse <- df.age %>% dplyr::select(f.eid, all_of(col_inst2))
df.age.instance2_collapse$age_date.inst2 <- apply(df.age.instance2_collapse[, 2:33], 1, function(x){paste(x, collapse=";")})
df.age.instance2_collapse <- df.age.instance2_collapse %>% dplyr::select(f.eid, age_date.inst2)
df.age.instance2_collapse <- df.age.instance2_collapse %>% left_join(df.birth_year, by="f.eid")
df.age.instance2_collapse$date.inst2 <- apply(df.age.instance2_collapse, 1, function(row){parse_date(row[2], row[3])})
df.age.instance2_collapse <- df.age.instance2_collapse %>% dplyr::select(f.eid, date.inst2)

col_inst3 <- paste0("f.92.3.", c(0:31))
df.age.instance3_collapse <- df.age %>% dplyr::select(f.eid, all_of(col_inst3))
df.age.instance3_collapse$age_date.inst3 <- apply(df.age.instance3_collapse[, 2:33], 1, function(x){paste(x, collapse=";")})
df.age.instance3_collapse <- df.age.instance3_collapse %>% dplyr::select(f.eid, age_date.inst3)
df.age.instance3_collapse <- df.age.instance3_collapse %>% left_join(df.birth_year, by="f.eid")
df.age.instance3_collapse$date.inst3 <- apply(df.age.instance3_collapse, 1, function(row){parse_date(row[2], row[3])})
df.age.instance3_collapse <- df.age.instance3_collapse %>% dplyr::select(f.eid, date.inst3)


df.merged <- df.instance0_collapse %>%
    left_join(df.age.instance0_collapse, by="f.eid") %>%
    left_join(df.instance1_collapse, by="f.eid") %>%
    left_join(df.age.instance1_collapse, by="f.eid") %>%
    left_join(df.instance2_collapse, by="f.eid") %>%
    left_join(df.age.instance2_collapse, by="f.eid") %>%
    left_join(df.instance3_collapse, by="f.eid") %>%
    left_join(df.age.instance3_collapse, by="f.eid")

saveRDS(df.merged, "UKBB.20004_92.Surgical_operation_and_age_year_merged.rds")
