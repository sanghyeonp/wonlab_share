library(dplyr)
library(tidyr)

######################################
### Preprocess Medication code ###
######################################
df <- readRDS("UKBB.6150.touchscreen_heart_problem.rds")
df[is.na(df)] <- ""

col_inst0 <- paste0("f.6150.0.", c(0:3))
df.instance0_collapse <- df %>% dplyr::select(f.eid, all_of(col_inst0))
df.instance0_collapse$touchscreen.inst0 <- apply(df.instance0_collapse[, 2:5], 1, function(x){paste(x, collapse=";")})
df.instance0_collapse <- df.instance0_collapse %>% dplyr::select(f.eid, touchscreen.inst0)

col_inst1 <- paste0("f.6150.1.", c(0:3))
df.instance1_collapse <- df %>% dplyr::select(f.eid, all_of(col_inst1))
df.instance1_collapse$touchscreen.inst1 <- apply(df.instance1_collapse[, 2:5], 1, function(x){paste(x, collapse=";")})
df.instance1_collapse <- df.instance1_collapse %>% dplyr::select(f.eid, touchscreen.inst1)

col_inst2 <- paste0("f.6150.2.", c(0:3))
df.instance2_collapse <- df %>% dplyr::select(f.eid, all_of(col_inst2))
df.instance2_collapse$touchscreen.inst2 <- apply(df.instance2_collapse[, 2:5], 1, function(x){paste(x, collapse=";")})
df.instance2_collapse <- df.instance2_collapse %>% dplyr::select(f.eid, touchscreen.inst2)

df.merged <- df.instance0_collapse %>%
    left_join(df.instance1_collapse, by="f.eid") %>%
    left_join(df.instance2_collapse, by="f.eid")

###### 
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
                parsed_d <- c(parsed_d, paste(as.integer(birth_year) + as.integer(v), 12, 31, sep="-"))
            }  
        }
    }
    return(paste(parsed_d, collapse=";"))
}

###### 
df.age.heart_attack <- readRDS("UKBB.3894.touchscreen_heartattack_age.rds")
df.age.heart_attack[is.na(df.age.heart_attack)] <- ""
df.age.heart_attack <- df.age.heart_attack %>%
    left_join(df.birth_year, by="f.eid")


df.age.heart.instance0 <- df.age.heart_attack %>% dplyr::select(f.eid, f.3894.0.0, f.34.0.0)
df.age.heart.instance0$date.heart_attack.inst0 <- apply(df.age.heart.instance0, 1, function(row){parse_date(row[2], row[3])})
df.age.heart.instance0 <- df.age.heart.instance0 %>% dplyr::select(f.eid, date.heart_attack.inst0)

df.age.heart.instance1 <- df.age.heart_attack %>% dplyr::select(f.eid, f.3894.1.0, f.34.0.0)
df.age.heart.instance1$date.heart_attack.inst1 <- apply(df.age.heart.instance1, 1, function(row){parse_date(row[2], row[3])})
df.age.heart.instance1 <- df.age.heart.instance1 %>% dplyr::select(f.eid, date.heart_attack.inst1)

df.age.heart.instance2 <- df.age.heart_attack %>% dplyr::select(f.eid, f.3894.2.0, f.34.0.0)
df.age.heart.instance2$date.heart_attack.inst2 <- apply(df.age.heart.instance2, 1, function(row){parse_date(row[2], row[3])})
df.age.heart.instance2 <- df.age.heart.instance2 %>% dplyr::select(f.eid, date.heart_attack.inst2)

df.merged <- df.merged %>%
    left_join(df.age.heart.instance0, by="f.eid") %>%
    left_join(df.age.heart.instance1, by="f.eid") %>%
    left_join(df.age.heart.instance2, by="f.eid")


###### 
df.age.angina <- readRDS("UKBB.3627.touchscreen_angina_age.rds")
df.age.angina[is.na(df.age.angina)] <- ""
df.age.angina <- df.age.angina %>%
    left_join(df.birth_year, by="f.eid")

df.age.angina.instance0 <- df.age.angina %>% dplyr::select(f.eid, f.3627.0.0, f.34.0.0)
df.age.angina.instance0$date.angina.inst0 <- apply(df.age.angina.instance0, 1, function(row){parse_date(row[2], row[3])})
df.age.angina.instance0 <- df.age.angina.instance0 %>% dplyr::select(f.eid, date.angina.inst0)

df.age.angina.instance1 <- df.age.angina %>% dplyr::select(f.eid, f.3627.1.0, f.34.0.0)
df.age.angina.instance1$date.angina.inst1 <- apply(df.age.angina.instance1, 1, function(row){parse_date(row[2], row[3])})
df.age.angina.instance1 <- df.age.angina.instance1 %>% dplyr::select(f.eid, date.angina.inst1)

df.age.angina.instance2 <- df.age.angina %>% dplyr::select(f.eid, f.3627.2.0, f.34.0.0)
df.age.angina.instance2$date.angina.inst2 <- apply(df.age.angina.instance2, 1, function(row){parse_date(row[2], row[3])})
df.age.angina.instance2 <- df.age.angina.instance2 %>% dplyr::select(f.eid, date.angina.inst2)

df.merged <- df.merged %>%
    left_join(df.age.angina.instance0, by="f.eid") %>%
    left_join(df.age.angina.instance1, by="f.eid") %>%
    left_join(df.age.angina.instance2, by="f.eid")


saveRDS(df.merged, "UKBB.6150_3894_3627.heart_problem_date_merged.rds")
