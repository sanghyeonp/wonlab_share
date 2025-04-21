library(dplyr)
library(tidyr)

df.assessment_date <- readRDS("UKBB.53.assessment_center_date.rds")


######################################
### Preprocess cause of death and date ###
######################################
# print("Preprocess cause of death and date")
# start_time <- Sys.time()
# df.cause_of_death <- readRDS("UKBB.40001.cause_of_death.rds")
# df.death_date <- readRDS("UKBB.40000.date_of_death.rds")

# df.cause_of_death.instance_collapse <- as.data.frame(df.cause_of_death %>%
#         dplyr::select(f.eid, "f.40001.0.0", "f.40001.1.0", "f.40001.2.0") %>%
#         rowwise() %>%
#         mutate(Cause_of_death = paste(c_across(-f.eid), collapse = ";")) %>%
#         ungroup() %>%
#         dplyr::select(f.eid, Cause_of_death))

# idate_columns <- sapply(df.death_date, inherits, "IDate")
# df.death_date.instance_collapse <- as.data.frame(df.death_date %>%
#     mutate(across(all_of(names(which(idate_columns))), as.character)) %>%
#     mutate_all(~replace_na(., "")) %>%
#     dplyr::select(f.eid, "f.40000.0.0", "f.40000.1.0", "f.40000.2.0") %>%
#     rowwise() %>%
#     mutate(Cause_of_death_date = paste(c_across(-f.eid), collapse = ";")) %>%
#     ungroup() %>%
#     dplyr::select(f.eid, Cause_of_death_date))

# df.cause_of_death_and_date <- merge(df.cause_of_death.instance_collapse, 
#                                     df.death_date.instance_collapse, by="f.eid", all=T)

# saveRDS(df.cause_of_death_and_date, "UKBB.40000_40001.Cause_of_death_and_date_merged.rds")

# end_time <- Sys.time()
# elapsed_time <- end_time - start_time
# print(elapsed_time)


######################################
### Preprocess cause of secondary death and date ###
#! The instances are separated by ";" and the values are separated by ","
######################################

print("Preprocess cause of secondary death and date")
start_time <- Sys.time()
df.cause_of_death2 <- readRDS("UKBB.40002.secondary_cause_of_death.rds")
df.death_date <- readRDS("UKBB.40000.date_of_death.rds")

df.cause_of_death2.instance0 <- as.data.frame(df.cause_of_death2 %>%
        dplyr::select(f.eid, all_of(paste0("f.40002.0.", 0:14))) %>%
        rowwise() %>%
        mutate(Secondary_cause_of_death.inst0 = paste(c_across(-f.eid), collapse = ",")) %>%
        ungroup() %>%
        dplyr::select(f.eid, Secondary_cause_of_death.inst0))
df.cause_of_death2.instance1 <- as.data.frame(df.cause_of_death2 %>%
        dplyr::select(f.eid, all_of(paste0("f.40002.1.", 0:14))) %>%
        rowwise() %>%
        mutate(Secondary_cause_of_death.inst1 = paste(c_across(-f.eid), collapse = ",")) %>%
        ungroup() %>%
        dplyr::select(f.eid, Secondary_cause_of_death.inst1))
df.cause_of_death2.instance2 <- as.data.frame(df.cause_of_death2 %>%
        dplyr::select(f.eid, all_of(paste0("f.40002.2.", 0:13))) %>%
        rowwise() %>%
        mutate(Secondary_cause_of_death.inst2 = paste(c_across(-f.eid), collapse = ",")) %>%
        ungroup() %>%
        dplyr::select(f.eid, Secondary_cause_of_death.inst2))

df.cause_of_death2.instance_collapse <- df.cause_of_death2.instance0 %>%
    left_join(df.cause_of_death2.instance1, by="f.eid") %>%
    left_join(df.cause_of_death2.instance2, by="f.eid") %>%
    mutate(Secondary_cause_of_death = paste(Secondary_cause_of_death.inst0,
                                            Secondary_cause_of_death.inst1,
                                            Secondary_cause_of_death.inst2,
                                            sep = ";")) %>%
    dplyr::select(f.eid, Secondary_cause_of_death)


idate_columns <- sapply(df.death_date, inherits, "IDate")
df.death_date.instance_collapse <- as.data.frame(df.death_date %>%
    mutate(across(all_of(names(which(idate_columns))), as.character)) %>%
    mutate_all(~replace_na(., "")) %>%
    dplyr::select(f.eid, "f.40000.0.0", "f.40000.1.0", "f.40000.2.0") %>%
    rowwise() %>%
    mutate(Secondary_cause_of_death_date = paste(c_across(-f.eid), collapse = ";")) %>%
    ungroup() %>%
    dplyr::select(f.eid, Secondary_cause_of_death_date))

df.secondary_cause_of_death_and_date <- merge(df.cause_of_death2.instance_collapse, 
                                    df.death_date.instance_collapse, by="f.eid", all=T)

saveRDS(df.secondary_cause_of_death_and_date, "UKBB.40002_40001.Secondary_cause_of_death_and_date_merged.rds")

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)
