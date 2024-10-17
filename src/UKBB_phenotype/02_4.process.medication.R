library(dplyr)
library(tidyr)

######################################
### Preprocess Medication code ###
######################################
df <- readRDS("/data1/sanghyeon/wonlab_contribute/combined/src/UKBB_phenotype/UKBB.20003.medication.rds")
df[is.na(df)] <- ""

col_inst0 <- paste0("f.20003.0.", c(0:47))
df.instance0_collapse <- df %>% dplyr::select(f.eid, all_of(col_inst0))
df.instance0_collapse$medication.inst0 <- apply(df.instance0_collapse[, 2:49], 1, function(x){paste(x, collapse=";")})
df.instance0_collapse <- df.instance0_collapse %>% dplyr::select(f.eid, medication.inst0)

col_inst1 <- paste0("f.20003.1.", c(0:47))
df.instance1_collapse <- df %>% dplyr::select(f.eid, all_of(col_inst1))
df.instance1_collapse$medication.inst1 <- apply(df.instance1_collapse[, 2:49], 1, function(x){paste(x, collapse=";")})
df.instance1_collapse <- df.instance1_collapse %>% dplyr::select(f.eid, medication.inst1)

col_inst2 <- paste0("f.20003.2.", c(0:47))
df.instance2_collapse <- df %>% dplyr::select(f.eid, all_of(col_inst2))
df.instance2_collapse$medication.inst2 <- apply(df.instance2_collapse[, 2:49], 1, function(x){paste(x, collapse=";")})
df.instance2_collapse <- df.instance2_collapse %>% dplyr::select(f.eid, medication.inst2)

col_inst3 <- paste0("f.20003.3.", c(0:47))
df.instance3_collapse <- df %>% dplyr::select(f.eid, all_of(col_inst3))
df.instance3_collapse$medication.inst3 <- apply(df.instance3_collapse[, 2:49], 1, function(x){paste(x, collapse=";")})
df.instance3_collapse <- df.instance3_collapse %>% dplyr::select(f.eid, medication.inst3)


df.merged <- df.instance0_collapse %>%
    left_join(df.instance1_collapse, by="f.eid") %>%
    left_join(df.instance2_collapse, by="f.eid") %>%
    left_join(df.instance3_collapse, by="f.eid")

saveRDS(df.merged, "UKBB.20003.medication_processed.rds")
