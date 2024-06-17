library(data.table)
library(dplyr)

file.case <- "UKBB.IA_case.txt"
file.control_exclude <- "UKBB.IA_IA_control_exclude.txt"

df.case <- fread(file.case, sep=" ", data.table=F)
df.control_exclude <- fread(file.control_exclude, sep=" ", data.table=F) %>%
    dplyr::select(FID, IID, IA_control_exclude)

df <- merge(df.case, df.control_exclude, by=c("FID", "IID"), all=T)

df <- df %>%
    mutate(IA.final = ifelse(is.na(IA), NA,
                            ifelse(IA == 1, 1, 
                            ifelse(IA == 0 & IA_control_exclude == 0, 0, NA))))

write.table(df,
            "UKBB.phenotype.IA.txt",
            sep=" ", quote=F, row.names=F)
