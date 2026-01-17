library(data.table)

f <- "UKBB.738.Household_income.rds"
df <- readRDS(f)

table(df$f.738.0.0)
### Data coding: https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=100294
# 1	Less than 18,000
# 2	18,000 to 30,999
# 3	31,000 to 51,999
# 4	52,000 to 100,000
# 5	Greater than 100,000
# -1	Do not know
# -3	Prefer not to answer

df1 <- df %>%
    dplyr::select(f.eid, f.738.0.0) %>%
    rename(Household_income = f.738.0.0) %>%
    mutate(Household_income = ifelse(Household_income %in% c(-1, -3), NA, Household_income))

table(df1$Household_income)

write.table(df1, "UKBB.738.instance_0.Household_income.txt", sep=" ", row.names=F, quote=F)
