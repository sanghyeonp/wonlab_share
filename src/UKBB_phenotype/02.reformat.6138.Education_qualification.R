library(data.table)

f <- "UKBB.6138.Education_qualification.rds"
df <- readRDS(f)

df1 <- df %>%
    dplyr::select(f.eid, f.6138.0.0, f.6138.0.1, f.6138.0.2, f.6138.0.3, f.6138.0.4, f.6138.0.5)

### Get maximum number of qualifications
df1$max_qual <- apply(df1[, -1], 1, function(x) {
    ifelse(max(x, na.rm = TRUE) %in% c(-Inf, -3), NA, max(x, na.rm = TRUE))
})

table(df1$max_qual, useNA = "always")

### Scoring
# Okbay (Nature, 2016; Supplementary Table 1.2 and 1.14)에서 정의된 기준에 따라 계산함.
# Okbay (Nature, 2016): Okbay, A., Beauchamp, J., Fontana, M. et al. Genome-wide association study identifies 74 loci associated with educational attainment. Nature 533, 539–542 (2016).

score_table <- c(
    "1" = 20,
    "2" = 13,
    "3" = 10, "4" = 10,
    "5" = 19,
    "6" = 15,
    "-7" = 7,
    "NA" = NA
)

df1$edu_year <- sapply(df1$max_qual, function(x){as.vector(score_table[as.character(x)])})

summary(df1$edu_year)

df2 <- df1 %>%
    dplyr::select(f.eid, edu_year)
write.table(df2, "UKBB.6138.instance_0.Education_years.txt", sep=" ", row.names=F, quote=F)
