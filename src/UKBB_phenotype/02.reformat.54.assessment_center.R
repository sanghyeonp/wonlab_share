library(data.table)
library(dplyr)

f <- "UKBB.54.assessment_center.rds"
df <- readRDS(f)

### Data coding: https://biobank.ctsu.ox.ac.uk/ukb/coding.cgi?id=10
# 10003  Stockport (pilot)
# 11001  Manchester
# 11002  Oxford
# 11003  Cardiff
# 11004  Glasgow
# 11005  Edinburgh
# 11006  Stoke
# 11007  Reading
# 11008  Bury
# 11009  Newcastle
# 11010  Leeds
# 11011  Bristol
# 11012  Barts
# 11013  Nottingham
# 11014  Sheffield
# 11016  Liverpool
# 11017  Middlesbrough
# 11018  Hounslow
# 11020  Croydon
# 11021  Birmingham
# 11022  Swansea
# 11023  Wrexham
# 11024  Cheadle (revisit)
# 11025  Cheadle (imaging)
# 11026  Reading (imaging)
# 11027  Newcastle (imaging)
# 11028  Bristol (imaging)


code_map <- c(
    "10003" = 1, "11001" = 2, "11002" = 3, "11003" = 4, "11004" = 5,
    "11005" = 6, "11006" = 7, "11007" = 8, "11008" = 9, "11009" = 10,
    "11010" = 11, "11011" = 12, "11012" = 13, "11013" = 14, "11014" = 15,
    "11016" = 16, "11017" = 17, "11018" = 18, "11020" = 19, "11021" = 20,
    "11022" = 21, "11023" = 22, "11024" = 23, "11025" = 24, "11026" = 25,
    "11027" = 26, "11028" = 27
)
df$assessment_center <- sapply(df$`f.54.0.0`, function(x){code_map[as.character(x)]})

df1 <- df %>% dplyr::select(f.eid, assessment_center)

write.table(df1, "UKBB.54.instance_0.Assessment_center.txt", sep=" ", row.names=F, quote=F)