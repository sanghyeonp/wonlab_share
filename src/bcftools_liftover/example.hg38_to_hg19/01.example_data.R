library(data.table)
library(dplyr)

f <- "/data1/sanghyeon/Projects/MetS_gSEM_EAS/data/gwas_EUR/BMI.EUR/EstBB_Abner2025/GCST90624699.tsv.gz"
df <- fread(f, data.table=F)

df$ID <- apply(df, 1, function(x){
    chr <- as.integer(x[["chromosome"]])
    bp <- as.integer(x[["base_pair_location"]])
    a1 <- x[["effect_allele"]]; a2 <- x[["other_allele"]]
    paste0(chr, ":", bp, ":", a1, ":", a2)
})

write.table(df, "GCST90624699.ID.tsv", sep="\t", row.names=F, quote=F)