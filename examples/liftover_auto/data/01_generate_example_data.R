library(data.table)
library(dplyr)

file <- "finngen_R8_T2D.autosome.10k"

df <- fread(file)

df1 <- df %>% select("#chrom", pos, ref, alt, rsids, beta, sebeta, pval)

write.table(df1, "finngen_R8_T2D.autosome.10k.subset",
            sep="\t", row.names=FALSE, quote=FALSE)

## chr:pos:ref:alt
df2 <- df1 %>% 
    rename(CHR := all_of("#chrom")) %>%
    mutate(variant=paste(CHR, pos, ref, alt, sep=":")) %>%
    select(-CHR, -pos, -ref, -alt, -rsids)

write.table(df2, "finngen_R8_T2D.autosome.10k.subset.infer",
            sep="\t", row.names=FALSE, quote=FALSE)
