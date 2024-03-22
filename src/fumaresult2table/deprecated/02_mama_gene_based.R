library(data.table)
library(dplyr)

dir_out <- dir_out <- "/data1/sanghyeon/Projects/MetabolicSyndrome/MetS_2022_08/nat_genet.revision_1/reviewer_2/R2Q1/FUMA/MetS_GWAS_COJO/out"
n_thread <- 20

df <- fread(paste0(dir_out, "/magma.genes.out"), sep="\t", data.table = F, nThread = n_thread)

bonferroni_threshold <- 0.05 / nrow(df)

df <- df %>%
    mutate(`Bonferroni significant` = P < bonferroni_threshold) %>%
    rename(`Ensembl ID` = GENE,
            Gene = SYMBOL,
            Chr = CHR,
            Start = START,
            End = STOP,
            `P-value` = P,
            `N SNP` = NSNPS,
            `N Param` = NPARAM,
            `Z-score` = ZSTAT)

df <- df %>%
    mutate(Gene = paste0("`", Gene))

write.table(df,
            "FUMA.MAGMA_gene_based.csv",
            sep = ",", row.names = F, quote = T)
