library(data.table)
library(dplyr)

df <- fread("FUMA.gene_mapping.combined.csv", sep=",",
            data.table = F, nThread = 20)

df.magma <- fread("FUMA.MAGMA_gene_based.csv", sep=",",
                data.table = F, nThread = 20) %>%
    filter(`Bonferroni significant` == TRUE) %>%
    dplyr::select(`Ensembl ID`) %>%
    mutate(MAGMA = "Yes")

df <- merge(df, df.magma, by = "Ensembl ID", all = T) %>%
    tidyr::replace_na(list(MAGMA = "No"))



write.table(df,
            "FUMA.gene_mapping.MAGMA_gene_based.csv",
            sep = ",", row.names = F, quote = T)

print(length(unique((df %>% filter(posMap == "Yes" &
                        eqtlMap == "Yes" &
                        ciMap == "Yes" &
                        MAGMA == "Yes"))$`Ensembl ID`)))
write.table(df %>% filter(posMap == "Yes" &
                        eqtlMap == "Yes" &
                        ciMap == "Yes" &
                        MAGMA == "Yes"),
            "FUMA.gene_mapping.MAGMA_gene_based.all_mapped.csv",
            sep = ",", row.names = F, quote = T)