library(data.table)
library(dplyr)

df.hg38 <- fread("/data1/sanghyeon/wonlab_contribute/combined/src/bcftools_liftover/full_example.hg38_to_hg19.no_del_multiallelic/GCST90624699.ID.tsv", data.table=F)
df.hg38.sub <- df.hg38 %>%
    dplyr::select(
        ID, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta
    ) %>%
    rename(
        chr.hg38=chromosome, bp.hg38=base_pair_location,
        ea.hg38=effect_allele, oa.hg38=other_allele, eaf.hg38=effect_allele_frequency,
        b.hg38=beta
    )
message("N SNP initial: ", nrow(df.hg38.sub))

df.hg38.sub.dup <- df.hg38.sub[duplicated(df.hg38.sub$ID) | duplicated(df.hg38.sub$ID, fromLast = TRUE), ]

df.hg19 <- fread("BMI_EUR.Abner2024.hg19.tsv", data.table=F)
df.hg19.sub <- df.hg19 %>%
    dplyr::select(
        ID, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta
    ) %>%
    rename(
        chr.hg19=chromosome, bp.hg19=base_pair_location,
        ea.hg19=effect_allele, oa.hg19=other_allele, eaf.hg19=effect_allele_frequency,
        b.hg19=beta
    )

df.hg19.sub.dup <- df.hg19.sub[duplicated(df.hg19.sub$ID) | duplicated(df.hg19.sub$ID, fromLast = TRUE), ]

## Allele swap, strand flip
df.swap <- fread("BMI_EUR.Abner2024.allele_swap.strand_flip.tsv", sep="\t", data.table=F)
# df.swap.dup <- df.swap[duplicated(df.swap$ID) | duplicated(df.swap$ID, fromLast = TRUE), ]

df.swap <- df.swap %>%
    dplyr::select(ID, SRC_REF_ALT, SWAP, FLIP)

df <- df.hg38.sub %>%
    inner_join(df.hg19.sub, by="ID") %>%
    left_join(df.swap, by="ID")
message("N SNPs lifted (hg38->hg19): ", nrow(df))

saveRDS(df, "df.rds")

## SWAP
df.swapped <- df %>%
    filter(SWAP=="1")
head(df.swapped)
