
# SNP genomic position annotation
```
library(biomaRt)
rsids <- df$rsID

# Connect to Ensembl GRCh37
mart37 <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp", GRCh = 37)

snp_annot <- getBM(
    attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end", "allele"),
    filters    = "snp_filter",
    values     = rsids,
    mart       = mart37
)

snp_annot <- snp_annot %>%
    filter(chr_name %in% as.character(1:22, "X", "Y"))

```

# skimr::skim
```
library(skimr)
df <- data.frame(A = c(1:10), B = c(rep("A", 5), rep("B", 5)))
skim(df)

```

# dput
```
df <- data.frame(A = c(1:10), B = c(rep("A", 5), rep("B", 5)))
dput(names(df))
# Output
"A", "B"
```
