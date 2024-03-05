
```
library(dplyr)
httr::set_config(httr::config(ssl_verifypeer = FALSE))
library(biomaRt)
```

# Step 1. Call Ensembl dataset
```
# GRCh38
ensembl <- useEnsembl(biomart = "genes", 
                    dataset = "hsapiens_gene_ensembl", 
                    mirror = "www")

# GRCh37
ensembl <- useEnsembl(biomart = "genes", 
                    dataset = "hsapiens_gene_ensembl", 
                    mirror = "www",
                    GRCh=37)
```

# Step 2. Get gene information for the specified genes (either gene symbol or Ensembl ID)
```
gene_list <- c()

df.gene_annot <- getBM(attributes = c("ensembl_gene_id",
                            "external_gene_name",
                            "gene_biotype",
                            "chromosome_name",
                            "band",
                            "start_position", 
                            "end_position",
                            "strand"),
                # if annotating based on gene symbol => "external_gene_name"
                # if annotating based on Ensembl ID => "ensembl_gene_id"
                filters = "external_gene_name",
                values = gene_list,
                mart = ensembl,
                uniqueRows = TRUE)

```

# Step 3. Filtering
```
df.gene_annot <- df.gene_annot %>%
    # Filter that doesn't have cytogenetic band
    # Reference: https://www.biostars.org/p/9571442/
    dplyr::filter(nchar(band) != 0) %>%
    # Filter invalid chromosome
    dplyr::filter(chromosome_name %in% as.character(c(1:22, "X", "Y")))
```
