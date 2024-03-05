
# bioMart: GRCh38
```
library(dplyr)
httr::set_config(httr::config(ssl_verifypeer = FALSE))
library(biomaRt)

snp_list <- c()

ensembl <- useEnsembl(biomart = "genes", 
                    dataset = "hsapiens_gene_ensembl", 
                    mirror = "www")

mapping <- getBM(attributes = c("ensembl_gene_id",
                            "external_gene_name",
                            "gene_biotype",
                            "chromosome_name",
                            "band",
                            "start_position", 
                            "end_position",
                            "strand"),
                filters = "external_gene_name",
                values = snp_list,
                mart = ensembl,
                uniqueRows = TRUE)

mapping <- mapping %>%
    # Filter that doesn't have cytogenetic band
    # Reference: https://www.biostars.org/p/9571442/
    dplyr::filter(nchar(band) != 0)
```

# bioMart: GRCh37
```
library(dplyr)
httr::set_config(httr::config(ssl_verifypeer = FALSE))
library(biomaRt)

snp_list <- c()

ensembl <- useEnsembl(biomart = "genes", 
                    dataset = "hsapiens_gene_ensembl", 
                    mirror = "www",
                    GRCh=37)

mapping <- getBM(attributes = c("ensembl_gene_id",
                            "external_gene_name",
                            "gene_biotype",
                            "chromosome_name",
                            "band",
                            "start_position", 
                            "end_position",
                            "strand"),
                filters = "external_gene_name",
                values = snp_list,
                mart = ensembl,
                uniqueRows = TRUE)

mapping <- mapping %>%
    # Filter that doesn't have cytogenetic band
    # Reference: https://www.biostars.org/p/9571442/
    dplyr::filter(nchar(band) != 0)
```
