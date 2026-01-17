library(biomaRt)
library(dplyr)

######## 1. biomaRt object 만들기
#: Ensembl version 및 host 확인 필요 (참고 1.1. Ensembl version 별, host 확인법)

# 아래 예제는 FUMA에서 제공하는 v102.
ensembl <- useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl",
    host = "https://nov2020.archive.ensembl.org"
)

######## 1.1. Ensembl version 별, host 확인법
my_version <- "102"

df.version_host <- biomaRt::listEnsemblArchives()
df.version_host %>% filter(version == my_version)


######## 2. 
df.attr <- listAttributes(ensembl)

df.mapper <- getBM(
    attributes = c(
        "ensembl_gene_id",
        "hgnc_symbol",
        "hgnc_id",
        "entrezgene_id",
        "external_gene_name",
        "external_gene_source",
        "external_synonym",
        "chromosome_name",
        "start_position",
        "end_position",
        "gene_biotype"
    ),
    mart = ensembl
) # nRow = 106,004

######## 3. HGNC prefix 제거
df.mapper1 <- df.mapper %>%
    mutate(hgnc_id = gsub("HGNC:", "", hgnc_id)) # nRow = 106,004

df.mapper1[df.mapper1 == "" | df.mapper1 == "NA"] <- NA

######## 4. autosome 및 sex chromosome만 남기기
# Filter chromosome 1:22, X, and Y
unique(df.mapper$chromosome_name)
df.mapper2 <- df.mapper1 %>%
    filter(chromosome_name %in% as.character(c(1:22, "X", "Y"))) # nRow = 91,811


######## 5. Unique한 유전자만 남기기.
dup.ensemblid1 <- df.mapper2[duplicated(df.mapper2$ensembl_gene_id) | 
                                 duplicated(df.mapper2$ensembl_gene_id, fromLast = TRUE), ] # nRow = 91,811
dup.hgnc_symbol1 <- df.mapper2[duplicated(df.mapper2$hgnc_symbol) | 
                                   duplicated(df.mapper2$hgnc_symbol, fromLast = TRUE), ] # nRow = 67,372
dup.hgnc_id1 <- df.mapper2[duplicated(df.mapper2$hgnc_id) | 
                               duplicated(df.mapper2$hgnc_id, fromLast = TRUE), ] # nRow = 67,372
dup.entrezgene_id1 <- df.mapper2[duplicated(df.mapper2$entrezgene_id) | 
                                     duplicated(df.mapper2$entrezgene_id, fromLast = TRUE), ] # nRow = 79,401
dup.external_gene_name1 <- df.mapper2[duplicated(df.mapper2$external_gene_name) | 
                                          duplicated(df.mapper2$external_gene_name, fromLast = TRUE), ] # nRow = 46,952

### 시도1: 모든 gene identifier를 group
# Note: 이렇게 했을 때, 아직 여전히 duplicates이 있음.
#   - 하나의 유전자 이름에 mapping되는 여러 gene identifier를 찾는 것이 목적이니, 유전자 이름을 기준으로 다 묶기.
df.mapper3 <- df.mapper2 %>%
    group_by(ensembl_gene_id, hgnc_symbol, hgnc_id, entrezgene_id, external_gene_name,
             chromosome_name, start_position, end_position) %>%
    reframe(external_synonym=paste(external_synonym, collapse=";")) %>%
    ungroup() %>%
    as.data.frame() # nRow = 60,729

dup.ensemblid2 <- df.mapper3[duplicated(df.mapper3$ensembl_gene_id) | 
                                 duplicated(df.mapper3$ensembl_gene_id, fromLast = TRUE), ] # nRow = 273
dup.hgnc_symbol2 <- df.mapper3[duplicated(df.mapper3$hgnc_symbol) | 
                                   duplicated(df.mapper3$hgnc_symbol, fromLast = TRUE), ] # nRow = 21,935
dup.hgnc_id2 <- df.mapper3[duplicated(df.mapper3$hgnc_id) | 
                               duplicated(df.mapper3$hgnc_id, fromLast = TRUE), ] # nRow = 21,935
dup.entrezgene_id2 <- df.mapper3[duplicated(df.mapper3$entrezgene_id) | 
                                     duplicated(df.mapper3$entrezgene_id, fromLast = TRUE), ] # nRow = 35,223
dup.external_gene_name2 <- df.mapper3[duplicated(df.mapper3$external_gene_name) | 
                                          duplicated(df.mapper3$external_gene_name, fromLast = TRUE), ] # nRow = 1,529

## external_gene_name를 기준으로 다 grouping.
# Note: ensembl version을 맞춘 ensembl ID를 가지고 gene symbol을 맵핑하는 것이 안전하다.
# Reference: https://www.researchgate.net/post/How-to-deal-with-multiple-ensemble-IDs-mapping-to-one-gene-symbol-in-a-RNA-Seq-dataset
df.mapper3 <- df.mapper2 %>%
    group_by(ensembl_gene_id, chromosome_name, start_position, end_position) %>%
    reframe(
        external_gene_name=paste(unique(external_gene_name), collapse=";"),
        hgnc_symbol=paste(unique(hgnc_symbol), collapse=";"), 
        hgnc_id=paste(unique(hgnc_id), collapse=";"), 
        entrezgene_id=paste(unique(entrezgene_id), collapse=";"),
        external_synonym=paste(unique(external_synonym), collapse=";"),
        gene_biotype=paste(unique(gene_biotype), collapse=";")
    ) %>%
    ungroup() %>%
    as.data.frame() # nRow = 60,579

dup.ensemblid2 <- df.mapper3[duplicated(df.mapper3$ensembl_gene_id) | 
                                          duplicated(df.mapper3$ensembl_gene_id, fromLast = TRUE), ] # nRow = 1,256

######## 5. GRCh37 붙이기.
ensembl_hg19 <- useMart(
    biomart = "ENSEMBL_MART_ENSEMBL",
    dataset = "hsapiens_gene_ensembl",
    host = "https://grch37.ensembl.org"
)

df.hg19 <- getBM(
    attributes = c(
        "ensembl_gene_id",
        "chromosome_name",
        "start_position",
        "end_position"
    ),
    mart = ensembl_hg19
) # nRow = 63,677

# Duplicated ensembl ID 있는지 확인.
dup.ensemblid3 <- df.hg19[duplicated(df.hg19$ensembl_gene_id) | 
                                 duplicated(df.hg19$ensembl_gene_id, fromLast = TRUE), ] # nRow = 0

# autosome 및 sex chromosome만 남기기
df.hg19_filtchr <- df.hg19 %>%
    filter(chromosome_name %in% as.character(c(1:22, "X", "Y"))) # nRow =57,736

# rename column before merging
df.hg19_filtchr_rename <- df.hg19_filtchr %>%
    dplyr::rename(chromosome_name.hg19=chromosome_name, start_poisition.hg19=start_position,
           end_position.hg19=end_position)

# merge
df.mapper4 <- df.mapper3 %>%
    left_join(df.hg19_filtchr_rename, by="ensembl_gene_id")

# Missing hg19
# 이런 missing인 유전자에 대한 hg19 position 정보가 필요로 하다면, 따로 liftOver를 해주기.
na_hg19 <- df.mapper4 %>% filter(is.na(chromosome_name.hg19) | is.na(start_poisition.hg19)) # nRow = 8,491

### Save
saveRDS(df.mapper4, "gene_id_info.Ensembl_v102.rds")
