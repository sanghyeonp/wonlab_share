library(data.table)
library(dplyr)
library(tidyr)


dir_out <- "/data1/sanghyeon/Projects/MetabolicSyndrome/MetS_2022_08/nat_genet.revision_1/reviewer_2/R2Q1/FUMA/MetS_GWAS_COJO/out"
posMap_max_distance <- 10
n_thread <- 20

### Function
table.FUMA.positional_mapping <- function(genes_filepath,
                    snps_filepath,
                    posMap_max_distance = 10,
                    n_thread = 1){
    # Require two files
    #   - genes.txt
    #   - snps.txt
    # posMap_max_distance (in Kb): maximum distance between SNP and gene start and end position for positional mapping

    df.snps <- fread(snps_filepath, sep="\t", data.table = F, nThread = n_thread) %>%
        dplyr::select(rsID, chr, pos, posMapFilt)

    df.posMap <- fread(genes_filepath, sep="\t", data.table = F, nThread = n_thread) %>%
        dplyr::select(ensg, symbol, type, chr, start, end, 
                    posMapSNPs, posMapMaxCADD) %>%
        dplyr::rename(`Ensembl ID` = ensg, Gene = symbol, `Gene type` = type,
                    Chr = chr, Start = start, End = end,
                    `N SNP posMap` = posMapSNPs, `Max CADD posMap` = posMapMaxCADD) %>%
        dplyr::mutate(posMap = ifelse(`N SNP posMap` > 0, "Yes", "No"))
                    
    df.posMap_SNP_mapped <- as.data.frame(df.posMap %>%
                                            inner_join(df.snps, by = c("Chr" = "chr"), relationship = "many-to-many") %>%
                                            filter(pos >= Start - posMap_max_distance * 1000 & 
                                                    pos <= End + posMap_max_distance * 1000 &
                                                    posMapFilt == 1) %>%
                                            group_by(`Ensembl ID`, Chr, Start, End) %>%
                                            summarise(`SNP posMap` = paste(rsID, collapse = ";")) %>%
                                            ungroup())
    
    df.posMap <- df.posMap %>%
                    left_join(df.posMap_SNP_mapped, by = c("Ensembl ID", "Chr", "Start", "End"))
    
    return (df.posMap)
}

table.FUMA.eQTL_mapping <- function(genes_filepath,
                    eqtl_filepath,
                    snps_filepath,
                    n_thread = 1){
    # Require two files
    #   - genes.txt
    #   - eqtl.txt
    #   - snps.txt
    
    df.snps <- fread(snps_filepath, sep="\t", data.table = F, nThread = n_thread) %>%
        dplyr::select(uniqID, rsID)

    df.eqtl <- as.data.frame(fread(eqtl_filepath, sep="\t", data.table = F, nThread = n_thread) %>%
            dplyr::select(uniqID, symbol, eqtlMapFilt) %>%
            filter(eqtlMapFilt == 1) %>%
            left_join(df.snps, by = "uniqID") %>%
            dplyr::select(-uniqID))

    df.eqtlMap <- fread(genes_filepath, sep="\t", data.table = F, nThread = n_thread) %>%
        dplyr::select(ensg, symbol, type, chr, start, end,
                    eqtlMapSNPs, eqtlMapminP, eqtlMapminQ, eqtlMapts) %>%
        dplyr::rename(`Ensembl ID` = ensg, Gene = symbol, `Gene type` = type,
                    Chr = chr, Start = start, End = end,
                    `N SNP eqtlMap` = eqtlMapSNPs, `Min P eqtlMap` = eqtlMapminP, 
                    `Tissue eqtlMap` = eqtlMapts) %>%
        dplyr::mutate(eqtlMap = ifelse(`N SNP eqtlMap` > 0, "Yes", "No"))

    df.eqtlMap_SNP_mapped <- as.data.frame(df.eqtlMap %>%
                                            inner_join(df.eqtl, by = c("Gene" = "symbol"), relationship = "many-to-many") %>%
                                            group_by(Gene) %>%
                                            summarise(`SNP eQTL` = paste(rsID, collapse = ";")) %>%
                                            ungroup())
    df.eqtlMap <- df.eqtlMap %>%
        left_join(df.eqtlMap_SNP_mapped, by = "Gene")
    
    return (df.eqtlMap)
}

table.FUMA.chromatin_interaction_mapping <- function(ci_filepath,
                                                    n_thread = 1){
    # Require one file
    #   - ci.txt
    # Lead SNP이 속한 region 1과 interaction이 있는 region 2에 overlap되는 유전자들이 mapping됨.
    # 그래서 하나의 lead SNP이 속한 region 1과 interaction이 있는 region 2와 overlap되는 유전자들이 2개 이상일 수도 있고,
    # 또, 여러 lead SNP이 같은 region 1에 속할 수도 있고,
    # 또, 다른 region1이 같은 region2와 interaction이 있을 수 있음.
    # 그래도 mapping된 unique 유전자에 어떤 lead SNP들, 그리고 tissue/cell에서 mapping 되었는지 정리함.

    df.ciMap <- as.data.frame(fread(ci_filepath, sep="\t", data.table = F, nThread = n_thread) %>%
        dplyr::select(SNPs, genes, `tissue/cell`, ciMapFilt) %>%
        filter(ciMapFilt == 1) %>%
        tidyr::separate_longer_delim(genes, delim = ":") %>%
        distinct(SNPs, genes, .keep_all = T) %>%
        group_by(genes) %>%
        summarise(`SNP ci` = paste0(SNPs, collapse = ";"),
                  `Tissue ci` = paste0(`tissue/cell`, collapse = ";")) %>%
        mutate(ciMap = "Yes") %>%
        dplyr::rename(`Ensembl ID` = genes))
    
    return (df.ciMap)
}
###


file_list <- list.files(dir_out,
                    full.names = TRUE)
potential_file_list <- c("genes.txt", "snps.txt",
                        "eqtl.txt",
                        "ci.txt")
potential_file_exists <- potential_file_list %in% file_list


for (filename in potential_file_list){
    print(paste0(filename, " : ", ifelse(file.exists(paste0(dir_out, "/", filename)), "exists", "does not exist")))
}

### Start with positional mapping since it is always performed by FUMA
df.pos <- table.FUMA.positional_mapping(genes_filepath = paste0(dir_out, "/genes.txt"),
                                    snps_filepath = paste0(dir_out, "/snps.txt"),
                                    posMap_max_distance = posMap_max_distance,
                                    n_thread = n_thread)

df_merged <- df.pos

### eQTL mapping
df.eqtl <- data.frame()
if (file.exists(paste0(dir_out, "/genes.txt")) &
    file.exists(paste0(dir_out, "/eqtl.txt")) &
    file.exists(paste0(dir_out, "/snps.txt"))){
    df.eqtl <- table.FUMA.eQTL_mapping(genes_filepath = paste0(dir_out, "/genes.txt"),
                                       eqtl_filepath = paste0(dir_out, "/eqtl.txt"),
                                       snps_filepath = paste0(dir_out, "/snps.txt"),
                                       n_thread = n_thread)
    
    df_merged <- merge(df_merged, df.eqtl, by = c("Ensembl ID", "Gene", "Gene type", "Chr", "Start", "End"), all = T)
}


### Chromatin interaction mapping
df.ci <- data.frame()
if (file.exists(paste0(dir_out, "/ci.txt"))){
    df.ci <- table.FUMA.chromatin_interaction_mapping(ci_filepath = paste0(dir_out, "/ci.txt"),
                                                      n_thread = n_thread)
    
    df_merged <- merge(df_merged, df.ci, by = "Ensembl ID", all = T) %>%
        replace_na(list(ciMap = "No"))
}

### Add ` to Gene name
df_merged <- df_merged %>%
    mutate(Gene = paste0("`", Gene))

### Save
write.table(df_merged, 
            "FUMA.gene_mapping.combined.csv",
            sep=",", row.names = F, quote = T)