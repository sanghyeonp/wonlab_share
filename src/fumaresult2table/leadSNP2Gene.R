options(scipen=999, error=traceback) # No scientific number
library(data.table)
library(dplyr)
library(tidyr)
library(argparse)

### Pre-defined function
remove_trailing_slash <- function(input_string) {
    # Replace a forward slash at the end of the string with an empty string
    return(sub("/$", "", input_string))
}

delim_map <- c("comma" = ",", "tab" = "\t", "whitespace" = " ")
delim_extension_map <- c("comma" = ".csv", "tab" = ".tsv", "whitespace" = ".txt")

### Parse command line arguments
parser <- ArgumentParser()
parser$add_argument("--dir-fuma", dest="dir_fuma", type = "character", required = TRUE,
                    help="Directory with FUMA results.")

parser$add_argument("--posMap-max-distance", dest = "posMap_max_distance", type = "integer", required = FALSE, default = 10,
                    help="Maximum distance used for position mapping.")
parser$add_argument("--thread", dest="n_thread", type = "integer", required = FALSE, default = 1, 
                    help="Number of threads to use.")

parser$add_argument("--remove-gene-placeholder", dest="remove_gene_placeholder", 
                    required = FALSE, action = "store_true", default = FALSE,
                    help="Remove gene placeholder.")

parser$add_argument("--dir-out", dest = "dir_out", type = "character", required = FALSE, default = ".",
                    help="Output directory. Default = current working directory.")
parser$add_argument("--name-out", dest = "name_out", type = "character", required = FALSE, 
                    default = "FUMA.gene_mapping.combined",
                    help="Output file prefix. Default = 'FUMA.gene_mapping.combined'")
parser$add_argument("--delim-out", dest = "delim_out", type = "character", required = FALSE, default = "comma",
                    help="Output file delimiter. Default = comma. Choices = comma, tab, whitespace.")

args <- parser$parse_args()
dir_fuma <- remove_trailing_slash(args$dir_fuma)

remove_gene_placeholder <- args$remove_gene_placeholder

posMap_max_distance <- args$posMap_max_distance
n_thread <- args$n_thread

dir_out <- remove_trailing_slash(args$dir_out)
name_out <- args$name_out
delim_out <- args$delim_out

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


file_list <- list.files(dir_fuma,
                    full.names = TRUE)
potential_file_list <- c("genes.txt", "snps.txt",
                        "eqtl.txt",
                        "ci.txt")
potential_file_exists <- potential_file_list %in% file_list


for (filename in potential_file_list){
    print(paste0(filename, " : ", ifelse(file.exists(paste0(dir_fuma, "/", filename)), "exists", "does not exist")))
}

### Start with positional mapping since it is always performed by FUMA
df.pos <- table.FUMA.positional_mapping(genes_filepath = paste0(dir_fuma, "/genes.txt"),
                                    snps_filepath = paste0(dir_fuma, "/snps.txt"),
                                    posMap_max_distance = posMap_max_distance,
                                    n_thread = n_thread)

df_merged <- df.pos

### eQTL mapping
df.eqtl <- data.frame()
if (file.exists(paste0(dir_fuma, "/genes.txt")) &
    file.exists(paste0(dir_fuma, "/eqtl.txt")) &
    file.exists(paste0(dir_fuma, "/snps.txt"))){
    df.eqtl <- table.FUMA.eQTL_mapping(genes_filepath = paste0(dir_fuma, "/genes.txt"),
                                       eqtl_filepath = paste0(dir_fuma, "/eqtl.txt"),
                                       snps_filepath = paste0(dir_fuma, "/snps.txt"),
                                       n_thread = n_thread)
    
    df_merged <- merge(df_merged, df.eqtl, by = c("Ensembl ID", "Gene", "Gene type", "Chr", "Start", "End"), all = T)
}


### Chromatin interaction mapping
df.ci <- data.frame()
if (file.exists(paste0(dir_fuma, "/ci.txt"))){
    df.ci <- table.FUMA.chromatin_interaction_mapping(ci_filepath = paste0(dir_fuma, "/ci.txt"),
                                                      n_thread = n_thread)
    
    df_merged <- merge(df_merged, df.ci, by = "Ensembl ID", all = T) %>%
        replace_na(list(ciMap = "No"))
}

### Add ` to Gene name
if (!remove_gene_placeholder){
    df_merged <- df_merged %>%
        mutate(Gene = paste0("`", Gene))
}


### Save
write.table(df_merged, 
            paste0(dir_out, "/", name_out, ".", delim_extension_map[delim_out]),
            sep=delim_map[delim_out], row.names = F, quote = T)