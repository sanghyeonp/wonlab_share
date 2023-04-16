
### Load packages from env_R.R
list.of.packages <- c("argparse", "tidyverse", "rstudioapi")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

require(argparse)

### Common files
variant_file <- "/data1/sanghyeon/wonlab_contribute/combined/data_common/variants.tsv.bgz"


### Define parser arguments
parser <- argparse::ArgumentParser(description=":: Reformat GWAS summary statistics ::", 
                                formatter_class="argparse.ArgumentDefaultsHelpFormatter")

parser$add_argument("--gwas", required=TRUE,
                    help="Path to the GWAS summary statistics.")
parser$add_argument("--retain-col", dest="retain_col", nargs = "*", default="NA",
                    help="Specify columns to retain. Options are [RSID, CHR, POS, ALT, REF, EAF, MAF, BETA, SE, PVAL, N]. These column names are common column names. If you want to retain other columns in the GWAS summary statistics, specify their original name.")
parser$add_argument("--rename-col", dest="rename_col", nargs = "*", default="NA",
                    help="Specify new column names specified in `--retain-col` with the same order.")

parser$add_argument("--maf-filter", dest="maf_filter", required=FALSE, default=0,
                    help="Specify MAF filter threshold. Default=0.")

parser$add_argument("--keep-no-rsid", dest="keep_no_rsid", action="store_true", 
                    help="Specify to retain the SNPs without mapped rsID. Default=FALSE.")
parser$add_argument("--keep-duplicates", dest="keep_duplicates", action="store_true", 
                    help="Specify to retain duplicated SNPs. Default=FALSE.")
parser$add_argument("--keep-na", dest="keep_na", action="store_true", 
                    help="Specify to retain rows with NA. Default=FALSE.")

parser$add_argument("--outf", required=FALSE, default="NA",
                    help="Specify the name of the output file. Default = reformat.<Input file name>")
parser$add_argument("--delim_out", required=FALSE, default="tab",
                    help="Specify the delimiter for the output file. Options = [tab, comma, whitespace]. Default = tab.")
parser$add_argument("--outd", required=FALSE, default="NA",
                    help="Specify the output directory. Default = current working directory.")

### Call packages
library(tidyverse)
getCurrentFileLocation <-  function()
{
    this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
    if (length(this_file)==0){
        this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
}

env_file <- paste(getCurrentFileLocation(), "env_R.R", sep="/")
source(env_file)

### Get parser arguments
args <- parser$parse_args()
gwas <- args$gwas
retain_col <- args$retain_col

maf_filter <- args$maf_filter
keep_no_rsid <- args$keep_no_rsid
keep_duplicates <- args$keep_duplicates
keep_na <- args$keep_na
rename_col <- args$rename_col

outf <- args$outf
delim_out <- args$delim_out
outd <- args$outd
if (outf == "NA"){
    outf <- paste0("reformat.", basename(gwas))
}
if (outd == "NA"){
    outd <- getwd()
}
out_path <- paste(outd, outf, sep="/")
################################################
log_file <- file(paste0(out_path, ".log"), open = "wt")
sink(log_file, type = "output")

if (retain_col[1] == "NA"){
    stop("Specify the columns to retain using `--retain-col`.")
}
################################################

### Read GWAS summary statistics
cat("\n::Run:: Reading input GWAS summary statistics")
df <- suppressMessages(suppressWarnings(readr::read_table(gwas, col_names = TRUE)))

### Read variant to rsID mapping file
cat("\n::Run:: Reading variant file")
df_variant <- suppressMessages(suppressWarnings(readr::read_table(variant_file, col_names = TRUE)))
df_variant <- subset(df_variant, select=c('variant', 'rsid'))
colnames(df_variant) <- c("variant", "RSID")

### Make column name renaming vector
col_map_dict <- setNames(c('variant', 'minor_allele', 'MAF', 'low_confidence_variant', 'N', 'AC', 'ytx', 'BETA', 'SE', 'tstat', 'PVAL'), 
                        c('variant', 'minor_allele', 'minor_AF', 'low_confidence_variant', 'n_complete_samples', 'AC', 'ytx', 'beta', 'se', 'tstat', 'pval'))

### Rename the columns
names(df) <- col_map_dict

### Obtain CHR, POS, ALT, REF, EAF
cat("\n::Run:: Extract CHR, POS, REF, ALT and EAF")
df <- df %>%
    separate(variant, c("CHR", "POS", "REF", "ALT"), sep = ":", convert = TRUE) %>%
    mutate(variant = paste0(CHR, ":", POS, ":", REF, ":", ALT)) %>%
    mutate(EAF = if_else(ALT == minor_allele, MAF, 1 - MAF)) 

### Map rsID
cat("\n::Run:: Mapping rsID")
df <- merge(df, df_variant, by='variant')
nsnp <- nrow(df)
cat(paste0("\n\tNumber of SNPs before: ", nsnp))
df2 <- df[grepl("^rs", df$RSID), ]
cat(paste0("\n\tNumber of SNPs without rsID: ", nsnp - nrow(df2)))
if (keep_no_rsid == FALSE){
    df <- df2
    rm(df2)
    cat(paste0("\n\tNumber of SNPs after filtering unmapped rsID: ", nrow(df)))
}

### MAF filter
if (maf_filter != 0){
    cat(paste0("\n::Run:: MAF filter with threshold: ", maf_filter))
    cat(paste0("\n\tNumber of SNPs before: ", nrow(df)))
    df <- df[df$MAF > maf_filter, ]
    cat(paste0("\n\tNumber of SNPs after: ", nrow(df)))
} else{
    cat(paste0("::Run:: No MAF filtering"))
}

### Remove duplicates
cat("\n::Run:: Handle duplicated SNP")
df2 <- df %>%
    distinct(RSID, .keep_all = TRUE)
cat(paste0("\n\tNumber of SNPs before: ", nrow(df)))
cat(paste0("\n\tNumber of duplicated SNPs: ", nrow(df) - nrow(df2)))
if (keep_duplicates == FALSE){
    cat("\n\tRemoving duplicates...")
    df <- df2
    rm(df2)
    cat(paste0("\n\tNumber of SNPs after: ", nrow(df)))
} else{
    cat("\n\tKeeping duplicates...")
}

### Remove rows with NA
cat("\n::Run:: Handle rows with NA")
df2 <- df[complete.cases(df), ]
cat(paste0("\n\tNumber of rows: ", nrow(df)))
cat(paste0("\n\tNumber of rows with NA: ", nrow(df) - nrow(df2)))
if (keep_na == FALSE){
    cat("\n\tRemoving rows with NA...")
    df <- df2
    rm(df2)
    cat(paste0("\n\tNumber of rows after: ", nrow(df)))
} else{
    cat("\n\tKeeping rows with NA...")
}


### Subset specified columns
cat("\n::Run:: Subset specified columns")
cat(paste0("\n\tColumns specified: ", paste(retain_col, collapse=", ")))

df <- subset(df, select=as.vector(retain_col))


### Make column name renaming vector
cat("\n::Run:: Rename the columns")
if (rename_col[1] == "NA"){
    cat("\n\t`--rename-col` not specified. No renaming done.")
} else{
    cat(paste0("\n\tNew column names: ", paste(rename_col, collapse=", ")))
    col_new_name <- setNames(rename_col, retain_col)
    names(df) <- col_new_name
}

### Sort by CHR and POS
# 여기서 df_sorted <- df[order(df$col1, df$col2), ]
# 이런식으로 sorting 하기

### Save the output
cat("\n::Run:: Save output GWAS summary statistics")
cat(paste0("\n\tSaved at: ", out_path))

delim_map <- list(tab = "\t", comma = ",", whitespace = " ")

write.table(df, 
            out_path, 
            sep = delim_map[[delim_out]],
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

sink(type = "output")
close(log_file)
