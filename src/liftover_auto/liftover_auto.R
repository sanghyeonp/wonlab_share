# :: LiftOver automated ::
# - Sanghyeon Park
# - 2024.02.17
# - First in-use project: 

options(scipen = 999) # No scientific number
library(data.table)
library(dplyr)
library(argparse)

# Parse command line arguments
parser <- ArgumentParser()
parser$add_argument("--file-in", dest="file_in", type = "character", required = TRUE,
                    help="Input file")
parser$add_argument("--delim", type = "character", required = TRUE,
                    help="Input file")

parser$add_argument("--snp-col", dest = "snp_col", type = "character", required = TRUE,
                    help="Input file")
parser$add_argument("--chr-col", dest = "chr_col", type = "character", required = TRUE,
                    help="Input file")
parser$add_argument("--pos-col", dest = "pos_col", type = "character", required = TRUE,
                    help="Input file")

parser$add_argument("--build-from", dest = "build_from", type = "integer", required = TRUE,
                    help="Input file")
parser$add_argument("--build-to", dest = "build_to", type = "integer", required = TRUE,
                    help="Input file")

parser$add_argument("--do-not-save-lifted-merged", dest = "do_not_save_lifted_merged", 
                    action = "store_true", required = FALSE, default = FALSE,
                    help="Output file")
parser$add_argument("--save-mapping-file", dest = "save_mapping_file", 
                    action = "store_true", required = FALSE, default = FALSE,
                    help="Output file")

parser$add_argument("--out-pref", dest = "out_pref", type = "character", required = FALSE, default = "NA",
                    help="Output file")
parser$add_argument("--out-dir", dest = "out_dir", type = "character", required = FALSE, default = "NA",
                    help="Output file")

parser$add_argument("--thread", type = "integer", required = FALSE, default = 1, 
                    help="Output file")


args <- parser$parse_args()
file_in <- args$file_in
delim <- args$delim

snp_col <- args$snp_col
chr_col <- args$chr_col
pos_col <- args$pos_col

build_from <- args$build_from
build_to <- args$build_to

do_not_save_lifted_merged <- args$do_not_save_lifted_merged
save_mapping_file <- args$save_mapping_file

out_pref <- args$out_pref
if (out_pref == "NA") {
    out_pref <- basename(file_in)
}
out_dir <- args$out_dir
if (out_dir == "NA") {
    out_dir <- "."
}
if (substr(out_dir, nchar(out_dir), nchar(out_dir)) == "/"){
    out_dir <- substr(out_dir, 1, nchar(out_dir) - 1)
}

n_thread <- args$thread

### Pre-defined variables
map_delim <- c("tab" = "\t", "comma" = ",", "semicolon" = ";", "whitespace" = " ")


### Read in the input file and make bed file
df.bed <- fread(file_in, sep=map_delim[delim], data.table = F, nThread = n_thread, na.strings = c("")) %>%
    # If SNP is NA, then make it as CHR:POS
    dplyr::mutate(!!as.name(snp_col) := ifelse(is.na(!!as.name(snp_col)), 
                                                paste0(!!as.name(chr_col), ":", !!as.name(pos_col)), 
                                                !!as.name(snp_col)),
                !!as.name(snp_col) := ifelse(grepl(" ", !!as.name(snp_col)), gsub(" ", "", !!as.name(snp_col)), !!as.name(snp_col))) %>%
    # Remove duplicated SNPs
    dplyr::distinct(!!as.name(snp_col), .keep_all = T) %>%
    dplyr::mutate(!!as.name(chr_col) := ifelse(grepl("chr", !!as.name(chr_col)), !!as.name(chr_col), paste("chr", !!as.name(chr_col), sep="")),
                `POS-1` = as.integer(!!as.name(pos_col) - 1),
                !!as.name(pos_col) := as.integer(!!as.name(pos_col))) %>%
    dplyr::select(!!as.name(chr_col), `POS-1`, !!as.name(pos_col), !!as.name(snp_col))

file.bed <- paste0(out_pref, ".bed")
write.table(df.bed,
            file.bed,
            sep="\t", row.names = F, col.names = F, quote = F)

### Perform LiftOver
LIFTOVER_DIR <- "/data1/sanghyeon/wonlab_contribute/combined/software/liftover" 
LIFTOVER_SOFTWARE <- paste0(LIFTOVER_DIR, "/liftOver")
LIFTOVER_CHAIN_DICT <- c("36:37" = paste(LIFTOVER_DIR, "chainfile", "hg18ToHg19.over.chain.gz", sep = "/"),
                        "36:38" = paste(LIFTOVER_DIR, "chainfile", "hg18ToHg38.over.chain.gz", sep = "/"),
                        "37:38" = paste(LIFTOVER_DIR, "chainfile", "hg19ToHg38.over.chain.gz", sep = "/"),
                        "38:37" = paste(LIFTOVER_DIR, "chainfile", "hg38ToHg19.over.chain.gz", sep = "/")
                        )


file_in.bed <- "gwas_catalog.accessed_20230131.b38.liftover_in.bed"


cmd <- paste0(LIFTOVER_SOFTWARE,
            " ", file.bed,
            " ", LIFTOVER_CHAIN_DICT[paste0(build_from, ":", build_to)],
            " ", out_pref, ".lifted", 
            " ", out_pref, ".unlifted"
            )
system(cmd, wait = TRUE)

### Merge the lifted file with the original file
df.lifted <- fread(paste0(out_pref, ".lifted"), sep="\t", data.table = F, nThread = n_thread,
                header = F, col.names = c(paste0("CHR_b", build_to), "POS-1", paste0("POS_b", build_to), "SNP")) %>%
            dplyr::select(-`POS-1`) %>%
            mutate(!!as.name(paste0("CHR_b", build_to)) := gsub("chr", "", !!as.name(paste0("CHR_b", build_to))))

df <- fread(file_in, sep=map_delim[delim], data.table = F, nThread = n_thread, na.strings = c("")) %>%
    dplyr::mutate(!!as.name(snp_col) := ifelse(is.na(!!as.name(snp_col)), 
                                        paste0(!!as.name(chr_col), ":", !!as.name(pos_col)), 
                                        !!as.name(snp_col)),
                !!as.name(snp_col) := ifelse(grepl(" ", !!as.name(snp_col)), gsub(" ", "", !!as.name(snp_col)), !!as.name(snp_col)))

df <- merge(df, df.lifted, by.x = snp_col, by.y = "SNP", all.x = T)

if (!do_not_save_lifted_merged){
    write.table(df, 
                paste0(out_dir, "/", out_pref, ".lifted_merged"),
                sep=map_delim[delim], row.names = F, quote = F)
}

if (save_mapping_file){
    write.table(df %>%
                    dplyr::select(!!as.name(snp_col), !!as.name(chr_col), !!as.name(pos_col), paste0("CHR_b", build_to), paste0("POS_b", build_to)),
                paste0(out_dir, "/", out_pref, ".mapping"),
                sep="\t", row.names = F, quote = F)
}
