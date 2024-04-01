# :: Multiple test correction for MAGMA from FUMA ::
# - Sanghyeon Park
# - 2024.03.23
# - First in-use project: GWAS-by-subtraction SCZ

options(scipen=999, error=traceback) # No scientific number
library(data.table)
library(dplyr)
library(tidyr)
library(stats)
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
parser$add_argument("--dir-fuma", dest="dir_fuma", type = "character", required = FALSE, default = "NA",
                    help="Directory with FUMA results.")
parser$add_argument("--file-magma", dest="file_magma", type = "character", required = FALSE, default = "NA",
                    help="File with MAGMA results.")

parser$add_argument("--remove-gene-placeholder", dest="remove_gene_placeholder", 
                    required = FALSE, action = "store_true", default = FALSE,
                    help="Remove gene placeholder.")

parser$add_argument("--thread", dest="n_thread", type = "integer", required = FALSE, default = 1, 
                    help="Number of threads to use.")

parser$add_argument("--dir-out", dest = "dir_out", type = "character", required = FALSE, default = ".",
                    help="Output directory. Default = current working directory.")
parser$add_argument("--name-out", dest = "name_out", type = "character", required = FALSE, 
                    default = "FUMA.MAGMA.gene_based.MTC",
                    help="Output file prefix. Default = 'FUMA.gene_mapping.combined'")
parser$add_argument("--delim-out", dest = "delim_out", type = "character", required = FALSE, default = "comma",
                    help="Output file delimiter. Default = comma. Choices = comma, tab, whitespace.")

args <- parser$parse_args()
dir_fuma <- remove_trailing_slash(args$dir_fuma)
file_magma <- args$file_magma
if (dir_fuma == "NA" & file_magma == "NA") {
    stop("Either --dir-fuma or --file-magma should be provided.")
}
if (dir_fuma != "NA"){
    file_magma <- paste0(dir_fuma, "/magma.genes.out")
}

remove_gene_placeholder <- args$remove_gene_placeholder

dir_out <- remove_trailing_slash(args$dir_out)
name_out <- args$name_out
delim_out <- args$delim_out
n_thread <- args$n_thread

###
df <- fread(file_magma, sep="\t", data.table=F, nThread=n_thread)
df <- df %>%
    mutate(P_FDR = stats::p.adjust(P, method="fdr"),
        P_Bonferroni = stats::p.adjust(P, method="bonferroni")) %>%
    arrange(P)

if (!remove_gene_placeholder){
    df <- df %>%
        mutate(SYMBOL = paste0("`", SYMBOL))
}

print(paste0("N genes FDR significant: ", sum(df$P_FDR < 0.05)))
print(paste0("N genes Bonferroni significant: ", sum(df$P_Bonferroni < 0.05)))

write.table(df,
            paste0(dir_out, "/", name_out, delim_extension_map[delim_out]),
            sep=delim_map[delim_out], row.names = F, quote = F)
