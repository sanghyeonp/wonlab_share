# :: SNP annotation using ANNOVAR ::
# - Sanghyeon Park
# - 2023.12.14
# - First in-use project: 

source("/data1/sanghyeon/wonlab_contribute/combined/src/annovar/make_annovar_input.R")
source("/data1/sanghyeon/wonlab_contribute/combined/src/annovar/run_annovar.R")
source("/data1/sanghyeon/wonlab_contribute/combined/src/annovar/merge_annovar_out.R")

library(data.table)
library(dplyr)
library(argparse)


### Define parser arguments
parser <- argparse::ArgumentParser(description=":: Run ANNOVAR to annotate rsID ::", 
                                formatter_class="argparse.ArgumentDefaultsHelpFormatter"
                                )

parser$add_argument("--gwas", required=TRUE, 
                    help="Path to GWAS summary statistics.")
parser$add_argument("--delim-in", dest="delim_in", required=FALSE, default="tab",
                help="Delimiter used in the input file. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.")

parser$add_argument("--chr-col", dest="chr_col", required=FALSE, default="CHR",
                help="Name of the chromosome column in the input file. Default=CHR.")
parser$add_argument("--pos-col", dest="pos_col", required=FALSE, default="POS",
                help="Name of the base position column in the input file. Default=POS.")
parser$add_argument("--ref-col", dest="ref_col", required=FALSE, default="REF",
                help="Name of the reference allele column in the input file. Default=REF.")
parser$add_argument("--alt-col", dest="alt_col", required=FALSE, default="ALT",
                help="Name of the alternative allele column in the input file. Default=ALT.")

parser$add_argument("--genome-build", dest="genome_build", required=FALSE, default=19, type="integer",
                help="Genome build of GWAS summary statistics. Default=19. Options = 19, 38.")
parser$add_argument("--dbgap-build", dest="dbgap_build", required=FALSE, default=150, type="integer",
                help="dbGAP build for rsID annotation. Default=150. Options = 142, 144, 147, 150.")

parser$add_argument("--nthread", required=FALSE, default=1, type="integer",
                help="Number of threads to use. Default=1.")

parser$add_argument("--rm-flip-col", dest="rm_flip_col", action="store_true",
                help="Specify to remove the 'flipped' column from the output. Default=FALSE.")
parser$add_argument("--no-gwas-annotation", dest="no_gwas_annotation", action="store_true",
                help="Specify to not save annotated GWAS summary statistics. Default=FALSE.")
parser$add_argument("--save-mapping-file", dest="save_mapping_file", action="store_true",
                help="Specify to save the rsID mapping file. Default=FALSE.")

parser$add_argument("--outf", required=FALSE, default="annovar_temp",
                help="Name of the output file.")

args <- parser$parse_args()

###
delim_map <- c("tab" = "\t", "comma" = ",", "whitespace" = " ")

### Make ANNOVAR input
output_list <- make_annovar_input(
    gwas = args$gwas, delim_in = delim_map[args$delim_in],
    chr_col = args$chr_col, pos_col = args$pos_col, ref_col = args$ref_col, alt_col = args$alt_col,
    outf = args$outf,
    nthread = args$nthread
)

filename_annovin <- output_list[[1]]
filename_annovin_flip <- output_list[[2]]
df <- output_list[[3]]

### Run ANNOVAR
filename_annovin_annot <- run_annovar(
    annov_in = filename_annovin,
    genome_build = args$genome_build, dbsnp_build = args$dbgap_build,
    nthread = args$nthread
)

filename_annovin_flip_annot <- run_annovar(
    annov_in = filename_annovin_flip,
    genome_build = args$genome_build, dbsnp_build = args$dbgap_build,
    nthread = args$nthread
)

### Merge ANNOVAR output
filename_annovin_annot_merged <- merge_annovar_out(
    filename_annovin_annot = filename_annovin_annot,
    filename_annovin_flip_annot = filename_annovin_flip_annot,
    nthread = args$nthread
)

### Remove temporary files
system(paste0("rm ", 
            filename_annovin, " ",
            filename_annovin_flip, " ",
            filename_annovin, ".hg", args$genome_build, "_avsnp", args$dbgap_build, "_dropped", " ",
            filename_annovin, ".hg", args$genome_build, "_avsnp", args$dbgap_build, "_filtered", " ",
            filename_annovin, ".invalid_input", " ",
            filename_annovin, ".log", " ",
            filename_annovin_flip, ".hg", args$genome_build, "_avsnp", args$dbgap_build, "_dropped", " ",
            filename_annovin_flip, ".hg", args$genome_build, "_avsnp", args$dbgap_build, "_filtered", " ",
            filename_annovin_flip, ".invalid_input", " ",
            filename_annovin_flip, ".log"
            ), wait = TRUE)


### Read annotated ANNOVAR output
df_annot <- as.data.frame(fread(filename_annovin_annot_merged, 
                                header = FALSE,
                                sep = "\t",
                                nThread = args$nthread,
                                showProgress = FALSE,
                                col.names = c("db",
                                    "rsid",
                                    "chr", 
                                    "bp", 
                                    "bp2", 
                                    "ref", 
                                    "alt",
                                    "flipped"
                                    )))

df_annot <- df_annot %>%
    mutate(new_ref = ifelse(flipped, alt, ref),
            new_alt = ifelse(flipped, ref, alt),
            variant_id = paste0(chr, ":", bp, ":", new_ref, ":", new_alt)) %>%
    select(variant_id, rsid, flipped)

if (args$rm_flip_col){
    df_annot <- df_annot %>% 
        select(-flipped)
}

if (!args$no_gwas_annotation){
    df <- df %>%
        mutate(variant_id = paste0(!!as.name(args$chr_col), ":", 
                                    !!as.name(args$pos_col), ":", 
                                    !!as.name(args$ref_col), ":", 
                                    !!as.name(args$alt_col)))
    
    df <- merge(df, df_annot, by = "variant_id", all.x = TRUE)

    df <- df %>% 
        mutate(rsid = ifelse(is.na(rsid), ".", rsid)) %>%
        select(-variant_id)

    if ("flipped" %in% colnames(df)){
        df <- df %>% 
            mutate(flipped = ifelse(is.na(flipped), ".", flipped))
    }

    write.table(df,
                paste0(args$outf, ".rsid_annotated"),
                sep = delim_map[args$delim_in],
                row.names = FALSE,
                col.names = TRUE,
                quote = FALSE)
}

if (args$save_mapping_file){
    write.table(df_annot,
                paste0(args$outf, ".rsid_mapper"),
                sep = "\t",
                row.names = FALSE,
                col.names = TRUE,
                quote = FALSE)
}