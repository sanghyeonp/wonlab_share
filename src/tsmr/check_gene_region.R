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

fnc_file1 <- paste(dirname(getCurrentFileLocation()), "gwas_plot", "fnc_manhattan_plot.R", sep="/")
fnc_file2 <- paste(dirname(getCurrentFileLocation()), "tsmr", "step_wise_fnc.R", sep="/")
source(fnc_file1)
source(fnc_file2)


library(data.table)

###############################
### Define parser arguments ###
###############################
parser <- argparse::ArgumentParser(description=":: Manhattan plot in the specified gene region ::", 
                                formatter_class="argparse.ArgumentDefaultsHelpFormatter")

parser$add_argument("--gwas", required=TRUE,
                    help="Path to the exposure GWAS summary statistics.")
parser$add_argument("--delim", required=TRUE,
                    help="Delimiter for the exposure GWAS. Options = ['tab', 'whitespace', 'comma'].")
parser$add_argument("--snp-col", dest="snp_col", required=TRUE,
                    help="Column name for SNP.")
parser$add_argument("--chr-col", dest="chr_col", required=TRUE,
                    help="Column name for the chromosome.")
parser$add_argument("--pos-col", dest="pos_col", required=TRUE,
                    help="Column name for the base position.")
parser$add_argument("--beta-col", dest="beta_col", required=TRUE,
                    help="Column name for the effect estimate.")
parser$add_argument("--se-col", dest="se_col", required=TRUE,
                    help="Column name for the SE.")
parser$add_argument("--p-col", dest="p_col", required=TRUE,
                    help="Column name for the p-value.")

## Filter by gene region
parser$add_argument("--gene-chr", dest="gene_chr", required=TRUE, type='integer',
                    help="Chromosome number of the specified gene.")
parser$add_argument("--gene-start", dest="gene_start", required=TRUE, type='integer',
                    help="Specify the gene start position.")
parser$add_argument("--gene-end", dest="gene_end", required=TRUE, type='integer',
                    help="Specify the gene end position.")
parser$add_argument("--gene-cis-window", dest="gene_cis_window", required=TRUE, type='integer',
                    help="Specify the cis-window around the gene body in kb. For example, `--gene-cis-window 500` means 500 kb around the gene body.")

### Manhattan plot options
# Annotation
parser$add_argument("--snps-to-annotate", dest="snps_to_annotate",  required=FALSE, default="NA",
                    help="Specify path containing list of SNPs to annotate.")
parser$add_argument("--color-annotate", dest="color_annotate", required=FALSE, default="darkorange1",
                    help="Specify the color for SNP annotation.")
# Color specification
parser$add_argument("--color1", required=FALSE, default="grey50",
                    help="Specify the first color.")
parser$add_argument("--color2", required=FALSE, default="grey",
                    help="Specify the second color.")
# Others
parser$add_argument("--chr-select", dest="chr_select", nargs = "*", default="NA",
                    help="Specify chromosomes to plot.")
# Save output
parser$add_argument("--img-type", dest="img_type", required=FALSE, default="png",
                    help="Specify the type of image extension. Options = ['png', 'pdf']. Default = 'png'")
parser$add_argument("--dpi", required=FALSE, default=300, type='integer',
                    help="Specify DPI for png image.")
parser$add_argument("--width", required=FALSE, default=180, type='double',
                    help="Specify the width of the plot. Default=180 mm")
parser$add_argument("--height", required=FALSE, default=100, type='double',
                    help="Specify the height of the plot. Default=100 mm")
parser$add_argument("--units", required=FALSE, default="mm",
                    help="Specify the units for the width and height of the plot. Default = 'mm'")

### Save option
parser$add_argument("--outf", required=TRUE, 
                    help="Specify the name of the output files.")
parser$add_argument("--outd", required=FALSE, default="NA",
                    help="Specify the output directory path. Default = `current working directory`.")
    

####
args <- parser$parse_args()

gwas <- args$gwas
delim <- args$delim
snp_col <- args$snp_col
chr_col <- args$chr_col
pos_col <- args$pos_col
beta_col <- args$beta_col
se_col <- args$se_col
p_col <- args$p_col

gene_chr <- args$gene_chr
gene_start <- args$gene_start
gene_end <- args$gene_end
gene_cis_window <- args$gene_cis_window

snps_to_annotate <- args$snps_to_annotate
color_annotate <- args$color_annotate
color1 <- args$color1
color2 <- args$color2
chr_select <- args$chr_select
img_type <- args$img_type
dpi <- args$dpi
width <- args$width
height <- args$height
units <- args$units

outf <- args$outf
outd <- args$outd

####

delim_list <- list(tab = "\t", whitespace = " ", comma = ",")
delim <- delim_list[[delim]]

if(snps_to_annotate != "NA"){
    snpsOfInterest <- readLines(snps_to_annotate)
} else{
    snpsOfInterest <- c("NA")
}

### Gene region extraction
extract_gene_region <- function(gwas, delim,
                                snp_col, chr_col, pos_col, beta_col, se_col, p_col,
                                gene_chr, gene_start, gene_end, gene_cis_window,
                                outf, outd='NA', save_table=FALSE){
    df <- fread(gwas, 
                sep=delim)

    df_filter <- df %>% 
            select(all_of(c(snp_col, chr_col, pos_col, beta_col, se_col, p_col))) %>%
            rename(SNP := all_of(snp_col),
                CHR := all_of(chr_col),
                POS := all_of(pos_col),
                BETA := all_of(beta_col),
                SE := all_of(se_col),
                P := all_of(p_col)
                ) %>%
            mutate(P = as.numeric(P)) %>%
            filter((CHR == gene_chr) &
                (POS > gene_start - (gene_cis_window * 1000)) &
                (POS <= gene_end + (gene_cis_window * 1000))
                )
    
    if (outd == 'NA'){
        outd <- getwd()
    }

    if (save_table){
        write.table(df_filter,
                    paste0(outd, "/", outf, ".tsv"),
                    sep="\t", row.names=FALSE, quote=FALSE
                    )
    }
    
    return (paste0(outd, "/", outf, ".tsv"))
}

gwas_gene_region_file <- extract_gene_region(gwas=gwas, delim=delim,
                        snp_col=snp_col, chr_col=chr_col, pos_col=pos_col, 
                        beta_col=beta_col, se_col=se_col, p_col=p_col,
                        gene_chr=gene_chr, gene_start=gene_start, gene_end=gene_end, gene_cis_window=gene_cis_window,
                        outf=outf, outd=outd, save_table=TRUE
                        )

### Plot Manhattan
plot_manhattan(gwas=gwas_gene_region_file, 
            snp_col='SNP', chr_col='CHR', pos_col='POS', p_col='P', 
            beta_col='BETA', se_col='SE',
            snps_to_annotate=snps_to_annotate, color_annotate=color_annotate,
            color1=color1, color2=color2, 
            chr_select=chr_select,
            img_type=img_type, dpi=dpi,
            outf=paste0("manhattan.", outf), outd=outd,
            width=width, height=height, units=units
            )
