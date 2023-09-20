
### Load packages from env_R.R
list.of.packages <- c("argparse", "tidyverse", "rstudioapi", "data.table", "dplyr", "ggplot2", "ggrepel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

suppressPackageStartupMessages({
    require(argparse)
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
})

### Define parser arguments
parser <- argparse::ArgumentParser(description=":: Manhattan plot ::", 
                                formatter_class="argparse.ArgumentDefaultsHelpFormatter")

# GWAS summary statistics
parser$add_argument("--gwas", required=TRUE,
                    help="Path to the GWAS summary statistics.")
# Column names
parser$add_argument("--snp-col", dest="snp_col",  required=TRUE,
                    help="Specify SNP column name.")
parser$add_argument("--chr-col", dest="chr_col",  required=TRUE,
                    help="Specify CHR column name.")
parser$add_argument("--pos-col", dest="pos_col",  required=TRUE,
                    help="Specify POS column name.")
parser$add_argument("--p-col", dest="p_col",  required=TRUE,
                    help="Specify P column name.")
parser$add_argument("--beta-col", dest="beta_col",  required=FALSE, default="NA",
                    help="Specify BETA column name.")
parser$add_argument("--se-col", dest="se_col",  required=FALSE, default="NA",
                    help="Specify SE column name.")
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
parser$add_argument("--dpi", required=FALSE, default=300, type="integer",
                    help="Specify DPI for png image.")
parser$add_argument("--outf", required=FALSE, default="NA",
                    help="Specify the name of the output file. Default = manhattan.<Input file name>")
parser$add_argument("--outd", required=FALSE, default="NA",
                    help="Specify the output directory. Default = current working directory.")

parser$add_argument("--scale", required=FALSE, default=1, type="double",
                    help="Specify the scale of the plot. Default=1")

parser$add_argument("--width", required=FALSE, default=180, type="double",
                    help="Specify the width of the plot. Default=180 mm")
parser$add_argument("--height", required=FALSE, default=100, type="double",
                    help="Specify the height of the plot. Default=100 mm")
parser$add_argument("--units", required=FALSE, default="mm",
                    help="Specify the units for the width and height of the plot. Default = 'mm'")


### Call packages
library(tidyverse)
getCurrentFileLocation <-  function()
{
    this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
    if(length(this_file)==0){
        this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
}

env_file <- paste(getCurrentFileLocation(), "fnc_manhattan_plot.R", sep="/")

source(env_file)

###########################################################################################################
### Handle parsers
# Call-in arguments
args <- parser$parse_args()
gwas <- args$gwas
snp_col <- args$snp_col
chr_col <- args$chr_col
pos_col <- args$pos_col
p_col <- args$p_col
beta_col <- args$beta_col
se_col <- args$se_col
snps_to_annotate <- args$snps_to_annotate
color_annotate <- args$color_annotate
color1 <- args$color1
color2 <- args$color2
chr_select <- args$chr_select
img_type <- args$img_type
dpi <- args$dpi
outf <- args$outf
outd <- args$outd

scale <- args$scale

width <- args$width
height <- args$height
units <- args$units

# Read annotation file if specified
if(snps_to_annotate != "NA"){
    snpsOfInterest <- readLines(snps_to_annotate)
} else{
    snpsOfInterest <- c("NA")
}
# outf
if (outf == "NA"){
    outf <- paste0("manhattan.", basename(gwas))
}
# outd
if (outd == "NA"){
    outd <- getwd()
}


plot_manhattan(gwas=gwas, 
            snp_col=snp_col, chr_col=chr_col, pos_col=pos_col, p_col=p_col, 
            beta_col=beta_col, se_col=se_col,
            snps_to_annotate=snps_to_annotate, color_annotate=color_annotate,
            color1=color1, color2=color2, 
            chr_select=chr_select,
            img_type=img_type, dpi=dpi,
            outf=outf, outd=outd,
            width=width, height=height, units=units,
            scale=scale
            )

