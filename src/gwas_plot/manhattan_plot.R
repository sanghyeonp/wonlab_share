
### Load packages from env_R.R
list.of.packages <- c("argparse", "tidyverse", "rstudioapi")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

require(argparse)

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
parser$add_argument("--dpi", required=FALSE, default=300,
                    help="Specify DPI for png image.")
parser$add_argument("--outf", required=FALSE, default="NA",
                    help="Specify the name of the output file. Default = manhattan.<Input file name>")
parser$add_argument("--outd", required=FALSE, default="NA",
                    help="Specify the output directory. Default = current working directory.")
parser$add_argument("--width", required=FALSE, default=180,
                    help="Specify the width of the plot. Default=180 mm")
parser$add_argument("--height", required=FALSE, default=100,
                    help="Specify the height of the plot. Default=180 mm")
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

env_file <- paste(getCurrentFileLocation(), "env_R.R", sep="/")

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
# oudd
if (outd == "NA"){
    outd <- getwd()
}
# dpi, width, height
dpi <- as.numeric(dpi)
width <- as.numeric(width)
height <- as.numeric(height)


### Read input file
df <- fread(gwas)

### Rename columns
if((beta_col != "NA") & (se_col != "NA")){
    df <- df %>%    
        select(all_of(c(snp_col, chr_col, pos_col, p_col, beta_col, se_col))) %>%
        rename("SNP" := snp_col,
                "CHR" := chr_col,
                "POS" := pos_col,
                "PVAL" := p_col,
                "BETA" := beta_col,
                "SE" := se_col
        )
} else{
    df <- df %>% 
        select(all_of(c(snp_col, chr_col, pos_col, p_col))) %>%
        rename("SNP" := snp_col,
                "CHR" := chr_col,
                "POS" := pos_col,
                "PVAL" := p_col
        )
}

### Filter chromosomes
if (chr_select[1] != "NA"){
    chr_select <- as.numeric(chr_select)
    df <- df %>% 
            filter(CHR %in% chr_select)
}

n_chr <- length(unique(df$CHR))

#####################
# # 여기 P-value를 새로 계산하고, given p-value랑 차이가 많이 나는지 확인하는 코드 작성하기.
# # 차이가 많이 안 난다면, 새로 계산한 P-value 이용. 
# if(beta_col != "NA" & se_col != "NA"){
#   df <- df %>% 
#     rename("BETA" := beta_col,
#            "SE" := se_col
#     )
#   df$PVAL_cal <- 2 * pnorm(-abs(df$BETA / df$SE))
#   
# }


#####################

### Select columns
df <- df %>%
    select(all_of(c("SNP", "CHR", "POS", "PVAL")))

### Check rows with NA
df_na <- df[!complete.cases(df), ]

if(nrow(df_na) != 0){
    cat(paste0("Number of removed SNPs with missing value: ", nrow(df_na)))

    write.table(df_na, "manhattan_SNPs_with_missing_value.txt",
                sep="\t", row.names=FALSE, col.names=TRUE
    )
}

### Remove rows with NA
df <- df[complete.cases(df), ]

### Prepare the dataset
df_plot <- df %>% 
    # Compute chromosome size
    group_by(CHR) %>%
    summarise(chr_len = max(POS)) %>%

    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len) %>%

    # Add this info to the initial dataset
    left_join(df, ., by=c("CHR"="CHR")) %>%

    # Add a cumulative position of each SNP
    arrange(CHR, POS) %>%
    mutate(BPcum = POS + tot) %>%

    # Add highlight and annotation information
    mutate(is_highlight = ifelse(SNP %in% snpsOfInterest, "yes", "no"))


### Prepare X-axis
axisdf <- df_plot %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


### Plot
p1 <- ggplot(df_plot, aes(x=BPcum, y=-log10(PVAL))) +

    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
    scale_color_manual(values = rep(c(color1, color2), n_chr)) +
    
    # Add horizontal line (Genome-wide significant)
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red", linewidth = 0.8) +
    
    # Add horizontal line (Suggestive)
    geom_hline(yintercept = -log10(1e-5), linetype = "dashed", color = "blue", linewidth = 0.8) +
    
    # custom X axis:
    scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(limits = c(0, -log10(min(df_plot$PVAL))), expand = c(0, 1) ) +     # remove space between plot area and x axis
    
    # Add highlighted points
    geom_point(data=subset(df_plot, is_highlight=="yes"), color=color_annotate, size=1.5) +

    # Label name
    labs(x = "Chromosome",
        y = expression(bold("-log"[10](italic(P))))
        ) +
    
    # Custom the theme:
    theme_bw() +
    theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = 14, face="bold"),
        axis.text = element_text(size = 14)
    )

if(img_type == "pdf"){
    ggsave(paste0(outf, ".", img_type), p1, 
    scale = 2, device = "pdf",
    width = width, height = height,
    units = units)
} else{
    ggsave(paste0(outf, ".", img_type), p1, 
    scale = 2, dpi = dpi,
    width = width, height = height,
    units = units)
}



