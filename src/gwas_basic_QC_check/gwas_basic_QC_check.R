# Sanghyeon Park
# 2024.06.04
# Project: Atherosclerosis GWAS

# Workflow
# 1. Subset specified column names
# 2. Generate SNP.tmp column
# 3. Check missing values
# 4. Check invalid values
# 5. Check SNP duplicates
# 6. Check multi-allelic SNPs
# 7. Check INDEL
# 8. Check monomorphic SNPs
# 9. Check non-automsomal SNPs
# 10. Check SNP AF < threshold
# 11. Check SNP INFO < threshold
# 12. Additional criteria given

#### Version control (Current version = 2)
# v1: Initial version
# v2: Allow input of addtional criteria by the users
####

source("/data1/sanghyeon/wonlab_contribute/combined/src/gwas_basic_QC_check/fnc.gwas_basic_QC_check.R")
library(data.table)
library(dplyr)
library(argparse)
options(error=traceback)

# Parse command line arguments
parser <- ArgumentParser()
parser$add_argument("--gwas", dest="file_gwas", type="character", required=TRUE,
                    help="Input file GWAS.")

parser$add_argument("--snp-col", dest="snp_col", type="character", required=FALSE, default="",
                    help="")
parser$add_argument("--chr-col", dest="chr_col", type="character", required=FALSE, default="",
                    help="")
parser$add_argument("--pos-col", dest="pos_col", type="character", required=FALSE, default="",
                    help="")
parser$add_argument("--a1-col", dest="a1_col", type="character", required=FALSE, default="",
                    help="")
parser$add_argument("--a2-col", dest="a2_col", type="character", required=FALSE, default="",
                    help="")
parser$add_argument("--af-col", dest="af_col", type="character", required=FALSE, default="",
                    help="")
parser$add_argument("--effect-col", dest="effect_col", type="character", required=FALSE, default="",
                    help="")
parser$add_argument("--se-col", dest="se_col", type="character", required=FALSE, default="",
                    help="")
parser$add_argument("--p-col", dest="p_col", type="character", required=FALSE, default="",
                    help="")
parser$add_argument("--info-col", dest="info_col", type="character", required=FALSE, default="",
                    help="")

parser$add_argument("--thres-af", dest="threshold_AF", type="numeric", required=FALSE, default=0.01,
                    help="")
parser$add_argument("--thres-info", dest="threshold_INFO", type="numeric", required=FALSE, default=0.9,
                    help="")

parser$add_argument("--additional-criteria", dest="additional_criteria", type="character", required=FALSE, default="",
                    help="")
parser$add_argument("--additional-col", dest="additional_col", type="character", required=FALSE, default="",
                    help="")

parser$add_argument("--pref", dest="prefix", type="character", required=FALSE, default="",
                    help="")
parser$add_argument("--save-snp-count", dest="save_snp_count", type="character", required=FALSE, default="",
                    help="")
parser$add_argument("--save-snp-list", dest="save_snp_list", type="character", required=FALSE, default="",
                    help="")

parser$add_argument("--n-thread", dest="n_thread", type="integer", required=FALSE, default=1, 
                    help="")

args <- parser$parse_args()

empty_to_NA <- function(x) {
    if (nchar(x) == 0) {
        return(NA)
    } else{
        return(x)
    }
}

args$snp_col <- empty_to_NA(args$snp_col)
args$chr_col <- empty_to_NA(args$chr_col)
args$pos_col <- empty_to_NA(args$pos_col)
args$a1_col <- empty_to_NA(args$a1_col)
args$a2_col <- empty_to_NA(args$a2_col)
args$af_col <- empty_to_NA(args$af_col)
args$effect_col <- empty_to_NA(args$effect_col)
args$se_col <- empty_to_NA(args$se_col)
args$p_col <- empty_to_NA(args$p_col)
args$info_col <- empty_to_NA(args$info_col)

args$additional_col <- empty_to_NA(args$additional_col)
args$additional_criteria <- empty_to_NA(args$additional_criteria)
args$save_snp_count <- empty_to_NA(args$save_snp_count)
args$save_snp_list <- empty_to_NA(args$save_snp_list)


df.gwas <- fread(args$file_gwas, data.table=F, nThread=args$n_thread)

basic_GWAS_filter_criteria(df_gwas=df.gwas, 
                        snp_col=args$snp_col, chr_col=args$chr_col, pos_col=args$pos_col, 
                        a1_col=args$a1_col, a2_col=args$a2_col, af_col=args$af_col, 
                        effect_col=args$effect_col, se_col=args$se_col, p_col=args$p_col, info_col=args$info_col,
                        threshold.AF=args$threshold_AF, threshold.INFO=args$threshold_INFO, 
                        additional_criteria=args$additional_criteria, additional_columns=args$additional_col,
                        prefix=args$prefix, 
                        save.snp_count=args$save_snp_count, save.snp_list=args$save_snp_list)
