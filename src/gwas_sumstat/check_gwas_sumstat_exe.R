list.of.packages <- c("argparse", "dplyr", "rstudioapi")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(rstudioapi)
    require(argparse)
    library(ggplot2)
})

### Define parser arguments
parser <- argparse::ArgumentParser(description=":: Check downloaded GWAS summary statistics ::", 
                                formatter_class="argparse.ArgumentDefaultsHelpFormatter")

parser$add_argument("--gwas", required=TRUE,
                    help="Path to the GWAS summary statistics.")
parser$add_argument("--delim", required=FALSE, default="tab",
                    help="Specify the delimiter for the output file. Options = [tab, comma, whitespace]. Default = tab.")

parser$add_argument("--chr-col", dest="chr_col", required=FALSE, default="NA",
                    help="Specify the column name for chromosome.")
parser$add_argument("--pos-col", dest="pos_col", required=FALSE, default="NA",
                    help="Specify the column name for base position.")
parser$add_argument("--ref-col", dest="ref_col", required=FALSE, default="NA",
                    help="Specify the column name for reference allele.")
parser$add_argument("--alt-col", dest="alt_col", required=FALSE, default="NA",
                    help="Specify the column name for alt allele.")
parser$add_argument("--freq-col", dest="freq_col", required=FALSE, default="NA",
                    help="Specify the column name for allele frequency.")
parser$add_argument("--snp-col", dest="snp_col", required=FALSE, default="NA",
                    help="Specify the column name for SNP ID.")
parser$add_argument("--beta-col", dest="beta_col", required=FALSE, default="NA",
                    help="Specify the column name for Effect.")
parser$add_argument("--se-col", dest="se_col", required=FALSE, default="NA",
                    help="Specify the column name for SE.")
parser$add_argument("--pval-col", dest="pval_col", required=FALSE, default="NA",
                    help="Specify the column name for P-value.")

parser$add_argument("--compute-pval", dest="compute_pval", action="store_true", 
                    help="Specify to compute P-value from BETA and SE. Default=FALSE.")

parser$add_argument("--name", required=FALSE, default="GWAS",
                    help="Specify the unique name of the specified GWAS summary statistics. Default = 'GWAS'.")
################## <<< Defined functions and variables
delim_map <- list(tab = "\t", comma = ",", whitespace = " ")
"%notin%" <- function(x,y) !("%in%"(x,y))
################## >>> Defined functions and variables

################## <<< Get parser arguments
args <- parser$parse_args()
gwas <- args$gwas
delim <- delim_map[[args$delim]]

chr_col <- args$chr_col
pos_col <- args$pos_col
ref_col <- args$ref_col
alt_col <- args$alt_col
freq_col <-args$freq_col
snp_col <- args$snp_col
beta_col <- args$beta_col
se_col <- args$se_col
pval_col <- args$pval_col

compute_pval <- args$compute_pval

uniq_name <- args$name
################## >>> Get parser arguments

################## <<< Read GWAS
cat("\n:: Reading GWAS summary statistics ::")
df <- fread(gwas, sep=delim)
colname_list <- colnames(df)

cat(paste0("\n\t Path: ", gwas))
cat(paste0("\n\t Columns: ", paste(colname_list, collapse=", ")))

# Check if all the specified columns are in GWAS
for (col in c(chr_col, pos_col, ref_col, alt_col, freq_col, snp_col, beta_col, se_col, pval_col)){
  if ((col != 'NA') & (col %notin% colname_list)){
    stop(paste0("Column '", col, "' is not present in the GWAS summary statistics."))
  }
}

cat("\n")
################## >>> Read GWAS

################## <<< Common modification
# Make index column
df$Index <- c(1:nrow(df))

# Capitalize alleles
df <- df %>%
  mutate(!!as.name(ref_col) := toupper(!!as.name(ref_col))) %>%
  mutate(!!as.name(alt_col) := toupper(!!as.name(alt_col)))
################## >>> Common modification

################## <<< Check rows to drop
### 여기에서는 NA 값 및 wrong data type 제거.
cat("\n:: CHECK - NULL and mismatching data type ::")

null_chr <- NULL
null_pos <- NULL
non_numeric_pos <- NULL
null_ref <- NULL
null_alt <- NULL
null_freq <- NULL
non_numeric_freq <- NULL
null_beta <- NULL
non_numeric_beta <- NULL
null_se <- NULL
non_numeric_se <- NULL

# Make a copy of dataframe to conduct checks
df_check <- df

# Chromosome column
null_chr <- which(is.na(df_check[[chr_col]]))

# Position column
null_pos <- which(is.na(df_check[[pos_col]]))
non_numeric_pos <- df_check[[pos_col]][!grepl("^\\d+$", df_check[[pos_col]])]

# Check allele columns
null_ref <- which(is.na(df_check[[ref_col]]))
null_alt <- which(is.na(df_check[[alt_col]]))

# Frequency column
if(freq_col != "NA"){
  null_freq <- which(is.na(df_check[[freq_col]]))
  non_numeric_freq <- which(!grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", df_check[[freq_col]]))
}

# Beta column
if(beta_col != "NA"){
  null_beta <- which(is.na(df_check[[beta_col]]))
  non_numeric_beta <- which(!grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", df_check[[beta_col]]))
}

# OR column
# 여기에 OR check

# SE column
null_se <- which(is.na(df_check[[se_col]]))
non_numeric_se <- which(!grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", df_check[[se_col]]))

# PVAL column
null_pval <- which(is.na(df_check[[pval_col]]))
non_numeric_pval <- which(!grepl("^[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?$", df_check[[pval_col]]))

#### Drop NA
# Merge vectors
temp_v <- NULL
for (ele in list(null_chr,
                  null_pos,
                  non_numeric_pos,
                  null_ref,
                  null_alt,
                  null_freq,
                  non_numeric_freq,
                  null_beta,
                  non_numeric_beta,
                  null_se,
                  non_numeric_se)){
                    if (!is.null(ele)){
                      temp_v <- c(temp_v, ele)
                    }
                  }
null_to_drop <- sort(unique(temp_v))


# Report summary
cat("\nREPORT SUMMARY:")
cat(paste0("\n\tNumber of rows with NULL Chromosome: ", length(null_chr)))
cat(paste0("\n\tNumber of rows with NULL Base position: ", length(null_pos)))
cat(paste0("\n\tNumber of rows with non-numeric Base position: ", length(non_numeric_pos)))
cat(paste0("\n\tNumber of rows with NULL Reference allele: ", length(null_ref)))
cat(paste0("\n\tNumber of rows with NULL Alternative allele: ", length(null_alt)))
cat(paste0("\n\tNumber of rows with NULL Allele frequency: ", length(null_freq)))
cat(paste0("\n\tNumber of rows with non-numeric Allele frequency: ", length(non_numeric_freq)))
cat(paste0("\n\tNumber of rows with NULL Effect estimate: ", length(null_beta)))
cat(paste0("\n\tNumber of rows with non-numeric Effect estimate: ", length(non_numeric_beta)))
cat(paste0("\n\tNumber of rows with NULL Standard error: ", length(null_se)))
cat(paste0("\n\tNumber of rows with non-numeric Standard error: ", length(non_numeric_se)))
cat(paste0("\n\tNumber of rows with NULL P-value: ", length(null_pval)))
cat(paste0("\n\tNumber of rows with non-numeric P-value: ", length(non_numeric_pval)))


# Remove null
if (length(null_to_drop) > 0){
  null_dropped_Index_ori <- df[null_to_drop, ]$Index
  df <- df[-null_to_drop, ]
}

cat("\n")
################## >>> Check rows to drop

################## <<< Check rows with insensible values
### 여기에서는 해당 column에 insensible 값 제거.
non_autosome_chr <- NULL
invalid_ref <- NULL
invalid_alt <- NULL
invalid_freq <- NULL
invalid_beta <- NULL
invalid_se <- NULL

cat("\n:: Check for invalid values ::")

# Make a copy of dataframe to conduct checks
df_check <- df

# Chromosome column
df_check[[chr_col]] <- as.character(df_check[[chr_col]])
unique_chr_list <- unique(df_check[[chr_col]])
cat(paste0("\n\t- Unique chromsomes: ", paste(unique_chr_list, collapse=", ")))

non_autosome_chr <- which(sum(df_check[[chr_col]] == as.character(c(1:22))) == 0)

# Allele columns
invalid_ref <- which(!grepl("^[ATCG]+$", as.character(df_check[[ref_col]])))
invalid_alt <- which(!grepl("^[ATCG]+$", as.character(df_check[[alt_col]])))

# Frequency column
if(freq_col != "NA"){
  invalid_freq <- which((df_check[[freq_col]] < 0) | (df_check[[freq_col]] > 1))
}

# Beta column
invalid_beta <- which(abs(df_check[[beta_col]]) == Inf)

# OR column
# if(or_col != "NA"){
#   pass
#   # CODE: 0 보다 작은건 제거.
# }

# SE column
invalid_se <- which((df_check[[se_col]] < 0) | (df_check[[se_col]] == Inf))

# PVAL column (여기서는 report만 하고, 실제로 drop하지 않음. 나중에 beta와 se로 compute 할 예정.)
invalid_pval <- which((df_check[[pval_col]] <= 0) | (df_check[[pval_col]] > 1))

#### Drop insensible values
# Merge vectors
temp_v <- NULL
for (ele in list(non_autosome_chr,
                  invalid_ref,
                  invalid_alt,
                  invalid_freq,
                  invalid_beta,
                  invalid_se)){
                    if (!is.null(ele)){
                      temp_v <- c(temp_v, ele)
                    }
                  }
invalid_to_drop <- sort(unique(temp_v))

# Report summary
cat("\nREPORT SUMMARY:")
cat(paste0("\n\tNumber of rows with non-autosomal Chromosome: ", length(non_autosome_chr)))
cat(paste0("\n\tNumber of rows with non-ATCG Reference allele: ", length(invalid_ref)))
cat(paste0("\n\tNumber of rows with non-ATCG Alternative allele: ", length(invalid_alt)))
cat(paste0("\n\tNumber of rows with invalid Allele frequency: ", length(invalid_freq)))
cat(paste0("\n\tNumber of rows with invalid Effect: ", length(invalid_beta)))
cat(paste0("\n\tNumber of rows with invalid Standard error: ", length(invalid_se)))
cat(paste0("\n\tNumber of rows with invalid P-value: ", length(invalid_pval)))

# Remove null
if (length(invalid_to_drop) > 0){
  invalid_dropped_Index_ori <- df[invalid_to_drop, ]$Index
  df <- df[-null_to_drop, ]
}

cat("\n")
################## >>> Check rows with insensible values

######################### <<< Statistics
cat("\n:: Statistics ::")
### SNPs
cat("\n- SNPs")
cat(paste0("\n\t- Number of SNPs: ", prettyNum(nrow(df), big.mark=",", scientific=FALSE)))

################## <<< Chromosomes
cat("\n- Chromosomes")
have_all_automsome <- TRUE
chr_not_present_list <- c()
unique_chr_list <- unique(df[[chr_col]])

for (chr in c(1:22)){
  # Check if autosome is present in the GWAS
  if (chr %notin% unique_chr_list){
    have_all_automsome <- FALSE
    chr_not_present_list <- c(chr_not_present_list, chr)
  }
}

if (have_all_automsome){
  cat("\n\t- All autosomes are present.")
} else{
  cat("\n\t- All autosomes are not present.")
  cat(paste0("\n\t- Absent chromsomes: ", paset(chr_not_present_list, collapse=", ")))
}



################## >>> Chromosomes

################## <<< Check BETA
cat("\n:: Check if BETA is sensible ::")
beta_median <- median(df[[beta_col]])
tolerance <- 0.1
expected_beta_median <- 0

if (abs(beta_median - expected_beta_median) > tolerance){
  cat(paste0("\nMedian BETA is ", beta_median, ", which seems insensible."))
} else{
  cat(paste0("\nMedian BETA is ", beta_median, ", which seems sensible."))
}

#! 여기 histogram
p_beta_hist <- ggplot(df, aes(x=!!as.name(beta_col))) +
  geom_histogram(aes(y=..density..), colour="black", fill="white") +
  geom_density(alpha=0.2, fill="#FF6666") +
  theme_minimal() +
  ylab("Density") +
  xlab("BETA") +
  theme(axis.text = element_text(size=13),
        axis.title = element_text(size=15, face="bold")
        )
ggsave(p_beta_hist, paste0("histogram_beta.", uniq_name, ".png"),
      device="png", dpi=300
      )
################## >>> Check BETA

cat("\n")
######################### >>> Statistics

################## <<< Recompute P-value from BETA and SE
new_pval <- 2*pnorm(abs(df[[beta_col]] / df[[se_col]]), lower.tail=FALSE)
################## >>> Recompute P-value from BETA and SE


cat("\n")