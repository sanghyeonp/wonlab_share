##### To do
# - OR 인 경우, log(OR)을 만들기
###########

library(data.table)
library(dplyr)

gwas <- "/data1/sanghyeon/Projects/mr_drug_repurposing/data/phenocode-709.3.tsv.gz"
delim <- "\t"

chr_col <- "chrom"
pos_col <- "pos"
ref_col <- "ref"
alt_col <- "alt"
freq_col <- "af" # default = "NA"
snp_col <- "rsids" # default = "NA", if NA make variant column
beta_col <- "beta"
se_col <- "sebeta"
pval_col <- "pval" # default = "NA", if NA, deduce from beta and se column

compute_p <- TRUE

### Function
"%notin%" <- function(x,y) !("%in%"(x,y))

### Read GWAS
df <- fread(gwas, sep=delim)

### Common modification
# Make index column
df$Index <- c(1:nrow(df))

# Capitalize alleles
df <- df %>%
  mutate(!!as.name(ref_col) := toupper(!!as.name(ref_col))) %>%
  mutate(!!as.name(alt_col) := toupper(!!as.name(alt_col)))


################## <<< Check rows to drop
### 여기에서는 NA 값 및 wrong data type 제거.
cat("\n:: CHECK - NULL and mismatching data type ::")

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
null_to_drop <- sort(unique(c(null_chr,
                              null_pos,
                              non_numeric_pos,
                              null_ref,
                              null_alt,
                              null_freq,
                              non_numeric_freq,
                              null_beta,
                              non_numeric_beta,
                              null_se,
                              non_numeric_se
                              )))


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

################## >>> Check rows to drop

################## <<< Check rows with insensible values
### 여기에서는 해당 column에 insensible 값 제거.

cat("\n:: Check for invalid values ::")

# Make a copy of dataframe to conduct checks
df_check <- df

# Chromosome column
df_check[[chr_col]] <- as.character(df_check[[chr_col]])
non_autosome_chr <- which(df_check[[chr_col]] %notin% as.character(c(1:22)))

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
if(or_col != "NA"){
  pass
  # CODE: 0 보다 작은건 제거.
}

# SE column
invalid_se <- which((df_check[[se_col]] < 0) | (df_check[[se_col]] == Inf))

# PVAL column (여기서는 report만 하고, 실제로 drop하지 않음. 나중에 beta와 se로 compute 할 예정.)
invalid_pval <- which((df_check[[pval_col]] <= 0) | (df_check[[pval_col]] > 1))

#### Drop insensible values
# Merge vectors
invalid_to_drop <- sort(unique(c(non_autosome_chr,
                                 invalid_ref,
                                 invalid_alt,
                                 invalid_freq,
                                 invalid_beta,
                                 invalid_se)))

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

################## >>> Check rows with insensible values

################## <<< Check if BETA is sensible
cat("\n:: Check if BETA is sensible ::")
beta_median <- median(df[[beta_col]])
tolerance <- 0.1
expected_beta_median <- 0

if (abs(beta_median - expected_beta_median) > tolerance){
  cat(paste0("\nMedian BETA is ", beta_median, ", which seems insensible."))
} else{
  cat(paste0("\nMedian BETA is ", beta_median, ", which seems sensible."))
}

################## >>> Check if BETA is sensible


################## <<< Recompute P-value from BETA and SE
new_pval <- 2*pnorm(abs(df[[beta_col]] / df[[se_col]], lower.tail=FALSE))
################## >>> Recompute P-value from BETA and SE


cat("\n")