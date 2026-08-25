# :: Independent signal determination ::
# - Sanghyeon Park
# - 2024.08.15
# - Project: Atherosclerosis GWAS
# Version 3 checked by ChatGPT

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(argparse)
})

options(error = traceback)

#######################
### Notes
# GWAS1 lead SNPs are tested against nearby significant GWAS2 lead SNPs.
# LD is computed chromosome by chromosome.
#######################

#######################
# Parser arguments
#######################
parser <- argparse::ArgumentParser(
  description = ":: SNP independence by P value and LD ::",
  formatter_class = "argparse.ArgumentDefaultsHelpFormatter"
)

parser$add_argument("--trait1", required = TRUE, help = "Trait name for GWAS1.")
parser$add_argument("--trait2", required = TRUE, help = "Trait name for GWAS2.")

# GWAS1 SNP list
parser$add_argument("--gwas1-snplist", dest = "file_gwas1_snplist", required = TRUE,
                    help = "Path to GWAS1 SNP list.")
parser$add_argument("--delim-gwas1-snplist", dest = "delim_gwas1_snplist", required = TRUE,
                    help = "Delimiter: tab, whitespace, comma")
parser$add_argument("--snp-col-gwas1-snplist", dest = "snp_gwas1_snplist", required = TRUE,
                    help = "GWAS1 SNP column name in SNP list.")

# GWAS2 SNP list
parser$add_argument("--gwas2-snplist", dest = "file_gwas2_snplist", required = TRUE,
                    help = "Path to GWAS2 SNP list.")
parser$add_argument("--delim-gwas2-snplist", dest = "delim_gwas2_snplist", required = TRUE,
                    help = "Delimiter: tab, whitespace, comma")
parser$add_argument("--snp-col-gwas2-snplist", dest = "snp_gwas2_snplist", required = TRUE,
                    help = "GWAS2 SNP column name in SNP list.")

# GWAS1 summary statistics
parser$add_argument("--gwas1", dest = "file_gwas1", required = TRUE,
                    help = "Path to GWAS1 summary statistics.")
parser$add_argument("--delim-gwas1", dest = "delim_gwas1", required = TRUE,
                    help = "Delimiter: tab, whitespace, comma")
parser$add_argument("--snp-col-gwas1", dest = "snp_gwas1", required = TRUE,
                    help = "GWAS1 SNP column name.")
parser$add_argument("--chr-col-gwas1", dest = "chr_gwas1", required = TRUE,
                    help = "GWAS1 chromosome column name.")
parser$add_argument("--pos-col-gwas1", dest = "pos_gwas1", required = TRUE,
                    help = "GWAS1 position column name.")
parser$add_argument("--pval-col-gwas1", dest = "p_gwas1", required = TRUE,
                    help = "GWAS1 P value column name.")

# GWAS2 summary statistics
parser$add_argument("--gwas2", dest = "file_gwas2", required = TRUE,
                    help = "Path to GWAS2 summary statistics.")
parser$add_argument("--delim-gwas2", dest = "delim_gwas2", required = TRUE,
                    help = "Delimiter: tab, whitespace, comma")
parser$add_argument("--snp-col-gwas2", dest = "snp_gwas2", required = TRUE,
                    help = "GWAS2 SNP column name.")
parser$add_argument("--chr-col-gwas2", dest = "chr_gwas2", required = TRUE,
                    help = "GWAS2 chromosome column name.")
parser$add_argument("--pos-col-gwas2", dest = "pos_gwas2", required = TRUE,
                    help = "GWAS2 position column name.")
parser$add_argument("--pval-col-gwas2", dest = "p_gwas2", required = TRUE,
                    help = "GWAS2 P value column name.")

# Reference panel
parser$add_argument("--reference-panel", dest = "reference_panel", required = FALSE, default = "1kG",
                    help = "Reference panel: 1kG, UKB_random_10K")

# Genome build
parser$add_argument("--genome-build", dest = "genome_build", required = FALSE, default = "GRCh37",
                    help = "Genome build. Currently only GRCh37 is supported.")

# Thresholds
parser$add_argument("--r2-threshold", dest = "r2_threshold", required = TRUE, type = "double",
                    help = "R2 threshold for non independence.")
parser$add_argument("--window", required = TRUE, type = "integer",
                    help = "One sided window in kb.")
parser$add_argument("--pval-threshold", dest = "p_threshold", required = TRUE, type = "double",
                    help = "P value threshold for GWAS2 candidate SNPs.")

# Parallel
parser$add_argument("--thread", required = FALSE, default = 1, type = "integer",
                    help = "Number of threads.")

# Output
parser$add_argument("--dir-out", dest = "dir_out", required = FALSE, default = "NA",
                    help = "Output directory. If NA, current working directory is used.")
parser$add_argument("--prefix-out", dest = "pref_out", required = FALSE, default = "NA",
                    help = "Output prefix.")
parser$add_argument("--delim-out", dest = "delim_out", required = FALSE, default = "comma",
                    help = "Output delimiter: tab, whitespace, comma")

args <- parser$parse_args()

#######################
# Utility
#######################
plink_exe <- "/data1/sanghyeon/wonlab_contribute/combined/software/plink/plink"

reference_panel_list <- c(
  "1kG" = "/data1/sanghyeon/wonlab_contribute/combined/data_common/reference_panel/1kGp3/EUR/CHR/reference.1kG.EUR.maf_0.005.geno_0.02.chr",
  "UKB_random_10K" = "/data1/sanghyeon/wonlab_contribute/combined/data_common/reference_panel/UKB_random10k_noQC/CHR/UKB_random10k_noqc_chr"
)

delim_mapper <- c(
  "tab" = "\t",
  "whitespace" = " ",
  "comma" = ","
)

ext_mapper <- c(
  "tab" = "tsv",
  "whitespace" = "txt",
  "comma" = "csv"
)

validate_choice <- function(x, choices, arg_name) {
  if (!x %in% choices) {
    stop(sprintf("Invalid value for %s: %s. Allowed values: %s",
                 arg_name, x, paste(choices, collapse = ", ")))
  }
}

check_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("%s does not exist: %s", label, path))
  }
}

check_cols <- function(df, cols, df_name) {
  missing_cols <- setdiff(cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns in %s: %s",
                 df_name, paste(missing_cols, collapse = ", ")))
  }
}

dedup_by_first <- function(df, snp_col) {
  if (anyDuplicated(df[[snp_col]]) > 0) {
    warning(sprintf("Duplicated SNP IDs found in %s. Keeping the first occurrence.", snp_col))
    df <- df %>% distinct(.data[[snp_col]], .keep_all = TRUE)
  }
  df
}

compute_snp_ld <- function(snp1, snp2, plink_exe, file_reference, chr, tmp_dir = tempdir()) {
  if (is.na(snp1) || is.na(snp2) || snp1 == "" || snp2 == "") {
    return(data.frame(ref_snp = snp1, other_snp = snp2, R2 = NA_real_, stringsAsFactors = FALSE))
  }

  if (identical(snp1, snp2)) {
    return(data.frame(ref_snp = snp1, other_snp = snp2, R2 = 1, stringsAsFactors = FALSE))
  }

  tmp_prefix <- file.path(
    tmp_dir,
    sprintf("plink_ld_chr%s_%s_%s_%s",
            chr, Sys.getpid(), gsub("[^A-Za-z0-9_.]", "_", snp1), gsub("[^A-Za-z0-9_.]", "_", snp2))
  )

  cmd <- paste(
    shQuote(plink_exe),
    "--bfile", shQuote(paste0(file_reference, chr)),
    "--ld", shQuote(snp1), shQuote(snp2),
    "--out", shQuote(tmp_prefix)
  )

  plink_out <- tryCatch(
    system(cmd, intern = TRUE, ignore.stderr = TRUE),
    error = function(e) character(0)
  )

  ld_result <- grep("R-sq", plink_out, value = TRUE)

  r2 <- NA_real_
  if (length(ld_result) > 0) {
    ld_line <- ld_result[1]
    ld_match <- regmatches(ld_line, regexpr("R-sq\\s+=\\s+([0-9.]+([eE][-+]?[0-9]+)?)", ld_line))
    if (length(ld_match) > 0 && nzchar(ld_match)) {
      r2 <- suppressWarnings(as.numeric(sub(".*=\\s+", "", ld_match)))
    }
  }

  unlink(paste0(tmp_prefix, c(".log", ".nosex")), force = TRUE)

  data.frame(ref_snp = snp1, other_snp = snp2, R2 = r2, stringsAsFactors = FALSE)
}

#######################
# Validate inputs
#######################
validate_choice(args$delim_gwas1_snplist, names(delim_mapper), "delim_gwas1_snplist")
validate_choice(args$delim_gwas2_snplist, names(delim_mapper), "delim_gwas2_snplist")
validate_choice(args$delim_gwas1, names(delim_mapper), "delim_gwas1")
validate_choice(args$delim_gwas2, names(delim_mapper), "delim_gwas2")
validate_choice(args$delim_out, names(delim_mapper), "delim_out")
validate_choice(args$reference_panel, names(reference_panel_list), "reference_panel")

if (args$genome_build != "GRCh37") {
  stop("Only GRCh37 is currently supported.")
}

check_file_exists(args$file_gwas1_snplist, "GWAS1 SNP list")
check_file_exists(args$file_gwas2_snplist, "GWAS2 SNP list")
check_file_exists(args$file_gwas1, "GWAS1 summary statistics")
check_file_exists(args$file_gwas2, "GWAS2 summary statistics")
check_file_exists(plink_exe, "PLINK executable")

#######################
# Parse arguments
#######################
trait1 <- args$trait1
trait2 <- args$trait2

file_gwas1_snplist <- args$file_gwas1_snplist
delim_gwas1_snplist <- args$delim_gwas1_snplist
snp_gwas1_snplist <- args$snp_gwas1_snplist

file_gwas2_snplist <- args$file_gwas2_snplist
delim_gwas2_snplist <- args$delim_gwas2_snplist
snp_gwas2_snplist <- args$snp_gwas2_snplist

file_gwas1 <- args$file_gwas1
delim_gwas1 <- args$delim_gwas1
snp_gwas1 <- args$snp_gwas1
chr_gwas1 <- args$chr_gwas1
pos_gwas1 <- args$pos_gwas1
p_gwas1 <- args$p_gwas1

file_gwas2 <- args$file_gwas2
delim_gwas2 <- args$delim_gwas2
snp_gwas2 <- args$snp_gwas2
chr_gwas2 <- args$chr_gwas2
pos_gwas2 <- args$pos_gwas2
p_gwas2 <- args$p_gwas2

reference_panel <- args$reference_panel
r2_threshold <- args$r2_threshold
window <- args$window
p_threshold <- args$p_threshold
n_thread <- args$thread
dir_out <- args$dir_out
pref_out <- args$pref_out
delim_out <- args$delim_out

if (dir_out == "NA") dir_out <- "."
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

if (pref_out == "NA") {
  pref_out <- paste0("snp_independence.", tools::file_path_sans_ext(basename(file_gwas1_snplist)))
}

out_file <- file.path(
  dir_out,
  paste0(pref_out, ".compare_", trait1, ".to_", trait2, ".", ext_mapper[[delim_out]])
)

#######################
# Read data
#######################
df_gwas1_snplist <- fread(
  file_gwas1_snplist,
  sep = delim_mapper[[delim_gwas1_snplist]],
  data.table = FALSE,
  nThread = n_thread,
  showProgress = FALSE
)
check_cols(df_gwas1_snplist, snp_gwas1_snplist, "df_gwas1_snplist")

df_gwas1_snplist <- df_gwas1_snplist %>%
  select(all_of(snp_gwas1_snplist)) %>%
  rename(leadSNP = all_of(snp_gwas1_snplist)) %>%
  distinct()

df_gwas1 <- fread(
  file_gwas1,
  sep = delim_mapper[[delim_gwas1]],
  data.table = FALSE,
  nThread = n_thread,
  showProgress = FALSE
)
check_cols(df_gwas1, c(snp_gwas1, chr_gwas1, pos_gwas1, p_gwas1), "df_gwas1")

df_gwas1 <- df_gwas1 %>%
  select(all_of(c(snp_gwas1, chr_gwas1, pos_gwas1, p_gwas1))) %>%
  rename(
    SNP = all_of(snp_gwas1),
    CHR = all_of(chr_gwas1),
    POS = all_of(pos_gwas1),
    P = all_of(p_gwas1)
  ) %>%
  dedup_by_first("SNP")

df_gwas1_snplist <- df_gwas1_snplist %>%
  left_join(df_gwas1, by = c("leadSNP" = "SNP"))

df_gwas2_snplist <- fread(
  file_gwas2_snplist,
  sep = delim_mapper[[delim_gwas2_snplist]],
  data.table = FALSE,
  nThread = n_thread,
  showProgress = FALSE
)
check_cols(df_gwas2_snplist, snp_gwas2_snplist, "df_gwas2_snplist")

df_gwas2_snplist <- df_gwas2_snplist %>%
  select(all_of(snp_gwas2_snplist)) %>%
  rename(leadSNP = all_of(snp_gwas2_snplist)) %>%
  distinct()

df_gwas2 <- fread(
  file_gwas2,
  sep = delim_mapper[[delim_gwas2]],
  data.table = FALSE,
  nThread = n_thread,
  showProgress = FALSE
)
check_cols(df_gwas2, c(snp_gwas2, chr_gwas2, pos_gwas2, p_gwas2), "df_gwas2")

df_gwas2 <- df_gwas2 %>%
  select(all_of(c(snp_gwas2, chr_gwas2, pos_gwas2, p_gwas2))) %>%
  rename(
    SNP = all_of(snp_gwas2),
    CHR = all_of(chr_gwas2),
    POS = all_of(pos_gwas2),
    P = all_of(p_gwas2)
  ) %>%
  dedup_by_first("SNP")

df_gwas2_snplist <- df_gwas2_snplist %>%
  left_join(df_gwas2, by = c("leadSNP" = "SNP"))

if (any(is.na(df_gwas1_snplist$CHR) | is.na(df_gwas1_snplist$POS))) {
  warning("Some GWAS1 lead SNPs could not be mapped to CHR/POS.")
}
if (any(is.na(df_gwas2_snplist$CHR) | is.na(df_gwas2_snplist$POS))) {
  warning("Some GWAS2 lead SNPs could not be mapped to CHR/POS.")
}

#######################
# LD
#######################
file_reference_panel <- reference_panel_list[[reference_panel]]

cl <- makeCluster(n_thread)
registerDoParallel(cl)
on.exit({
  try(stopCluster(cl), silent = TRUE)
}, add = TRUE)

all_results <- list()
res_idx <- 1L

for (chr in 1:22) {
  message(sprintf("Processing chromosome %d ...", chr))

  df_gwas1_chr <- df_gwas1_snplist %>%
    filter(CHR == chr, !is.na(POS), !is.na(leadSNP))

  df_gwas2_chr <- df_gwas2_snplist %>%
    filter(CHR == chr, !is.na(POS), !is.na(leadSNP), !is.na(P))

  if (nrow(df_gwas1_chr) == 0) next

  unique_leadsnp <- unique(df_gwas1_chr$leadSNP)

  for (ref_snp in unique_leadsnp) {
    ref_pos <- df_gwas1_chr %>%
      filter(leadSNP == ref_snp) %>%
      pull(POS) %>%
      unique()

    if (length(ref_pos) != 1L || is.na(ref_pos)) {
      warning(sprintf("Skipping %s on chr %d because POS is missing or duplicated.", ref_snp, chr))
      all_results[[res_idx]] <- data.frame(
        SNP_gwas1 = ref_snp,
        SNP_gwas2_tested = "NA",
        SNP_gwas2_tested_pval = "NA",
        SNP_gwas2_tested_R2 = "NA",
        In_LD = NA,
        SNP_gwas2_in_LD = "NA",
        SNP_gwas2_in_LD_pval = "NA",
        SNP_gwas2_in_LD_R2 = "NA",
        LeadSNP_gwas2_in_ld = "NA",
        stringsAsFactors = FALSE
      )
      res_idx <- res_idx + 1L
      next
    }

    df_candidates <- df_gwas2_chr %>%
      filter(
        POS >= ref_pos - (window * 1000),
        POS <= ref_pos + (window * 1000),
        P < p_threshold
      ) %>%
      distinct(leadSNP, .keep_all = TRUE)

    if (nrow(df_candidates) == 0) {
      all_results[[res_idx]] <- data.frame(
        SNP_gwas1 = ref_snp,
        SNP_gwas2_tested = "NA",
        SNP_gwas2_tested_pval = "NA",
        SNP_gwas2_tested_R2 = "NA",
        In_LD = FALSE,
        SNP_gwas2_in_LD = "NA",
        SNP_gwas2_in_LD_pval = "NA",
        SNP_gwas2_in_LD_R2 = "NA",
        LeadSNP_gwas2_in_ld = "NA",
        stringsAsFactors = FALSE
      )
      res_idx <- res_idx + 1L
      next
    }

    leadsnp_to_compare <- unique(df_candidates$leadSNP)

    plink_ld_result <- foreach(
      leadsnp_gwas2 = leadsnp_to_compare,
      .combine = bind_rows,
      .packages = c("dplyr")
    ) %dopar% {
      compute_snp_ld(
        snp1 = ref_snp,
        snp2 = leadsnp_gwas2,
        plink_exe = plink_exe,
        file_reference = file_reference_panel,
        chr = chr
      )
    }

    df_tmp <- df_candidates %>%
      select(leadSNP, P) %>%
      left_join(
        plink_ld_result %>% select(other_snp, R2),
        by = c("leadSNP" = "other_snp")
      )

    snp_gwas2_tested <- df_tmp$leadSNP
    snp_gwas2_tested_pval <- df_tmp$P
    snp_gwas2_tested_r2 <- as.numeric(df_tmp$R2)

    idx_in_ld <- which(!is.na(snp_gwas2_tested_r2) & snp_gwas2_tested_r2 > r2_threshold)

    if (length(idx_in_ld) == 0) {
      is_in_LD <- FALSE
      snp_gwas2_in_ld <- "NA"
      snp_gwas2_in_ld_pval <- "NA"
      snp_gwas2_in_ld_r2 <- "NA"
      snp_gwas2_in_ld_leadsnp <- "NA"
    } else {
      is_in_LD <- TRUE
      snp_gwas2_in_ld <- paste(snp_gwas2_tested[idx_in_ld], collapse = ";")
      snp_gwas2_in_ld_pval <- paste(snp_gwas2_tested_pval[idx_in_ld], collapse = ";")
      snp_gwas2_in_ld_r2 <- paste(snp_gwas2_tested_r2[idx_in_ld], collapse = ";")
      snp_gwas2_in_ld_leadsnp <- paste(unique(snp_gwas2_tested[idx_in_ld]), collapse = ";")
    }

    all_results[[res_idx]] <- data.frame(
      SNP_gwas1 = ref_snp,
      SNP_gwas2_tested = paste(snp_gwas2_tested, collapse = ";"),
      SNP_gwas2_tested_pval = paste(snp_gwas2_tested_pval, collapse = ";"),
      SNP_gwas2_tested_R2 = paste(snp_gwas2_tested_r2, collapse = ";"),
      In_LD = is_in_LD,
      SNP_gwas2_in_LD = snp_gwas2_in_ld,
      SNP_gwas2_in_LD_pval = snp_gwas2_in_ld_pval,
      SNP_gwas2_in_LD_R2 = snp_gwas2_in_ld_r2,
      LeadSNP_gwas2_in_ld = snp_gwas2_in_ld_leadsnp,
      stringsAsFactors = FALSE
    )
    res_idx <- res_idx + 1L
  }

  df_ld_result_gwas1 <- bind_rows(all_results)

  write.table(
    df_ld_result_gwas1,
    file = out_file,
    sep = delim_mapper[[delim_out]],
    row.names = FALSE,
    quote = FALSE
  )
}

df_ld_result_gwas1 <- bind_rows(all_results)

write.table(
  df_ld_result_gwas1,
  file = out_file,
  sep = delim_mapper[[delim_out]],
  row.names = FALSE,
  quote = FALSE
)