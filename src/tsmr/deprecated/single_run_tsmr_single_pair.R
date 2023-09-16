setwd("/data1/sanghyeon/Projects/SGLT2_CMC/src/tsmr_230913")

#########################################################
# Two-sample Mendelian randomization analysis AUTO
# Sanghyeon Park
# 2023.09.14
#########################################################

library(data.table)
library(dplyr)
library(TwoSampleMR)

### Exposure information
exposure_gwas <- "/data1/sanghyeon/Projects/SGLT2_CMC/gwas/ukb/annovar/Glycated_haemoglobin_HbA1c/annovarmapped.Glycated_haemoglobin_HbA1c.imp.gz"
exposure_delim <- "\t"
exposure_name <- "HbA1c"
exposure_cohort <- "UKB"
exposure_type <- "quantitative"
exposure_snp <- "rsID_annov"
exposure_beta <- "Effect"
exposure_se <- "StdErr"
exposure_eaf <- "MAF"
exposure_ea <- "ALT"
exposure_oa <- "REF"
exposure_p <- "P-value"
exposure_chr <- "#CHROM"
exposure_pos <- "POS"
## If exposure is binary.
# exposure_ncase <- 0
# exposure_ncontrol <- 0
## If exposure is quantitative.
exposure_n <- 363228

### Outcome information
outcome_gwas <- "/data1/sanghyeon/Projects/SGLT2_CMC/gwas/DIAGRAM/Mahajan.NatGenet2018b.T2D-noUKBB.European.ANNOVAR.txt"
outcome_delim <- "\t"
outcome_name <- "T2D"
outcome_cohort <- "DIAGRAM"
outcome_type <- "binary"
outcome_snp <- "rsID_annov"
outcome_beta <- "Beta"
outcome_se <- "SE"
outcome_eaf <- "EAF"
outcome_ea <- "EA"
outcome_oa <- "NEA"
outcome_p <- "Pvalue"
outcome_chr <- "Chr"
outcome_pos <- "Pos"
## If outcome is binary.
#  Total (74,124 T2D cases and 824,006 controls ) - UKB (19,119 T2D cases, 423,698 T2D controls) = 55005 T2D cases, 400308 T2D controls = 455313
outcome_ncase <- 55005
outcome_ncontrol <- 400308
## If outcome is quantitative.
# outcome_n <- 0

### Gene information
slc5a2_chr <- 16
slc5a2_start <- 31494444
slc5a2_end <- 31502090
window <- 500

### Clumping information
r2 <- 0.1
clump_window <- 500
pval_thres <- 5e-8

#################################################################################################################################


if (!file.exists(paste0(exposure_name, "_h_", outcome_name, "_GeneWindow", window, "kb.RDS"))){
  ##########################################
  # 1. Exposure data preparation
  ##########################################
  df_exp <- fread(exposure_gwas, sep=exposure_delim)

  ## Filter SLC5A2 region
  exp_dat <- df_exp %>% 
    select(all_of(c(exposure_snp, exposure_chr, exposure_pos, exposure_ea, exposure_oa, exposure_eaf, exposure_beta, exposure_se, exposure_p))) %>%
    rename(SNP := all_of(exposure_snp),
           CHR := all_of(exposure_chr),
           POS := all_of(exposure_pos),
           A1 := all_of(exposure_ea),
           A2 := all_of(exposure_oa),
           EAF := all_of(exposure_eaf),
           BETA := all_of(exposure_beta),
           SE := all_of(exposure_se),
           PVAL := all_of(exposure_p)
           ) %>%
    mutate(PVAL = as.numeric(PVAL)) %>%
    filter((CHR == slc5a2_chr) & (POS > slc5a2_start - (window * 1000)) & (POS <= slc5a2_end + (window * 1000)))

  write.table(exp_dat, paste0(exposure_name, "_SLC5A2_GeneWindow", window, "kb.tsv"),
              sep="\t", row.names=FALSE, quote = FALSE
  )

  ## Perform clumping
  exp_dat <- filter(exp_dat, PVAL < pval_thres)

  # Initialize an empty DataFrame
  exp_iv_dat <- data.frame()
  connection_trial <- 1

  # Keep running the code until there's at least one row in the DataFrame
  repeat {
    print(paste0("Clumping connection trial: ", connection_trial))
    exp_iv_dat <- clump_data(format_data(dat=exp_dat,
                                  type = "exposure",
                                  header = TRUE,
                                  snp_col = "SNP",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  effect_allele_col = "A1",
                                  other_allele_col = "A2",
                                  eaf_col = "EAF",
                                  pval_col = "PVAL",
                                  chr_col = "CHR",
                                  pos_col = "POS"), 
                    clump_r2 =r2,
                    clump_kb = clump_window,
                    pop = "EUR"
                    )
    connection_trial <- connection_trial + 1
    if (!is.data.frame(exp_iv_dat) || nrow(exp_iv_dat) > 0) {
      break
    }
  }

  exp_iv_dat$exposure <-paste0(exposure_name, "(", exposure_cohort, ")")
  if (exposure_type == "binary"){
    exp_iv_dat$ncase.exposure <- exposure_ncase
    exp_iv_dat$ncontrol.exposure <- exposure_ncontrol
  } else{
    exp_iv_dat$samplesize.exposure <- exposure_n
  }

  exp_iv_dat

  ##########################################
  # 2. IV validity check
  ##########################################

  ## F-statistics
  exp_iv_dat$F <- (exp_iv_dat$beta.exposure / exp_iv_dat$se.exposure)^2

  exp_iv_dat[exp_iv_dat$F <= 10, ] # F < 10 인 IV는 없는 것 확인.

  F_hba1c <- sum(exp_iv_dat$F / nrow(exp_iv_dat))
  F_hba1c

  exp_iv_dat

  ##########################################
  # 3. Outcome data preparation
  ##########################################
  # cad_gwas <- ""
  # ckd_gwas <- ""
  # ad_gwas <- ""

  out_dat <- read_outcome_data(filename=outcome_gwas,
                              snps = exp_dat$SNP,
                              sep = outcome_delim,
                              snp_col = outcome_snp,
                              beta_col = outcome_beta,
                              se_col = outcome_se,
                              eaf_col = outcome_eaf,
                              effect_allele_col = outcome_ea,
                              other_allele_col = outcome_oa,
                              pval_col = outcome_p,
                              chr_col = outcome_chr,
                              pos_col = outcome_pos
                            )


  out_dat$outcome <- paste0(outcome_name, "(", outcome_cohort, ")")

  if (outcome_type == "binary"){
    out_dat$ncase.outcome <- outcome_ncase
    out_dat$ncontrol.outcome <- outcome_ncontrol
  } else{
    out_dat$samplesize.outcome <- outcome_n
  }

  out_dat

  ##########################################
  # 4. Harmonization
  ##########################################

  exp_h_out_dat <- harmonise_data(exp_iv_dat, out_dat, action=1)
  exp_h_out_dat

  saveRDS(exp_h_out_dat, paste0(exposure_name, "_h_", outcome_name, "_GeneWindow", window, "kb.RDS"))
}


exp_h_out_dat <- readRDS(paste0(exposure_name, "_h_", outcome_name, "_GeneWindow", window, "kb.RDS"))

##########################################
# 5. Steiger test
##########################################
## Steiger test
if (exposure_type == "binary"){
  r.exposure <- get_r_from_lor(
      lor=exp_h_out_dat$beta.exposure,
      af=exp_h_out_dat$eaf.exposure,
      ncase=exp_h_out_dat$ncase.exposure,
      ncontrol=exp_h_out_dat$ncontrol.exposure,
      prevalence=exp_h_out_dat$ncase.exposure / (exp_h_out_dat$ncase.exposure + exp_h_out_dat$ncontrol.exposure),
      model = "logit",
      correction = FALSE
  )
}
if (outcome_type == "binary"){
  r.outcome <- get_r_from_lor(
      lor=exp_h_out_dat$beta.outcome,
      af=exp_h_out_dat$eaf.outcome,
      ncase=exp_h_out_dat$ncase.outcome,
      ncontrol=exp_h_out_dat$ncontrol.outcome,
      prevalence=exp_h_out_dat$ncase.outcome / (exp_h_out_dat$ncase.outcome + exp_h_out_dat$ncontrol.outcome),
      model = "logit",
      correction = FALSE
      )
  exp_h_out_dat$r.outcome <- r.outcome
}

out_steiger <- directionality_test(exp_h_out_dat)
out_steiger


write.table(out_steiger, paste0("steiger_", exposure_name, "_", outcome_name, "_GeneWindow", window, "kb.csv"), 
            sep=",", row.names=FALSE, quote=FALSE)
##########################################
# 6. TSMR
##########################################
res <- try(mr(exp_h_out_dat , method_list = c('mr_ivw',
                                              'mr_wald_ratio', 
                                              'mr_egger_regression',
                                              'mr_weighted_median',
                                              'mr_weighted_mode'
))) 

res <- generate_odds_ratios(res)
res

write.table(res, paste0("tsmr_", exposure_name, "_", outcome_name, "_GeneWindow", window, "kb.csv"), 
            sep=",", row.names=FALSE, quote=FALSE)

##########################################
# 7. Pleiotropy test
##########################################
mr_egger_intercept <- mr_pleiotropy_test(exp_h_out_dat)
mr_egger_intercept

write.table(mr_egger_intercept, paste0("mr_egger_intercept_", exposure_name, "_", outcome_name, "_GeneWindow", window, "kb.csv"), 
            sep=",", row.names=FALSE, quote=FALSE)

##########################################
# 8. Heterogeneity test
##########################################
hetero <- mr_heterogeneity(exp_h_out_dat, method_list = c("mr_egger_regression", "mr_ivw"))
hetero
write.table(hetero, paste0("hetero_", exposure_name, "_", outcome_name, "_GeneWindow", window, "kb.csv"), 
            sep=",", row.names=FALSE, quote=FALSE)
