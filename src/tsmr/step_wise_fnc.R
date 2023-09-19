library(data.table)
library(dplyr)


exp_dat_prep <- function(exposure_gwas, exposure_delim,
                        exposure_snp, exposure_chr, exposure_pos, exposure_ea, exposure_oa, exposure_eaf, exposure_beta, exposure_se, exposure_p,
                        prefix, outd){

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
        mutate(PVAL = as.numeric(PVAL))

    return (exp_dat)
}

exp_dat_gene_region_prep <- function(exp_dat, 
                                    gene_chr, gene_start, gene_cis_window, gene_end,
                                    prefix, outd){
    exp_dat <- exp_dat %>% 
        filter((CHR == gene_chr) & (POS > gene_start - (gene_cis_window * 1000)) & (POS <= gene_end + (gene_cis_window * 1000)))
    
    write.table(exp_dat, 
            paste0(outd, "/", "exposure_gene_region.", prefix, ".tsv"),
            sep="\t", row.names=FALSE, quote = FALSE
    )

    return (exp_dat)
}


exp_iv_dat_prep <- function(exp_dat,
                            clump_p, clump_r2, clump_window,
                            exposure_name, exposure_cohort, exposure_type, 
                            exposure_ncase, exposure_ncontrol, exposure_n,
                            verbose,
                            F_thres, perform_Fstat=FALSE){
    ## Perform clumping
    exp_dat <- filter(exp_dat, PVAL < clump_p)

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
                        clump_r2 =clump_r2,
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

    if (verbose){
        print(">>> Exposure clumping <<<")
        print(paste0("Number of IV after clumping: ", nrow(exp_iv_dat)))
        print(exp_iv_dat[, c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure")])
    }

    if (perform_Fstat){
        exp_iv_dat <- exp_iv_F_statistics(exp_iv_dat, F_thres, verbose)
    }

    return (exp_iv_dat)
}


exp_iv_F_statistics <- function(exp_iv_dat, F_thres, verbose){
    ##########################################
    # 2. IV validity check
    ##########################################

    ## F-statistics
    exp_iv_dat$F <- (exp_iv_dat$beta.exposure / exp_iv_dat$se.exposure)^2

    exp_iv_dat_F_false <- exp_iv_dat[exp_iv_dat$F <= F_thres, ] # F < 10 인 IV는 없는 것 확인.
    exp_iv_dat <- exp_iv_dat[exp_iv_dat$F > F_thres, ]
    
    F_hba1c <- sum(exp_iv_dat$F / nrow(exp_iv_dat))

    if (verbose){
        print(">>> F-statistics for IV <<<")
        print(paste0("Number of IV: ", nrow(exp_iv_dat)))
        print(paste0("Number of IV with F <= ", F_thres, " : ", nrow(exp_iv_dat_F_false)))
        print(paste0("Average F-statistics for IVs: ", F_hba1c))
        print(exp_iv_dat[, c("SNP", "F")])
    }

    return (exp_iv_dat)
}

##########################################
# 3. Outcome data preparation
##########################################
outcome_dat_prep <- function(exp_iv_dat,
                            outcome_gwas, outcome_delim, 
                            outcome_name, outcome_cohort, outcome_type,
                            outcome_snp, outcome_beta, outcome_se, outcome_eaf, outcome_ea, outcome_oa, outcome_p, outcome_chr, outcome_pos,
                            outcome_ncase, outcome_ncontrol, outcome_n                            
                            ){

    out_dat <- read_outcome_data(filename=outcome_gwas,
                                snps = exp_iv_dat$SNP,
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

    return (out_dat)
}


harmonization <- function(exp_iv_dat, out_dat,
                        prefix, outd, verbose){
    exp_h_out_dat <- harmonise_data(exp_iv_dat, out_dat, action=1)

    saveRDS(exp_h_out_dat, 
            paste0(outd, "/", "harmonized.", prefix, ".RDS")
            )
    write.table(exp_h_out_dat, 
                paste0(outd, "/", "harmonized.", prefix, ".csv"),
                sep=",", row.names=FALSE, quote=FALSE)

    if (verbose){
        print(">>> Harmonization <<<")
        print(paste0("Number of IV after harmonization: ", nrow(exp_h_out_dat)))
        print(exp_h_out_dat)
    }

    return (exp_h_out_dat)

}

perform_steiger_test <- function(exp_h_out_dat,
                                exposure_type, outcome_type,
                                prefix, outd, verbose){
    ## Steiger test
    if (exposure_type == "binary"){
        r.exposure <- get_r_from_lor(lor=exp_h_out_dat$beta.exposure,
                                    af=exp_h_out_dat$eaf.exposure,
                                    ncase=exp_h_out_dat$ncase.exposure,
                                    ncontrol=exp_h_out_dat$ncontrol.exposure,
                                    prevalence=exp_h_out_dat$ncase.exposure / (exp_h_out_dat$ncase.exposure + exp_h_out_dat$ncontrol.exposure),
                                    model = "logit",
                                    correction = FALSE
        )
    }
    if (outcome_type == "binary"){
        r.outcome <- get_r_from_lor(lor=exp_h_out_dat$beta.outcome,
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

    if (verbose){
        print(">>> Steiger test <<<")
        print(out_steiger)
    }

    write.table(out_steiger, 
                paste0(outd, "/", "steiger_test.", prefix, ".csv"),
                sep=",", row.names=FALSE, quote=FALSE)
    return (out_steiger)
}

perform_TSMR <- function(exp_h_out_dat, prefix, outd, verbose){
    res <- try(mr(exp_h_out_dat , method_list = c('mr_ivw',
                                            'mr_wald_ratio', 
                                            'mr_egger_regression',
                                            'mr_weighted_median',
                                            'mr_weighted_mode'
                ))) 

    res <- generate_odds_ratios(res)

    if (verbose){
        print(">>> Two-sample Mendelian randomization <<<")
        print(res)
    }

    write.table(res, 
                paste0(outd, "/", "tsmr.", prefix, ".csv"),
                sep=",", row.names=FALSE, quote=FALSE)

    return (res)
}


perform_pleiotropy_test <- function(exp_h_out_dat,
                                    prefix, outd, verbose){
    mr_egger_intercept <- mr_pleiotropy_test(exp_h_out_dat)

    if (verbose){
        print(">>> MR-Egger intercept test <<<")
        print(mr_egger_intercept)
    }

    write.table(mr_egger_intercept, 
                paste0(outd, "/", "mr_egger_intercept_test.", prefix, ".csv"),
                sep=",", row.names=FALSE, quote=FALSE)
    return (mr_egger_intercept)

}

perform_heterogneity_test <- function(exp_h_out_dat, prefix, outd, verbose){
    hetero <- mr_heterogeneity(exp_h_out_dat, method_list = c("mr_egger_regression", "mr_ivw"))

    if (verbose){
        print(">>> Heterogeneity test <<<")
        print(hetero)
    }

    write.table(hetero, 
                paste0(outd, "/", "heterogeniety_test.", prefix, ".csv"),
                sep=",", row.names=FALSE, quote=FALSE)

    return (hetero)
}


##################################################################
# Exposure의 gene body +- X kb 내에서 IV 추출.
##################################################################
run_single_gene_target_tsmr <- function(
    exposure_gwas, exposure_delim, 
    exposure_name, exposure_cohort, exposure_type, 
    exposure_snp, exposure_chr, exposure_pos, 
    exposure_ea, exposure_oa, exposure_eaf, 
    exposure_beta, exposure_se, exposure_p, 
    exposure_ncase, exposure_ncontrol, exposure_n, 
    outcome_gwas, outcome_delim, 
    outcome_name, outcome_cohort, outcome_type, 
    outcome_snp, outcome_chr, outcome_pos, 
    outcome_ea, outcome_oa, outcome_eaf, 
    outcome_beta, outcome_se, outcome_p, 
    outcome_ncase, outcome_ncontrol, outcome_n, 
    gene, gene_chr, gene_start, gene_end, gene_cis_window, 
    clump_r2, clump_window, clump_p, 
    F_thres, 
    verbose, 
    outd
    ){
    
    cat(paste0(":: TSMR analysis ::\n",
                "\tExposure: ", exposure_name, " (", exposure_cohort, ")\n",
                "\tOutcome: ", outcome_name, " (", outcome_cohort, ")\n")
        )
    
    prefix <- paste0(exposure_name, ".", gsub(" ", "_", exposure_cohort), ".", 
                    outcome_name, ".", gsub(" ", "_", outcome_cohort), ".",
                    gene, "_GeneWindow", gene_cis_window, "kb"
                    )

    if (!file.exists(paste0(outd, "/", "harmonized.", prefix, ".RDS"))){
        ##########################################
        # 1. Exposure data preparation
        ##########################################
        cat(">>> Preparing exposure <<<\n")
        exp_dat <- exp_dat_prep(exposure_gwas, exposure_delim,
                                exposure_snp, exposure_chr, exposure_pos, exposure_ea, exposure_oa, exposure_eaf, exposure_beta, exposure_se, exposure_p,
                                prefix, outd
                                )

        exp_dat <- exp_dat_gene_region_prep(exp_dat, 
                                            gene_chr, gene_start, gene_cis_window, gene_end,
                                            prefix, outd
                                            )

        ##########################################
        # 2. Exposure data clumping
        ##########################################
        cat(">>> Clumping <<<\n")
        exp_iv_dat <- exp_iv_dat_prep(exp_dat,
                                    clump_p, clump_r2, clump_window,
                                    exposure_name, exposure_cohort, exposure_type, 
                                    exposure_ncase, exposure_ncontrol, exposure_n,
                                    verbose,
                                    F_thres, perform_Fstat=TRUE
                                    )

        ##########################################
        # 3. Outcome data preparation
        ##########################################
        cat(">>> Preparing outcome <<<\n")
        out_dat <- outcome_dat_prep(exp_iv_dat,
                                    outcome_gwas, outcome_delim, 
                                    outcome_name, outcome_cohort, outcome_type,
                                    outcome_snp, outcome_beta, outcome_se, outcome_eaf, outcome_ea, outcome_oa, outcome_p, outcome_chr, outcome_pos,
                                    outcome_ncase, outcome_ncontrol, outcome_n

        )

        ##########################################
        # 4. Harmonization
        ##########################################
        cat(">>> Running harmonization <<<\n")
        exp_h_out_dat <- harmonization(exp_iv_dat, out_dat,
                                        prefix, outd, verbose)
    }

    ### Read in harmonized data if present
    exp_h_out_dat <- readRDS(paste0(outd, "/", "harmonized.", prefix, ".RDS"))

    ##########################################
    # 5. Steiger test
    ##########################################
    cat(">>> Running Steiger test <<<\n")
    out_steiger <- perform_steiger_test(exp_h_out_dat,
                                    exposure_type, outcome_type,
                                    prefix, outd, verbose)


    ##########################################
    # 6. TSMR
    ##########################################
    cat(">>> Running two-sample MR <<<\n")
    out_tsmr <- perform_TSMR(exp_h_out_dat, prefix, outd, verbose)

    ##########################################
    # 7. Pleiotropy test
    ##########################################
    cat(">>> Running pleiotropy test <<<\n")
    out_pleiotropy <- perform_pleiotropy_test(exp_h_out_dat,
                                            prefix, outd, verbose)


    ##########################################
    # 8. Heterogeneity test
    ##########################################
    cat(">>> Running heterogeneity test <<<\n")
    out_heterogeneity <- perform_heterogneity_test(exp_h_out_dat, prefix, outd, verbose)
    
    return (out_tsmr)
}


##################################################################
# :: IV list에서 다로 추출. ::
##################################################################

run_iv_specified_tsmr <- function(
    exposure_gwas, exposure_delim, 
    exposure_name, exposure_cohort, exposure_type, 
    exposure_snp, exposure_chr, exposure_pos, 
    exposure_ea, exposure_oa, exposure_eaf, 
    exposure_beta, exposure_se, exposure_p, 
    exposure_ncase, exposure_ncontrol, exposure_n,
    iv_list_path, iv_exp_p_thres,
    outcome_gwas, outcome_delim, 
    outcome_name, outcome_cohort, outcome_type, 
    outcome_snp, outcome_chr, outcome_pos, 
    outcome_ea, outcome_oa, outcome_eaf, 
    outcome_beta, outcome_se, outcome_p, 
    outcome_ncase, outcome_ncontrol, outcome_n, 
    F_thres, 
    verbose, 
    outd
    ){
    
    cat(paste0(":: TSMR analysis ::\n",
                "\tExposure: ", exposure_name, " (", exposure_cohort, ")\n",
                "\tOutcome: ", outcome_name, " (", outcome_cohort, ")\n")
        )
    
    prefix <- paste0(exposure_name, ".", gsub(" ", "_", exposure_cohort), ".", 
                    outcome_name, ".", gsub(" ", "_", outcome_cohort), "."
                    )

    if (!file.exists(paste0(outd, "/", "harmonized.", prefix, ".RDS"))){
        ##########################################
        # 1. Exposure data preparation
        ##########################################
        cat(">>> Preparing exposure <<<\n")
        exp_dat <- exp_dat_prep(exposure_gwas, exposure_delim,
                                exposure_snp, exposure_chr, exposure_pos, exposure_ea, exposure_oa, exposure_eaf, exposure_beta, exposure_se, exposure_p,
                                prefix, outd
                                )
        
        # Read IV list file
        iv_list <- fread(iv_list_path, header=FALSE)[[1]]
        exp_dat <- exp_dat %>%
                        filter(SNP %in% iv_list)
        
        df_iv_exclude_list <- NULL
        df_iv_missing <- NULL
        df_iv_p <- NULL

        if (length(iv_list) != nrow(exp_dat)){
            iv_missing <- setdiff(iv_list, exp_dat$SNP)
            n_iv_missing <- length(iv_missing)
            df_iv_missing <- data.frame(SNP = iv_missing,
                                        Exclusion = rep("Missing", n_iv_missing))
            df_iv_exclude_list <- bind_rows(df_iv_exclude_list, df_iv_missing)
        }

        if (nrow(exp_dat[exp_dat$PVAL > iv_exp_p_thres, ]) > 0){
            df_iv_p <- exp_dat[exp_dat$PVAL > iv_exp_p_thres, ]
            df_iv_p$Exclusion <- paste0("P-value > ", iv_exp_p_thres)
            df_iv_exclude_list <- bind_rows(df_iv_exclude_list, df_iv_p)

            exp_dat <- exp_dat[exp_dat$PVAL <= iv_exp_p_thres, ]
        }

        if (!is.null(df_iv_exclude_list)){
            df_iv_exclude_list <- df_iv_exclude_list %>%
                                    mutate(CHR = coalesce(CHR, -9),
                                            POS = coalesce(POS, -9),
                                            A1 = coalesce(A1, 'NA'),
                                            A2 = coalesce(A2, 'NA'),
                                            EAF = coalesce(EAF, -9),
                                            BETA = coalesce(BETA, -9),
                                            SE = coalesce(SE, -9),
                                            PVAL = coalesce(PVAL, -9)
                                            )
            write.table(df_iv_exclude_list, 
                        paste0(outd, "/", "iv_excluded.", prefix, ".txt"),
                        sep=" ", row.names=FALSE, quote=FALSE
                        )
        }

        if (nrow(exp_dat) == 0){
            stop()
        }

        exp_iv_dat <- format_data(dat=exp_dat,
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
                                pos_col = "POS")

    exp_iv_dat$exposure <-paste0(exposure_name, "(", exposure_cohort, ")")
    if (exposure_type == "binary"){
        exp_iv_dat$ncase.exposure <- exposure_ncase
        exp_iv_dat$ncontrol.exposure <- exposure_ncontrol
    } else{
        exp_iv_dat$samplesize.exposure <- exposure_n
    }

        ##########################################
        # 3. Outcome data preparation
        ##########################################
        cat(">>> Preparing outcome <<<\n")
        out_dat <- outcome_dat_prep(exp_iv_dat,
                                    outcome_gwas, outcome_delim, 
                                    outcome_name, outcome_cohort, outcome_type,
                                    outcome_snp, outcome_beta, outcome_se, outcome_eaf, outcome_ea, outcome_oa, outcome_p, outcome_chr, outcome_pos,
                                    outcome_ncase, outcome_ncontrol, outcome_n

        )

        ##########################################
        # 4. Harmonization
        ##########################################
        cat(">>> Running harmonization <<<\n")
        exp_h_out_dat <- harmonization(exp_iv_dat, out_dat,
                                        prefix, outd, verbose)
    }

    ### Read in harmonized data if present
    exp_h_out_dat <- readRDS(paste0(outd, "/", "harmonized.", prefix, ".RDS"))

    ##########################################
    # 5. Steiger test
    ##########################################
    cat(">>> Running Steiger test <<<\n")
    out_steiger <- perform_steiger_test(exp_h_out_dat,
                                    exposure_type, outcome_type,
                                    prefix, outd, verbose)


    ##########################################
    # 6. TSMR
    ##########################################
    cat(">>> Running two-sample MR <<<\n")
    out_tsmr <- perform_TSMR(exp_h_out_dat, prefix, outd, verbose)

    ##########################################
    # 7. Pleiotropy test
    ##########################################
    cat(">>> Running pleiotropy test <<<\n")
    out_pleiotropy <- perform_pleiotropy_test(exp_h_out_dat,
                                            prefix, outd, verbose)


    ##########################################
    # 8. Heterogeneity test
    ##########################################
    cat(">>> Running heterogeneity test <<<\n")
    out_heterogeneity <- perform_heterogneity_test(exp_h_out_dat, prefix, outd, verbose)
    
    return (out_tsmr)
}