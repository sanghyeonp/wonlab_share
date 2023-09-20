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

fnc_file <- paste(getCurrentFileLocation(), "step_wise_fnc.R", sep="/")
source(fnc_file)

##################################################################################################################
# Two-sample Mendelian randomization analysis for a target gene
# - Exposure: Drug-targeting gene region within the biomarker (endophenotype) GWAS
#   - Drug가 target 하는 유전자 영역을 exposure GWAS에서 filter 후, clumping 진행하여 instrument variables 추출.
# - Outcome: Any quantitative or binary trait

# Sanghyeon Park
# 2023.09.15
##################################################################################################################

tryCatch({
    ### Load packages from env_R.R
    list.of.packages <- c("argparse", "data.table", "rstudioapi", "dplyr", "foreach", "doParallel")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

    require(argparse)

    library(data.table)
    library(dplyr)
    library(TwoSampleMR)
    library(foreach)
    library(doParallel)

    #################################################################################################################################
    ###############################
    ### Define parser arguments ###
    ###############################
    parser <- argparse::ArgumentParser(description=":: Run two-sample Mendelian randomization for a target gene ::", 
                                    formatter_class="argparse.ArgumentDefaultsHelpFormatter")

    ### Exposure list
    parser$add_argument("--exposure-list-file", dest="exposure_list_file", required=TRUE, 
                        help="Specify exposure list file path.")
    parser$add_argument("--exposure-list-delim", dest="exposure_list_file_delim", required=TRUE, 
                        help="Specify exposure list file delimiter.")

    ### Outcome list
    parser$add_argument("--outcome-list-file", dest="outcome_list_file", required=TRUE, 
                        help="Specify outcome list file path.")
    parser$add_argument("--outcome-list-delim", dest="outcome_list_file_delim", required=TRUE, 
                        help="Specify outcome list file delimiter.")

    ### Exposure filter
    ## Filter by gene region
    parser$add_argument("--gene", dest="gene", required=TRUE,
                        help="Name of the gene of interest.")
    parser$add_argument("--gene-chr", dest="gene_chr", required=TRUE, type='integer',
                        help="Chromosome number of the specified gene.")
    parser$add_argument("--gene-start", dest="gene_start", required=TRUE, type='integer',
                        help="Specify the gene start position.")
    parser$add_argument("--gene-end", dest="gene_end", required=TRUE, type='integer',
                        help="Specify the gene end position.")
    parser$add_argument("--gene-cis-window", dest="gene_cis_window", required=TRUE, type='integer',
                        help="Specify the cis-window around the gene body in kb. For example, `--gene-cis-window 500` means 500 kb around the gene body.")

    ### Clumping
    parser$add_argument("--clump-r2", dest="clump_r2", required=TRUE, type='double',
                        help="Specify the r2 threshold for clumping.")
    parser$add_argument("--clump-window", dest="clump_window", required=TRUE, type='integer',
                        help="Specify the window size for clumping.")
    parser$add_argument("--clump-p", dest="clump_p", required=TRUE, type='double',
                        help="Specify the P-value threshold for clumping.")

    ### F-statistics
    parser$add_argument("--F-thres", dest="F_thres", required=FALSE, type='double', default=10,
                        help="Specify the F-statistics threshold to filter IV. Default = 10.")

    ### Causal effect reverse-code
    parser$add_argument("--reverse-effect", dest="reverse_effect", action="store_true",
                        help="Specify to reverse-code the result of TSMR.")

    ### Parallelization
    parser$add_argument("--n-cores", dest="n_cores", required=FALSE, default=1, type='integer',
                        help="Specify the number of cores.")

    ### Others
    parser$add_argument("--verbose", action="store_true", 
                        help="Print extra output.")

    ### Save option
    parser$add_argument("--outd", required=FALSE, default="NA",
                        help="Specify the output directory path. Default = `current working directory`.")

    #################################################################################################################################
    #########################
    ### Read in arguments ###
    #########################
    args <- parser$parse_args()

    gene <- args$gene
    gene_chr <- args$gene_chr
    gene_start <- args$gene_start
    gene_end <- args$gene_end
    gene_cis_window <- args$gene_cis_window
    clump_r2 <- args$clump_r2
    clump_window <- args$clump_window
    clump_p <- args$clump_p

    F_thres <- args$F_thres

    exposure_list_file <- args$exposure_list_file
    exposure_list_file_delim <- args$exposure_list_file_delim

    outcome_list_file <- args$outcome_list_file
    outcome_list_file_delim <- args$outcome_list_file_delim

    reverse_effect <- args$reverse_effect

    n_cores <- args$n_cores

    verbose <- args$verbose

    outd <- args$outd

    #################################################################################################################################
    #########################
    ### Modify arguments ####
    #########################

    ## Delimiter setting
    delim_list <- list(tab = "\t", whitespace = " ", comma = ",")
    exposure_list_file_delim <- delim_list[[exposure_list_file_delim]]
    outcome_list_file_delim <- delim_list[[outcome_list_file_delim]]

    ## Output directory setting
    if (outd == "NA"){
        outd <- getwd()
    }

    #################################################################################################################################

    df_exposure_list <- fread(exposure_list_file,
                            sep=exposure_list_file_delim
                            )
    df_exposure_list <- df_exposure_list %>% 
                        rename('exp.trait' = 'trait',
                            'exp.cohort' = 'cohort',
                            'exp.type' = 'type',
                            'exp.gwas' = 'gwas',
                            'exp.delim' = 'delim',
                            'exp.snp_col' = 'snp_col',
                            'exp.chr_col' = 'chr_col',
                            'exp.pos_col' = 'pos_col',
                            'exp.ea_col' = 'ea_col',
                            'exp.oa_col' = 'oa_col',
                            'exp.eaf_col' = 'eaf_col',
                            'exp.beta_col' = 'beta_col',
                            'exp.se_col' = 'se_col',
                            'exp.p_col' = 'p_col',
                            'exp.ncase' = 'ncase',
                            'exp.ncontrol' = 'ncontrol',
                            'exp.ntotal' = 'ntotal'
                            )

    df_outcome_list <- fread(outcome_list_file, 
                            sep=outcome_list_file_delim
                            )
    df_outcome_list <- df_outcome_list %>% 
                        rename('op.trait' = 'trait',
                            'op.cohort' = 'cohort',
                            'op.type' = 'type',
                            'op.gwas' = 'gwas',
                            'op.delim' = 'delim',
                            'op.snp_col' = 'snp_col',
                            'op.chr_col' = 'chr_col',
                            'op.pos_col' = 'pos_col',
                            'op.ea_col' = 'ea_col',
                            'op.oa_col' = 'oa_col',
                            'op.eaf_col' = 'eaf_col',
                            'op.beta_col' = 'beta_col',
                            'op.se_col' = 'se_col',
                            'op.p_col' = 'p_col',
                            'op.ncase' = 'ncase',
                            'op.ncontrol' = 'ncontrol',
                            'op.ntotal' = 'ntotal'
                            )

    df_combined <- NULL
    for (idx1 in 1:nrow(df_exposure_list)){
        row1 <- df_exposure_list[idx1, ]
        
        for (idx2 in 1:nrow(df_outcome_list)){
            row2 <- df_outcome_list[idx2, ]
            new_row <- cbind(row1, row2)
            df_combined <- rbind(df_combined, new_row)
        }
    }

    df_combined <- df_combined %>%
                mutate(dir_name = paste0(outd, "/", gsub(" ", "_", exp.trait), ".", gsub(" ", "_", exp.cohort),
                            ".", gsub(" ", "_", op.trait), ".", gsub(" ", "_", op.cohort)
                        )) %>%
                mutate(prefix = paste0(exp.trait, ".", gsub(" ", "_", exp.cohort), ".", 
                    op.trait, ".", gsub(" ", "_", op.cohort), ".",
                    gene, "_GeneWindow", gene_cis_window, "kb"
                    ))

    print(df_combined)

    registerDoParallel(n_cores)

    foreach (idx=1:nrow(df_combined), .combine=c) %dopar% {
        row <- df_combined[idx, ]

        exposure_gwas <- row$exp.gwas
        exposure_delim <- delim_list[[row$exp.delim]]
        exposure_name <- row$exp.trait
        exposure_cohort <- row$exp.cohort
        exposure_type <- row$exp.type
        exposure_snp <- row$exp.snp_col
        exposure_chr <- row$exp.chr_col
        exposure_pos <- row$exp.pos_col
        exposure_ea <- row$exp.ea_col
        exposure_oa <- row$exp.oa_col
        exposure_eaf <- row$exp.eaf_col
        exposure_beta <- row$exp.beta_col
        exposure_se <- row$exp.se_col
        exposure_p <- row$exp.p_col
        exposure_ncase <- row$exp.ncase
        exposure_ncontrol <- row$exp.ncontrol
        exposure_n <- row$exp.ntotal
    
        outcome_gwas <- row$op.gwas
        outcome_delim <- delim_list[[row$op.delim]]
        outcome_name <- row$op.trait
        outcome_cohort <- row$op.cohort
        outcome_type <- row$op.type
        outcome_snp <- row$op.snp_col
        outcome_chr <- row$op.chr_col
        outcome_pos <- row$op.pos_col
        outcome_ea <- row$op.ea_col
        outcome_oa <- row$op.oa_col
        outcome_eaf <- row$op.eaf_col
        outcome_beta <- row$op.beta_col
        outcome_se <- row$op.se_col
        outcome_p <- row$op.p_col
        outcome_ncase <- row$op.ncase
        outcome_ncontrol <- row$op.ncontrol
        outcome_n <- row$op.ntotal

        outd1 <- paste0(outd, "/", gsub(" ", "_", exposure_name), ".", gsub(" ", "_", exposure_cohort),
                            ".", gsub(" ", "_", outcome_name), ".", gsub(" ", "_", outcome_cohort)
                        )

        if (!dir.exists(outd1)) {
            dir.create(outd1)
        }

        run_single_gene_target_tsmr(exposure_gwas, exposure_delim, 
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
                                reverse_effect,
                                verbose, 
                                outd1)

    }

    stopImplicitCluster()

    df_tsmr_combined <- NULL
    df_harmonized_combined <- NULL
    df_steiger_combined <- NULL
    df_mregger_intercept_combined <- NULL
    df_hetero_combined <- NULL

    for (idx in 1:nrow(df_combined)){
        row <- df_combined[idx, ]

        # Save TSMR results
        file <- paste0(row$dir_name, "/", "tsmr.", row$prefix, ".csv")

        if (file.exists(file)){
            temp <- tryCatch({fread(file)
                                }, error = function(e) {
                                temp <- NULL
                                })
            if (!is.null(temp)){
                df_tsmr_combined <- rbind(df_tsmr_combined, temp)
            }
        }

        # Save harmonized table
        file <- paste0(row$dir_name, "/", "harmonized.", row$prefix, ".csv")
        if (file.exists(file)){
            temp <- tryCatch({fread(file)
                                }, error = function(e) {
                                temp <- NULL
                                })
            if (!is.null(temp)){
                df_harmonized_combined <- rbind(df_harmonized_combined, temp)
            }
        }

        # Save Steiger test table
        file <- paste0(row$dir_name, "/", "steiger_test.", row$prefix, ".csv")
        if (file.exists(file)){
            temp <- tryCatch({fread(file)
                                }, error = function(e) {
                                temp <- NULL
                                })
            if (!is.null(temp)){
                df_steiger_combined <- rbind(df_steiger_combined, temp)
            }
        }

        # Save MR-Egger intercept test table
        file <- paste0(row$dir_name, "/", "mr_egger_intercept_test.", row$prefix, ".csv")
        if (file.exists(file)){
            temp <- tryCatch({fread(file)
                                }, error = function(e) {
                                temp <- NULL
                                })
            if (!is.null(temp)){
                df_mregger_intercept_combined <- rbind(df_mregger_intercept_combined, temp)
            }
        }

        # Save heterogeneity test table
        file <- paste0(row$dir_name, "/", "heterogeniety_test.", row$prefix, ".csv")
        if (file.exists(file)){
            temp <- tryCatch({fread(file)
                                }, error = function(e) {
                                temp <- NULL
                                })
            if (!is.null(temp)){
                df_hetero_combined <- rbind(df_hetero_combined, temp)
            }
        }

    }

    write.table(df_tsmr_combined, 
                paste0(outd, "/", "tsmr.combined.csv"),
                sep=",", row.names=FALSE, quote=FALSE
                )
    
    write.table(df_harmonized_combined, 
            paste0(outd, "/", "harmonized.combined.csv"),
            sep=",", row.names=FALSE, quote=FALSE
            )

    write.table(df_steiger_combined, 
            paste0(outd, "/", "steiger_test.combined.csv"),
            sep=",", row.names=FALSE, quote=FALSE
            )

    write.table(df_mregger_intercept_combined, 
            paste0(outd, "/", "mr_egger_intercept_test.combined.csv"),
            sep=",", row.names=FALSE, quote=FALSE
            )

    write.table(df_hetero_combined, 
            paste0(outd, "/", "heterogeneity.combined.csv"),
            sep=",", row.names=FALSE, quote=FALSE
            )

    combine_results(save_ivw_mregger_only=TRUE,
                    save_wr_only=TRUE)

}, error = function(e){
    # Print the error message
    cat("Error: ", conditionMessage(e), "\n")

    # Get the line number where the error occurred
    traceback_lines <- traceback()
    if (!is.null(traceback_lines)) {
        # Extract the line number from the traceback
        line_number <- as.numeric(gsub("[^0-9]", "", traceback_lines[[1]]))
        cat("Error occurred on line:", line_number, "\n")
    } else {
        cat("Line number not available in traceback.\n")
    }
})
