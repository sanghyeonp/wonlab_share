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

    ### Exposure 
    parser$add_argument("--exp-gwas", dest="exposure_gwas", required=TRUE,
                        help="Path to the exposure GWAS summary statistics.")
    parser$add_argument("--exp-delim", dest="exposure_delim", required=TRUE,
                        help="Delimiter for the exposure GWAS. Options = ['tab', 'whitespace', 'comma'].")
    parser$add_argument("--exp-name", dest="exposure_name", required=TRUE,
                        help="Name of the exposure.")
    parser$add_argument("--exp-cohort", dest="exposure_cohort", required=TRUE,
                        help="Name of the exposure GWAS cohort.")
    parser$add_argument("--exp-type", dest="exposure_type", required=TRUE,
                        help="Type of the exposure. Choices = ['binary', 'quantitative']")
    parser$add_argument("--exp-snp", dest="exposure_snp", required=TRUE,
                        help="Column name for SNP.")
    parser$add_argument("--exp-chr", dest="exposure_chr", required=TRUE,
                        help="Column name for the chromosome.")
    parser$add_argument("--exp-pos", dest="exposure_pos", required=TRUE,
                        help="Column name for the base position.")
    parser$add_argument("--exp-ea", dest="exposure_ea", required=TRUE,
                        help="Column name for the effect allele.")
    parser$add_argument("--exp-oa", dest="exposure_oa", required=TRUE,
                        help="Column name for the other allele.")
    parser$add_argument("--exp-eaf", dest="exposure_eaf", required=TRUE,
                        help="Column name for the effect allele frequency.")
    parser$add_argument("--exp-beta", dest="exposure_beta", required=TRUE,
                        help="Column name for the beta.")
    parser$add_argument("--exp-se", dest="exposure_se", required=TRUE,
                        help="Column name for the standard error.")
    parser$add_argument("--exp-p", dest="exposure_p", required=TRUE,
                        help="Column name for the p-value.")
    parser$add_argument("--exp-ncase", dest="exposure_ncase", required=FALSE, default=0, type='integer',
                        help="The number of case if exposure is binary. Default = 0.")
    parser$add_argument("--exp-ncontrol", dest="exposure_ncontrol", required=FALSE, default=0, type='integer',
                        help="The number of control if exposure is binary. Default = 0.")
    parser$add_argument("--exp-n", dest="exposure_n", required=TRUE, type='integer',
                        help="The total number of samples.")

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

    ### Outcome
    parser$add_argument("--outcome-list-file", dest="outcome_list_file", required=TRUE, 
                        help="Specify outcome list file path.")
    parser$add_argument("--outcome-list-delim", dest="outcome_list_file_delim", required=TRUE, 
                        help="Specify outcome list file delimiter.")

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

    exposure_gwas <- args$exposure_gwas
    exposure_delim <- args$exposure_delim
    exposure_name <- args$exposure_name
    exposure_cohort <- args$exposure_cohort
    exposure_type <- args$exposure_type
    exposure_snp <- args$exposure_snp
    exposure_chr <- args$exposure_chr
    exposure_pos <- args$exposure_pos
    exposure_ea <- args$exposure_ea
    exposure_oa <- args$exposure_oa
    exposure_eaf <- args$exposure_eaf
    exposure_beta <- args$exposure_beta
    exposure_se <- args$exposure_se
    exposure_p <- args$exposure_p
    exposure_ncase <- args$exposure_ncase
    exposure_ncontrol <- args$exposure_ncontrol
    exposure_n <- args$exposure_n

    gene <- args$gene
    gene_chr <- args$gene_chr
    gene_start <- args$gene_start
    gene_end <- args$gene_end
    gene_cis_window <- args$gene_cis_window
    clump_r2 <- args$clump_r2
    clump_window <- args$clump_window
    clump_p <- args$clump_p

    F_thres <- args$F_thres

    outcome_list_file <- args$outcome_list_file
    outcome_list_file_delim <- args$outcome_list_file_delim

    n_cores <- args$n_cores

    verbose <- args$verbose

    outd <- args$outd

    #################################################################################################################################
    #########################
    ### Modify arguments ####
    #########################

    ## Delimiter setting
    delim_list <- list(tab = "\t", whitespace = " ", comma = ",")
    exposure_delim <- delim_list[[exposure_delim]]
    outcome_list_file_delim <- delim_list[[outcome_list_file_delim]]

    ## Output directory setting
    if (outd == "NA"){
        outd <- getwd()
    }

    #################################################################################################################################

    df_outcome_list <- fread(outcome_list_file, 
                            sep=outcome_list_file_delim
                            )
    
    df_outcome_list <- df_outcome_list %>%
                    mutate(dir_name = paste0(outd, "/", gsub(" ", "_", trait), ".", gsub(" ", "_", cohort))) %>%
                    mutate(prefix =  paste0(expsoure_name, ".", gsub(" ", "_", expsoure_cohort), ".", 
                                            trait, ".", gsub(" ", "_", cohort), ".",
                                            gene, "_GeneWindow", gene_cis_window, "kb")
                            )

    registerDoParallel(n_cores)

    foreach (idx=1:nrow(df_outcome_list), .combine=c) %dopar% {
        row <- df_outcome_list[idx, ]

        print(paste0("Running TSMR between ", exposure_name, " and ", row$trait, "..."))

        outcome_gwas <- row$gwas
        outcome_delim <- delim_list[[row$delim]]
        outcome_name <- row$trait
        outcome_cohort <- row$cohort
        outcome_type <- row$type
        outcome_snp <- row$snp_col
        outcome_chr <- row$chr_col
        outcome_pos <- row$pos_col
        outcome_ea <- row$ea_col
        outcome_oa <- row$oa_col
        outcome_eaf <- row$eaf_col
        outcome_beta <- row$beta_col
        outcome_se <- row$se_col
        outcome_p <- row$p_col
        outcome_ncase <- row$ncase
        outcome_ncontrol <- row$ncontrol
        outcome_n <- row$ntotal

        outd1 <- paste0(outd, "/", gsub(" ", "_", outcome_name), ".", gsub(" ", "_", outcome_cohort))
        
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
                            verbose, 
                            outd)
    }

    stopImplicitCluster()

    df_tsmr_combined <- NULL
    df_harmonized_combined <- NULL
    df_steiger_combined <- NULL
    df_mregger_intercept_combined <- NULL
    df_hetero_combined <- NULL

    for (idx in 1:nrow(df_outcome_list)){
        row <- df_outcome_list[idx, ]

        # Save TSMR results
        file <- paste0(row$dir_name, "/", "tsmr.", row$prefix, ".csv")
        if (file.exists(file)){
            temp <- fread(file)
            df_tsmr_combined <- rbind(df_tsmr_combined, temp)
        }

        # Save harmonized table
        file <- paste0(row$dir_name, "/", "harmonized.", row$prefix, ".csv")
        if (file.exists(file)){
            temp <- fread(file)
            df_harmonized_combined <- rbind(df_harmonized_combined, temp)
        }

        # Save Steiger test table
        file <- paste0(row$dir_name, "/", "steiger_test.", row$prefix, ".csv")
        if (file.exists(file)){
            temp <- fread(file)
            df_steiger_combined <- rbind(df_steiger_combined, temp)
        }

        # Save MR-Egger intercept test table
        file <- paste0(row$dir_name, "/", "mr_egger_intercept_test.", row$prefix, ".csv")
        if (file.exists(file)){
            temp <- fread(file)
            df_mregger_intercept_combined <- rbind(df_mregger_intercept_combined, temp)
        }

        # Save heterogeneity test table
        file <- paste0(row$dir_name, "/", "heterogeniety_test.", row$prefix, ".csv")
        if (file.exists(file)){
            temp <- fread(file)
            df_hetero_combined <- rbind(df_hetero_combined, temp)
        }

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
