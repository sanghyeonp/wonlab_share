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
    list.of.packages <- c("argparse", "data.table", "rstudioapi", "dplyr")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

    require(argparse)

    library(data.table)
    library(dplyr)
    library(TwoSampleMR)

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

    ### Causal effect reverse-code
    parser$add_argument("--reverse-effect", dest="reverse_effect", action="store_true",
                        help="Specify to reverse-code the result of TSMR.")

    ### Outcome
    parser$add_argument("--op-gwas", dest="outcome_gwas", required=TRUE,
                        help="Path to the outcome GWAS summary statistics.")
    parser$add_argument("--op-delim", dest="outcome_delim", required=TRUE,
                        help="Delimiter for the outcome GWAS. Options = ['tab', 'whitespace', 'comma'].")
    parser$add_argument("--op-name", dest="outcome_name", required=TRUE,
                        help="Name of the outcome.")
    parser$add_argument("--op-cohort", dest="outcome_cohort", required=TRUE,
                        help="Name of the outcome GWAS cohort.")
    parser$add_argument("--op-type", dest="outcome_type", required=TRUE,
                        help="Type of the outcome. Choices = ['binary', 'quantitative']")
    parser$add_argument("--op-snp", dest="outcome_snp", required=TRUE,
                        help="Column name for SNP.")
    parser$add_argument("--op-chr", dest="outcome_chr", required=TRUE,
                        help="Column name for the chromosome.")
    parser$add_argument("--op-pos", dest="outcome_pos", required=TRUE,
                        help="Column name for the base position.")
    parser$add_argument("--op-ea", dest="outcome_ea", required=TRUE,
                        help="Column name for the effect allele.")
    parser$add_argument("--op-oa", dest="outcome_oa", required=TRUE,
                        help="Column name for the other allele.")
    parser$add_argument("--op-eaf", dest="outcome_eaf", required=TRUE,
                        help="Column name for the effect allele frequency.")
    parser$add_argument("--op-beta", dest="outcome_beta", required=TRUE,
                        help="Column name for the beta.")
    parser$add_argument("--op-se", dest="outcome_se", required=TRUE,
                        help="Column name for the standard error.")
    parser$add_argument("--op-p", dest="outcome_p", required=TRUE,
                        help="Column name for the p-value.")
    parser$add_argument("--op-ncase", dest="outcome_ncase", required=FALSE, default=0, type='integer',
                        help="The number of case if outcome is binary. Default = 0.")
    parser$add_argument("--op-ncontrol", dest="outcome_ncontrol", required=FALSE, default=0, type='integer',
                        help="The number of control if outcome is binary. Default = 0.")
    parser$add_argument("--op-n", dest="outcome_n", required=TRUE, type='integer',
                        help="The total number of samples.")

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

    reverse_effect <- args$reverse_effect

    outcome_gwas <- args$outcome_gwas
    outcome_delim <- args$outcome_delim
    outcome_name <- args$outcome_name
    outcome_cohort <- args$outcome_cohort
    outcome_type <- args$outcome_type
    outcome_snp <- args$outcome_snp
    outcome_chr <- args$outcome_chr
    outcome_pos <- args$outcome_pos
    outcome_ea <- args$outcome_ea
    outcome_oa <- args$outcome_oa
    outcome_eaf <- args$outcome_eaf
    outcome_beta <- args$outcome_beta
    outcome_se <- args$outcome_se
    outcome_p <- args$outcome_p
    outcome_ncase <- args$outcome_ncase
    outcome_ncontrol <- args$outcome_ncontrol
    outcome_n <- args$outcome_n

    verbose <- args$verbose

    outd <- args$outd

    #################################################################################################################################
    #########################
    ### Modify arguments ####
    #########################

    ## Delimiter setting
    delim_list <- list(tab = "\t", whitespace = " ", comma = ",")
    exposure_delim <- delim_list[[exposure_delim]]
    outcome_delim <- delim_list[[outcome_delim]]

    ## Output directory setting
    if (outd == "NA"){
        outd <- getwd()
    } 
    outd <- paste0(outd, "/", gsub(" ", "_", outcome_name), ".", gsub(" ", "_", outcome_cohort))
    if (!dir.exists(outd)) {
        dir.create(outd)
    }

    #################################################################################################################################

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
                                outd)

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
