#------------------------------------------------------------
# parse_proxy
#------------------------------------------------------------
# This function identifies proxy SNPs in linkage disequilibrium
# with a given sentinel SNP using PLINK.
#
# The function calls PLINK to calculate pairwise r2 values
# between the target SNP and nearby variants within a
# user-defined genomic window.
#
# Arguments:
#   snp1            Character string that specifies the sentinel SNP ID
#   plink_exe       Path to the PLINK executable
#   reference_panel Prefix of the PLINK reference panel files
#   chr             Chromosome identifier appended to the reference panel
#   window          Window size in kilobases for LD calculation
#   r2              Minimum r2 threshold for proxy SNP selection
#
# Output:
#   PLINK output files with the prefix equal to the sentinel SNP ID
#------------------------------------------------------------
parse_proxy <- function(snp1,
                        plink_exe,
                        reference_panel,
                        chr,
                        window,
                        r2) {
    
    # Construct the PLINK command for LD-based proxy SNP extraction
    cmd <- paste0(
        plink_exe,
        " --bfile ", reference_panel, chr,
        " --r2",
        " --ld-snp ", snp1,
        " --ld-window-kb ", window,
        " --ld-window-r2 ", r2,
        " --out ", snp1
    )
    
    # Execute the PLINK command and wait until completion
    system(cmd, wait = TRUE)
}


find_credible_set <- function(credible_idx,
                              df.leadsnp,
                              df_gwas,
                              credible_window,
                              credible_r2,
                              plink_exe, reference_panel
){
    leadsnp <- df.leadsnp[credible_idx, ]$SNP
    chr <- df_gwas[df_gwas$SNP == leadsnp, ]$CHR
    pos <- df_gwas[df_gwas$SNP == leadsnp, ]$POS
    pval <- df_gwas[df_gwas$SNP == leadsnp, ]$p_trait1

    df_gwas.window <- df_gwas %>%
        # p_trait1 < 5e-100 to avoid
        filter(CHR == chr &
               POS >= pos - credible_window * 1000 &
               POS <= pos + credible_window * 1000 &
                   (p_trait1 < pval * 100 | p_trait1 < 5e-100))

    if (nrow(df_gwas.window) > 0){
        snplist <- df_gwas.window$SNP
        df.r2 <- NULL
        parse_proxy(snp1 = leadsnp,
                     plink_exe = plink_exe,
                     reference_panel = reference_panel,
                     chr = chr,
                     window = credible_window,
                     r2 = credible_r2)
        df.r2 <- fread(paste0(leadsnp, ".ld"), data.table=F)

        df.r2 <- df.r2 %>%
            mutate(leadSNP=SNP_A) %>%
            filter(!is.na(R2)) %>%
            mutate(other_snp=SNP_B) %>%
            dplyr::select(leadSNP, other_snp, R2)
        
        system(paste0("rm ", leadsnp, ".ld ", leadsnp, ".log"), wait=TRUE)
        if (nrow(df.r2) > 0){
            df_gwas.window.r2 <- merge(df_gwas.window, df.r2,
                                    by.x = "SNP", by.y = "other_snp", all.x = T) %>%
                # Make sure only R2 > credible_r2 is retained
                filter(R2 >= credible_r2) %>%
                mutate(credible_set = credible_idx)
            df_gwas.window.r2$SNP_type <- ifelse(df_gwas.window.r2$SNP == leadsnp, "lead", "proxy")
            return (df_gwas.window.r2)
        } else{
            return (data.frame())
        }
    } else{
        return (data.frame())
    }

}
#------------------------------------------------------------
# find_credible_set
#------------------------------------------------------------
# This function constructs a credible set for a given lead SNP
# based on physical proximity and linkage disequilibrium criteria.
#
# For a selected lead SNP, the function first defines a genomic
# window and applies a p-value filter relative to the lead SNP.
# It then identifies proxy SNPs that meet an r2 threshold using PLINK
# and returns variants that satisfy both positional and LD constraints.
#
# Arguments:
#   credible_idx     Integer index that identifies the lead SNP
#                     in df.leadsnp
#   df.leadsnp       Data frame that contains lead SNP information
#   df_gwas          GWAS summary statistics data frame
#   credible_window  Window size in kilobases around the lead SNP
#   credible_r2      Minimum r2 threshold for credible set inclusion
#   plink_exe        Path to the PLINK executable
#   reference_panel  Prefix of the PLINK reference panel files
#
# Output:
#   A data frame that contains SNPs in the credible set with
#   LD information and SNP type annotation. An empty data frame
#   is returned if no SNP satisfies the criteria.
#------------------------------------------------------------
find_credible_set <- function(credible_idx,
                            df.leadsnp,
                            df_gwas,
                            credible_window,
                            credible_r2,
                            plink_exe,
                            reference_panel) {
    
    # Extract the lead SNP identifier
    leadsnp <- df.leadsnp[credible_idx, ]$SNP
    
    # Retrieve chromosome, position, and p-value of the lead SNP
    chr  <- df_gwas[df_gwas$SNP == leadsnp, ]$CHR
    pos  <- df_gwas[df_gwas$SNP == leadsnp, ]$POS
    pval <- df_gwas[df_gwas$SNP == leadsnp, ]$p_trait1
    
    # Subset GWAS variants within the genomic window
    # and apply a relative p-value threshold
    df_gwas.window <- df_gwas %>%
        filter(
            CHR == chr,
            POS >= pos - credible_window * 1000,
            POS <= pos + credible_window * 1000,
            p_trait1 < pval * 100
        )
    
    # Exit early if no variants are present in the window
    if (nrow(df_gwas.window) == 0) {
        return(data.frame())
    }
    
    # Run PLINK to calculate LD between the lead SNP and nearby variants
    parse_proxy(
        snp1            = leadsnp,
        plink_exe       = plink_exe,
        reference_panel = reference_panel,
        chr             = chr,
        window          = credible_window,
        r2              = credible_r2
    )
    
    # Read PLINK LD output
    df.r2 <- fread(paste0(leadsnp, ".ld"), data.table = FALSE)
    
    # Process LD results and retain relevant columns
    df.r2 <- df.r2 %>%
        mutate(leadSNP = SNP_A,
            other_snp = SNP_B) %>%
        filter(!is.na(R2)) %>%
        dplyr::select(leadSNP, other_snp, R2)
    
    # Remove temporary PLINK output files
    system(paste0("rm ", leadsnp, ".ld ", leadsnp, ".log ", leadsnp, ".nosex"), wait = TRUE)
    
    # Exit if no LD pairs satisfy the criteria
    if (nrow(df.r2) == 0) {
        return(data.frame())
    }
    
    # Merge LD information with GWAS variants
    df_gwas.window.r2 <- merge(
        df_gwas.window,
        df.r2,
        by.x = "SNP",
        by.y = "other_snp",
        all.x = TRUE
    ) %>%
        filter(R2 >= credible_r2) %>%
        mutate(
            credible_set = credible_idx,
            SNP_type = ifelse(SNP == leadsnp, "lead", "proxy")
        )
    
    return(df_gwas.window.r2)
}


power_for_quan_Huang <- function(alpha, df.credible.final, N_trait2){
    ## URL: https://github.com/Nsallah1/GH_Manuscript/blob/main/power_calc_quant.R
    th <- qchisq(alpha, df=1, lower.tail=F) # Significance threshold for chi-square, corresponding to P-value threshold
    beta <- df.credible.final$beta_trait1 #beta in europeans
    maf <- df.credible.final$freq_trait2 # frequency in south asians
    n <- N_trait2

    pow <- matrix(NA, nrow=length(beta), ncol=length(maf)); # Pre-allocalte matrix to store power values
    rownames(pow) <- beta; colnames(pow) <- maf; # Assign row and column names

    q2 <- 2*maf*(1-maf)*(beta^2); ncp <- n*q2/(1-q2); # Calculate qsq and then NCP parameter
    pow <- pchisq(th, df=1, lower.tail=F, ncp=ncp)
    return(pow)
}

power_for_binary_Huang <- function(alpha, df.credible.final, Ncase_trait2, Ncontrol_trait2){
    ## URL: https://github.com/Nsallah1/GH_Manuscript/blob/main/power_calc_binary.R
    n <- Ncase_trait2+Ncontrol_trait2
    # ratio cases controls
    phi <- Ncase_trait2/n
    f <- df.credible.final$freq_trait2
    # beta europeans
    b <- df.credible.final$beta_trait1

    pow <- pchisq(qchisq(alpha,df=1,lower = F), df=1, ncp = 2*f*(1-f)*n*phi*(1-phi)*b^2, lower = F)
    return(pow)
}