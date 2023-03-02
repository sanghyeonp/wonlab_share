main <- function(ref_gwas, ref_gwas_delim, ref_gwas_snp, ref_gwas_beta,
                alt_gwas, alt_gwas_delim, alt_gwas_snp, alt_gwas_beta,
                alt_rev_beta,
                snplist_, ref_name="NA", alt_name="NA"){
  # snplist_ : a single string with , separated list of SNPs

  logs <- c()
  log_ <- paste0("Computing Cohen's kappa coefficient:\n", "\tReference: ", ref_gwas, "\n\tAlternative: ", alt_gwas); logs <- save_log(logs, log_)

  ### Read files, rename columns, and select only the necessary columns
  df_ref <- fread(ref_gwas, sep=ref_gwas_delim, header=TRUE, stringsAsFactors=FALSE)
  df_alt <- fread(alt_gwas, sep=alt_gwas_delim, header=TRUE, stringsAsFactors=FALSE)

  names(df_ref)[names(df_ref) == ref_gwas_snp] <- 'snp'
  names(df_ref)[names(df_ref) == ref_gwas_beta] <- 'beta'
  names(df_alt)[names(df_alt) == alt_gwas_snp] <- 'snp'
  names(df_alt)[names(df_alt) == alt_gwas_beta] <- 'beta'

  df_ref <- df_ref[, c('snp', 'beta')]
  df_alt <- df_alt[, c('snp', 'beta')]

  ### Select only the specified SNPs
  snplist_ <- strsplit(snplist_, split = ",")[[1]]

  df_ref <- df_ref[df_ref$snp %in% c(snplist_), ]
  df_alt <- df_alt[df_alt$snp %in% c(snplist_), ]
  log_ <- paste0("Number of SNPs in reference: ", nrow(df_ref)); logs <- save_log(logs, log_)
  log_ <- paste0("Number of SNPs in alternative: ", nrow(df_alt)); logs <- save_log(logs, log_)

  ### Reverse effect direction
  if (alt_rev_beta){
    df_alt$beta <- df_alt$beta * -1
  }

  ### Merge two data
  df_ <- merge(df_ref, df_alt, by = "snp", all = FALSE, suffixes=c('.ref', '.alt'))
  log_ <- paste0("Number of SNPs in both reference and alternative: ", nrow(df_)); logs <- save_log(logs, log_)

  ### Categorize based on effect direction
  df_ <- df_ %>% 
    mutate(beta.ref.c = if_else(beta.ref > 0, 1, 0),
          beta.alt.c = if_else(beta.alt > 0, 1, 0),
          concordance = if_else(beta.ref.c == beta.alt.c, TRUE, FALSE)
    )

  ### Run Cohen's kappa test
  cohen <- fmsb::Kappa.test(x=df_$beta.ref.c, y=df_$beta.alt.c, conf.level=0.95)

  ### Result as dataframe
  if (ref_name == "NA"){
    ref_name <- basename(ref_gwas)
  }
  if (alt_name == "NA"){
    alt_name <- basename(alt_gwas)
  }
  result <- data.frame (ref  = basename(ref_name),
                        alt = basename(alt_name),
                        cohenkappa = cohen$Result$estimate,
                        Z = cohen$Result$statistic,
                        P = pnorm(q=cohen$Result$statistic, lower.tail=FALSE) * 2,
                        CI_95L = cohen$Result$conf.int[1],
                        CI_95U = cohen$Result$conf.int[2],
                        judgement = cohen$Judgement
                        )

  return(list(cohen, logs, result))
}