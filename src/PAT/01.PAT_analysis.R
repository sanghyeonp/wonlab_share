# 1. Discovery GWAS에서 lead SNP 추리기.
# **Observed number of transferable loci 구하기.**
# 2. 추린 lead SNP 기준으로 proxy SNPs을 찾아서, lead SNP 별로 credible set 만들기.
# - R2 >= 0.8 그리고 window 50kb를 잡고, 그 안에 포함된 SNP 중 P-value < 100 * PleadSNP  
# 3. Credible set 별로, 해당 locus가 transferable 하다 안하다를 판단하고, credible set 별로 하나의 SNP 추리기.
# - Credible set에 포함되는 variants 중, 하나라도 다른 인종 GWAS에서 P-value < 0.05 이고, effect direction이 matching 할 경우. -> 여기서 transferable loci라고 분류된 것은 “# of observed transferable loci”
# - (Report 하는 용도) 만약 하나의 locus (i.e., credible set)에서 여러 variant가 위 경우에 만족한다면, 그 중 제일 p-value가 significant한 variant 하나를 기준으로 잡음.
# **Expected number of transferable loci 구하기**
# 4. Lead SNP 마다 alpha = 0.05로 잡고, power를 계산.
# 5. 계산된 power들을 다 summation하고 rounding하면 “expected # of loci to be significantly associated 랑 같음.”
# **PAT calculation**
# 6. PAT = Observed loci / Expected loci


#============================================================
# Load required libraries and utility functions
#============================================================
library(data.table)
library(dplyr)
library(foreach)
library(doParallel)

source("util_PAT.R")


#============================================================
# 0. Input files and column definitions
#============================================================

## Discovery GWAS (Trait 1)
gwas.trait1 <- "/data1/sanghyeon/Projects/MetS_gSEM_EAS/src_EUR/gSEM/06.Qsnp/v1m1/usergwas.v1.m1.MetS_EUR.assoc.Qsnp.tsv"

col_snp.trait1    <- "SNP"
col_chr.trait1    <- "CHR"
col_pos.trait1    <- "BP"
col_a1.trait1     <- "A1"
col_a1freq.trait1 <- "MAF"
col_b.trait1      <- "est"
col_se.trait1     <- "SE"
col_p.trait1      <- "Pval_Estimate"


## Validation GWAS (Trait 2)
gwas.trait2 <- "/data1/sanghyeon/Projects/MetS_gSEM_EAS/src_EAS/gSEM/06.Qsnp/v1m2/usergwas.v1.m2.MetS_EAS.assoc.Qsnp.tsv"

col_snp.trait2    <- "SNP"
col_chr.trait2    <- "CHR"
col_pos.trait2    <- "BP"
col_a1.trait2     <- "A1"
col_a1freq.trait2 <- "MAF"
col_b.trait2      <- "est"
col_se.trait2     <- "SE"
col_p.trait2      <- "Pval_Estimate"

## Lead SNP file from COJO
f.leadsnp <- "/data1/sanghyeon/Projects/MetS_gSEM_EAS/src_EUR/COJO/v1m1.MetS_EUR/cojo_out_jma.v1m1.MetS_EUR.csv"
col_snp.leadsnp <- "SNP"
col_p.leadsnp   <- "p"


## Reference panel and PLINK settings
reference_panel_chr <- "/data1/sanghyeon/wonlab_contribute/combined/data_common/reference_panel/UKB_random10k_noQC/CHR/UKB_random10k_noqc_chr"
plink_exe <- "/data1/sanghyeon/wonlab_contribute/combined/software/plink/plink"

credible_window <- 50
credible_r2     <- 0.8

n_thread <- 20

## Power calculation
alpha <- 0.05
N_trait2 <- 1334596

is_binary_trait2 <- FALSE
# If the trait 2 is binary, specify case and control. Otherwise, set it NA.
Ncase_trait2 <- NA
Ncontrol_trait2 <- NA


#============================================================
# Setup log
log_file <- "PAT_analysis.log"
log_con  <- file(log_file, open = "a")
sink(log_con, type = "message")
on.exit({
    sink(type = "message")
    close(log_con)
}, add = TRUE)

#============================================================
# 1. Load and preprocess discovery GWAS summary statistics
#============================================================

df.gwas.trait1 <- fread(gwas.trait1, data.table = FALSE) %>%
    dplyr::select(
        all_of(c(
            col_snp.trait1,
            col_chr.trait1,
            col_pos.trait1,
            col_a1.trait1,
            col_a1freq.trait1,
            col_b.trait1,
            col_se.trait1,
            col_p.trait1
        ))
    ) %>%
    rename(
        SNP          = all_of(col_snp.trait1),
        CHR          = all_of(col_chr.trait1),
        POS          = all_of(col_pos.trait1),
        a1_trait1    = all_of(col_a1.trait1),
        freq_trait1  = all_of(col_a1freq.trait1),
        beta_trait1  = all_of(col_b.trait1),
        se_trait1    = all_of(col_se.trait1),
        p_trait1     = all_of(col_p.trait1)
    )


#============================================================
# 2. Extract lead SNPs from discovery GWAS
#============================================================

df.leadsnp <- fread(f.leadsnp, data.table = FALSE) %>%
    dplyr::select(all_of(c(col_snp.leadsnp, col_p.leadsnp))) %>%
    rename(
        SNP      = all_of(col_snp.leadsnp),
        p_trait1 = all_of(col_p.leadsnp)
    ) %>%
    mutate(credible_set = seq_len(nrow(.)))


#============================================================
# 3. Construct credible sets using LD information
#============================================================

## Register parallel backend
cl <- makeCluster(n_thread)
registerDoParallel(cl)

df.credible <- foreach(
    credible_idx = df.leadsnp$credible_set,
    .export = c(
        "df.leadsnp",
        "df.gwas.trait1",
        "credible_window",
        "credible_r2",
        "plink_exe",
        "reference_panel_chr",
        "parse_proxy",
        "find_credible_set"
    ),
    .packages = c("dplyr", "data.table"),
    .combine  = "bind_rows"
) %dopar% {

    find_credible_set(
        credible_idx    = credible_idx,
        df.leadsnp      = df.leadsnp,
        df_gwas         = df.gwas.trait1,
        credible_window = credible_window,
        credible_r2     = credible_r2,
        plink_exe       = plink_exe,
        reference_panel = reference_panel_chr
    )
}

stopCluster(cl)

write.table(
    df.credible,
    "df_credible.csv",
    sep = ",",
    row.names = FALSE,
    quote = FALSE
)

saveRDS(df.credible, "df_credible.rds")

#============================================================
# 4. Merge validation GWAS information
#============================================================

df.gwas.trait2 <- fread(gwas.trait2, sep = "\t", data.table = FALSE) %>%
    dplyr::select(
        all_of(c(
            col_snp.trait2,
            col_chr.trait2,
            col_pos.trait2,
            col_a1.trait2,
            col_a1freq.trait2,
            col_b.trait2,
            col_se.trait2,
            col_p.trait2
        ))
    ) %>%
    rename(
        SNP          = all_of(col_snp.trait2),
        CHR_trait2 = all_of(col_chr.trait2),
        POS_trait2 = all_of(col_pos.trait2),
        a1_trait2    = all_of(col_a1.trait2),
        freq_trait2  = all_of(col_a1freq.trait2),
        beta_trait2  = all_of(col_b.trait2),
        se_trait2    = all_of(col_se.trait2),
        p_trait2     = all_of(col_p.trait2)
    )

df.credible2 <- merge(df.credible, df.gwas.trait2, by = "SNP", all.x = TRUE) %>%
    mutate(
        a1_trait2_new   = ifelse(a1_trait1 == a1_trait2, a1_trait2, a1_trait1),
        freq_trait2_new = ifelse(a1_trait1 == a1_trait2, freq_trait2, 1 - freq_trait2),
        beta_trait2_new = ifelse(a1_trait1 == a1_trait2, beta_trait2, -beta_trait2)
    ) %>%
    dplyr::select(-a1_trait2, -freq_trait2, -beta_trait2) %>%
    rename(
        a1_trait2   = a1_trait2_new,
        freq_trait2 = freq_trait2_new,
        beta_trait2 = beta_trait2_new
    )

write.table(
    df.credible2,
    "df_credible.merged_aligned.csv",
    sep = ",",
    row.names = FALSE,
    quote = FALSE
)

df.credible3 <- na.omit(df.credible2)

write.table(
    df.credible3,
    "df_credible.merged_aligned.noNA.csv",
    sep = ",",
    row.names = FALSE,
    quote = FALSE
)


#============================================================
# 5. Select transferable variants per credible set
#============================================================

df.credible.final <- df.credible3 %>%
    mutate(
        transferable.per_snp = p_trait2 < 0.05 & ((beta_trait1 < 0 & beta_trait2 < 0) | (beta_trait1 > 0 & beta_trait2 > 0))
    )


######### Power calculation: Code from Huang
if (is_binary_trait2){
    pow <- power_for_binary_Huang(alpha, df.credible.final, Ncase_trait2, Ncontrol_trait2)
} else{
    pow <- power_for_quan_Huang(alpha, df.credible.final, N_trait2)
}

df.credible.final$power <- pow
df.credible.final$powered <- pow >= 0.8

write.table(
    df.credible.final,
    "transferability_and_power_per_snp.csv",
    sep = ",",
    row.names = FALSE,
    quote = FALSE
)
################

### Define transferable loci per credible set
df.transferable_loci <- data.frame()
for (cs in unique(df.credible.final$credible_set)){
    df.cs <- df.credible.final %>%
        filter(credible_set == cs)
    
    leadsnp <- df.leadsnp %>%
        filter(credible_set == cs) %>%
        pull(SNP)
    df.gwas.trait1_lead <- df.gwas.trait1 %>%
        filter(SNP == leadsnp)
    pos <- df.gwas.trait1_lead$POS

    df.val_cis <- df.gwas.trait2 %>%
        filter(CHR_trait2 == df.cs$CHR[1] &
               POS_trait2 >= pos - credible_window * 1000 &
               POS_trait2 <= pos + credible_window * 1000) %>%
        filter(p_trait2 < 1e-3)

    if (any(df.cs$transferable.per_snp)){
        transferability <- "transferable"
    } else if (any(df.cs$powered) & !any(df.cs$p_trait2 < 0.05) & nrow(df.val_cis) == 0){
        transferability <- "non-transferable"
    } else{
        transferability <- "underpowered"
    }
    df.transferable_loci <- rbind(df.transferable_loci,
                                data.frame(credible_set = cs,
                                            transferability = transferability))
}

message("Number of loci analyzed: ", nrow(df.transferable_loci))
message("Number of transferable loci: ", nrow(df.transferable_loci %>% filter(transferability == "transferable")))
message("Number of non-transferable loci: ", nrow(df.transferable_loci %>% filter(transferability == "non-transferable")))
message("Number of underpowered loci: ", nrow(df.transferable_loci %>% filter(transferability == "underpowered")))

##
df.transferable_loci.final <- df.leadsnp %>% dplyr::select(-p_trait1) %>%
    left_join(df.transferable_loci, by="credible_set") %>%
    left_join(df.credible.final %>% filter(SNP_type=="lead") %>% dplyr::select(credible_set, power, powered), by="credible_set") %>%
    left_join(df.gwas.trait1, by="SNP")

write.table(df.transferable_loci.final, "transferability_per_locus.csv", sep = ",", row.names = FALSE, quote = FALSE)

########### Power-adjusted transferability (PAT) calculation
power_sum <- sum(df.transferable_loci.final$power, na.rm = TRUE) #sum up power and divide by number of loci/variants for total expected power
message("Number of expected transferable loci: ", power_sum)
PAT_ratio <- nrow(df.transferable_loci.final %>% filter(transferability == "transferable"))/power_sum
message("PAT ratio: ", PAT_ratio)
