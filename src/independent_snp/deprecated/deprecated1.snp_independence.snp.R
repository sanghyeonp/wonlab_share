library(data.table)
library(dplyr)
library(tidyr)
library(foreach)
library(doParallel)
require(argparse)
options(error=traceback)

#######################
### Notes
# GWAS 1 is the reference SNP
# Compute per chromosome
# Compare the given SNP list
#######################

#######################
# Parser arguments
#######################
### Define parser arguments
parser <- argparse::ArgumentParser(description=":: SNP independence by P-value and LD ::", 
                                formatter_class="argparse.ArgumentDefaultsHelpFormatter")

parser$add_argument("--trait1", required=TRUE, 
                    help="Specify the trait name for GWAS 1.")
parser$add_argument("--trait2", required=TRUE, 
                    help="Specify the trait name for GWAS 2.")

# GWAS1 SNP list
parser$add_argument("--gwas1-snplist", dest = "file_gwas1_snplist", required=TRUE, 
                    help="Path to the SNP list. Must have Chr, Pos, and SNP information.")
parser$add_argument("--delim-gwas1-snplist", dest="delim_gwas1_snplist",  required=TRUE,
                    help="Specify delimiter for GWAS summary statistics. Options = [tab, whitespace, comma]")
parser$add_argument("--snp-col-gwas1-snplist", dest="snp_gwas1_snplist",  required=TRUE,
                    help="Specify SNP column name.")
parser$add_argument("--secondary-snp-col-gwas1-snplist", dest="secondary_snp_gwas1_snplist",  required=TRUE,
                    help="Specify secondary SNP column name. The SNP should be separated by semicolon.")

# GWAS2 SNP list
parser$add_argument("--gwas2-snplist", dest = "file_gwas2_snplist", required=TRUE, 
                    help="Path to the SNP list. Must have Chr, Pos, and SNP information.")
parser$add_argument("--delim-gwas2-snplist", dest="delim_gwas2_snplist",  required=TRUE,
                    help="Specify delimiter for GWAS summary statistics. Options = [tab, whitespace, comma]")
parser$add_argument("--snp-col-gwas2-snplist", dest="snp_gwas2_snplist",  required=TRUE,
                    help="Specify SNP column name.")
parser$add_argument("--secondary-snp-col-gwas2-snplist", dest="secondary_snp_gwas2_snplist",  required=TRUE,
                    help="Specify secondary SNP column name. The SNP should be separated by semicolon.")

# GWAS1 summary statistics
parser$add_argument("--gwas1", dest = "file_gwas1", required=TRUE, 
                    help="Path to the 'other' GWAS summary statistics.")
parser$add_argument("--delim-gwas1", dest="delim_gwas1",  required=TRUE,
                    help="Specify delimiter for GWAS summary statistics. Options = [tab, whitespace, comma]")
parser$add_argument("--snp-col-gwas1", dest="snp_gwas1",  required=TRUE,
                    help="Specify SNP column name.")
parser$add_argument("--chr-col-gwas1", dest="chr_gwas1",  required=TRUE,
                    help="Specify chromosome column name.")
parser$add_argument("--pos-col-gwas1", dest="pos_gwas1",  required=TRUE,
                    help="Specify position column name.")
parser$add_argument("--pval-col-gwas1", dest="p_gwas1",  required=TRUE,
                    help="Specify P-value column name.")

# GWAS2 summary statistics
parser$add_argument("--gwas2", dest = "file_gwas2", required=TRUE, 
                    help="Path to the 'other' GWAS summary statistics.")
parser$add_argument("--delim-gwas2", dest="delim_gwas2",  required=TRUE,
                    help="Specify delimiter for GWAS summary statistics. Options = [tab, whitespace, comma]")
parser$add_argument("--snp-col-gwas2", dest="snp_gwas2",  required=TRUE,
                    help="Specify SNP column name.")
parser$add_argument("--chr-col-gwas2", dest="chr_gwas2",  required=TRUE,
                    help="Specify chromosome column name.")
parser$add_argument("--pos-col-gwas2", dest="pos_gwas2",  required=TRUE,
                    help="Specify position column name.")
parser$add_argument("--pval-col-gwas2", dest="p_gwas2",  required=TRUE,
                    help="Specify P-value column name.")

### Reference panel
parser$add_argument("--reference-panel", dest="reference_panel",  required=FALSE, default = "1kG",
                    help="Specify the reference panel to compute LD. Options = [1kg, UKB_random_10K]")

### Genome build
parser$add_argument("--genome-build", dest="genome_build",  required=FALSE, default = "GRCh37",
                    help="Specify the genome build. Currently available = [GRCh37]")

### Independence threshold
parser$add_argument("--r2-threshold", dest="r2_threshold",  required=TRUE, type = "numeric",
                    help="Specify R2 threshold to determine LD. If SNP from GWAS 1 has larger R2 than R2 threshold, then this SNP is considered as not independent from the signals from GWAS 2.")
parser$add_argument("--window",  required=TRUE, type = "integer",
                    help="Specify the one-sided window to search for SNPs from GWAS 2. Unit = kb.")
parser$add_argument("--pval-threshold", dest="p_threshold",  required=TRUE, type = "numeric",
                    help="Specify the P-value threshold to select SNPs from GWAS 2. SNPs from GWAS 2 below this threshold will be considered as candidate SNPs to compute LD.")

### Parallel processing
parser$add_argument("--n-thread", dest="n_thread",  required=FALSE, default = 1, type = "integer",
                    help="Specify the number of parallelization.")


### Output arguments
parser$add_argument("--dir-out", dest="dir_out",  required=FALSE, default = "NA",
                    help="Specify the directory for the output. If NA, then it will be saved in the current working directory.")
parser$add_argument("--prefix-out", dest="pref_out",  required=FALSE, default = "NA",
                    help="Specify the prefix of the output. If NA, then the prefix will be like 'snp_independence.<filename of gwas1 snplist>'.")
parser$add_argument("--delim-out", dest="delim_out",  required=FALSE, default = "comma",
                    help="Specify delimiter for the output file. Options = [tab, whitespace, comma]")


#### Read parser
args <- parser$parse_args()

trait1 <- args$trait1
trait2 <- args$trait2

file_gwas1_snplist <- args$file_gwas1_snplist
delim_gwas1_snplist <- args$delim_gwas1_snplist
snp_gwas1_snplist <- args$snp_gwas1_snplist
secondary_snp_gwas1_snplist <- args$secondary_snp_gwas1_snplist

file_gwas2_snplist <- args$file_gwas2_snplist
delim_gwas2_snplist <- args$delim_gwas2_snplist
snp_gwas2_snplist <- args$snp_gwas2_snplist
secondary_snp_gwas2_snplist <- args$secondary_snp_gwas2_snplist

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
genome_build <- args$genome_build
r2_threshold <- args$r2_threshold
window <- args$window
p_threshold <- args$p_threshold
n_thread <- args$n_thread
dir_out <- args$dir_out
pref_out <- args$pref_out
delim_out <- args$delim_out


#######################
# Utility
#######################
### PLINK software
plink_exe <- "/data1/sanghyeon/wonlab_contribute/combined/software/plink/plink"

### Reference panel for LD computation
# Reference panel splitted into chromosomes
reference_panel_list <- c("1kG" = "/data1/sanghyeon/wonlab_contribute/combined/data_common/reference_panel/OpenGWAS/EUR_CHR/EUR_chr",
                          "UKB_random_10K" = "/data1/sanghyeon/wonlab_contribute/combined/data_common/reference_panel/UKB_random10k_noQC/CHR/UKB_random10k_noqc_chr"
)

### Delimiter mapper
delim_mapper <- c("tab" = "\t",
                  "whitespace" = " ",
                  "comma" = ",")

### Extension mapper
exe_mapper <- c("tab" = "tsv",
                "whitespace" = "txt",
                "comma" = "csv")

### Parallel function
compute_snp_ld <- function(snp1, snp2, plink_exe, file_reference, chr){
  plink_out <- system(paste0(plink_exe, 
                             " --bfile ", file_reference, chr, 
                             " --ld ", snp1, " ", snp2),
                      intern = TRUE,
                      wait = TRUE,
                      ignore.stderr = TRUE)
  ld_result <- grep("R-sq", plink_out, value = TRUE)
  
  if (identical(ld_result, character(0))){
    r2 <- "NA"
  } else{
    ld_result <- ld_result[1] # Refer to important note 1
    if (nchar(ld_result) != 0){
      ld_result_query <- regmatches(ld_result, regexpr("R-sq\\s+=\\s+([0-9.]+([eE][-+]?[0-9]+)?)", ld_result))

      r2 <- as.character(as.numeric(strsplit(ld_result_query, " ")[[1]][length(strsplit(ld_result_query, " ")[[1]])]))
      
    } else{
      r2 <- "NA"
    }
  }

  return (data.frame(ref_snp = snp1,
                    other_snp = snp2,
                    R2 = r2))
}


##### << delete
# file_gwas1_snplist <- "/data1/sanghyeon/Projects/MetabolicSyndrome/MetS_2022_08/nat_genet.revision_1/reviewer_2/R2Q1/SNP_independence.snp_based.secondary_sig/02_query_cojo_from_clumping/MetS/MetS.clumped.cojo_snp_query"
# delim_gwas1_snplist <- "tab"
# snp_gwas1_snplist <- "SNP"
# secondary_snp_gwas1 <- "clumped_SNP"

# file_gwas2_snplist <- "/data1/sanghyeon/Projects/MetabolicSyndrome/MetS_2022_08/nat_genet.revision_1/reviewer_2/R2Q1/SNP_independence.snp_based.secondary_sig/02_query_cojo_from_clumping/FG/FG.clumped.cojo_snp_query"
# delim_gwas2_snplist <- "tab"
# snp_gwas2_snplist <- "SNP"
# secondary_snp_gwas2 <- "clumped_SNP"

# file_gwas1 <- "/data1/sanghyeon/Projects/MetabolicSyndrome/MetS_2022_08/results/8_GWAS/MetS12/QC.MetS12_GWAS_MetS.withQsnp.csv"
# delim_gwas1 <- "comma"
# snp_gwas1 <- "SNP"
# chr_gwas1 <- "CHR"
# pos_gwas1 <- "BP"
# p_gwas1 <- "Pval_Estimate"

# file_gwas2 <- "/data1/sanghyeon/Projects/MetabolicSyndrome/GWASsumstat/QC_final/CLEANED.FG_MAGIC.txt"
# delim_gwas2 <- "tab"
# snp_gwas2 <- "SNP"
# chr_gwas2 <- "CHR"
# pos_gwas2 <- "POS"
# p_gwas2 <- "Pval"

# reference_panel <- "1kG"
# genome_build <- 37
# r2_threshold <- 0.01
# window <- 500
# p_threshold <- 5e-8
# n_thread <- 10
##### >> delete

#######################
# Read data
#######################

df_gwas1_snplist <- fread(file_gwas1_snplist, 
                          sep = delim_mapper[[delim_gwas1_snplist]],
                          data.table = F,
                          nThread = n_thread,
                          showProgress = F)

df_gwas1_snplist <- df_gwas1_snplist %>%
    dplyr::select(!!as.name(snp_gwas1_snplist), !!as.name(secondary_snp_gwas1_snplist)) %>%
    tidyr::separate_rows(clumped_SNP, sep=";") %>%
    rename(leadSNP = !!as.name(snp_gwas1_snplist),
           secondarySNP = !!as.name(secondary_snp_gwas1_snplist))
# Make sure leadSNP itself is in the clumped SNP list
gwas1.leadSNP <- unique(df_gwas1_snplist$leadSNP)
gwas1.leadSNP_to_add <- gwas1.leadSNP[!(gwas1.leadSNP %in% df_gwas1_snplist$secondarySNP)]
if (length(gwas1.leadSNP_to_add) > 0){
    df_gwas1_snplist <- rbind(df_gwas1_snplist,
                              data.frame(leadSNP = gwas1.leadSNP_to_add,
                                         secondarySNP = gwas1.leadSNP_to_add))
}
# Map genomic position and P-value
df_gwas1 <- fread(file_gwas1, 
                  sep = delim_mapper[[delim_gwas1]],
                  data.table = F,
                  nThread = n_thread,
                  showProgress = F)
df_gwas1 <- df_gwas1 %>%
    dplyr::select(!!as.name(snp_gwas1), !!as.name(chr_gwas1), !!as.name(pos_gwas1), !!as.name(p_gwas1)) %>%
    rename(CHR = !!as.name(chr_gwas1),
           POS = !!as.name(pos_gwas1),
           P = !!as.name(p_gwas1))

df_gwas1_snplist <- merge(df_gwas1_snplist, df_gwas1, by.x = "secondarySNP", by.y = snp_gwas1, all.x = T)

df_gwas1_leadsnplist <- df_gwas1_snplist %>%
    dplyr::select(leadSNP) %>%
    distinct()
df_gwas1_leadsnplist <- merge(df_gwas1_leadsnplist, df_gwas1, by.x = "leadSNP", by.y = snp_gwas1, all.x = T)

## For GWAS2
df_gwas2_snplist <- fread(file_gwas2_snplist, 
                          sep = delim_mapper[[delim_gwas2_snplist]],
                          data.table = F,
                          nThread = n_thread,
                          showProgress = F)

df_gwas2_snplist <- df_gwas2_snplist %>%
    dplyr::select(!!as.name(snp_gwas2_snplist), !!as.name(secondary_snp_gwas2_snplist)) %>%
    tidyr::separate_rows(clumped_SNP, sep=";") %>%
    rename(leadSNP = !!as.name(snp_gwas2_snplist),
           secondarySNP = !!as.name(secondary_snp_gwas2_snplist))
# Make sure leadSNP itself is in the clumped SNP list
gwas2.leadSNP <- unique(df_gwas2_snplist$leadSNP)
gwas2.leadSNP_to_add <- gwas2.leadSNP[!(gwas2.leadSNP %in% df_gwas2_snplist$secondarySNP)]
if (length(gwas2.leadSNP_to_add) > 0){
    df_gwas2_snplist <- rbind(df_gwas2_snplist,
                              data.frame(leadSNP = gwas2.leadSNP_to_add,
                                         secondarySNP = gwas2.leadSNP_to_add))
}
# Map genomic position and P-value
df_gwas2 <- fread(file_gwas2, 
                  sep = delim_mapper[[delim_gwas2]],
                  data.table = F,
                  nThread = n_thread,
                  showProgress = F)
df_gwas2 <- df_gwas2 %>%
    dplyr::select(!!as.name(snp_gwas2), !!as.name(chr_gwas2), !!as.name(pos_gwas2), !!as.name(p_gwas2)) %>%
    rename(CHR = !!as.name(chr_gwas2),
           POS = !!as.name(pos_gwas2),
           P = !!as.name(p_gwas2))

df_gwas2_snplist <- merge(df_gwas2_snplist, df_gwas2, by.x = "secondarySNP", by.y = snp_gwas2, all.x = T)

df_gwas2_leadsnplist <- df_gwas2_snplist %>%
    dplyr::select(leadSNP) %>%
    distinct()
df_gwas2_leadsnplist <- merge(df_gwas2_leadsnplist, df_gwas2, by.x = "leadSNP", by.y = snp_gwas2, all.x = T)


#######################
# LD
#######################

file_reference_panel <- reference_panel_list[[reference_panel]]


quote <- F
if (delim_out == "comma"){
  quote <- T
}

if (dir_out == "NA"){
  dir_out <- "."
}
if (substr(dir_out, nchar(dir_out), nchar(dir_out)) != "/") {
  dir_out <- paste0(dir_out, "/")
}

#######################
# LD - GWAS1 lead SNP
#######################

df_ld_result.gwas1 <- NULL

for (chr in 1:22){
  print(paste0("Processing chromosome ", chr, "..."))
  df_gwas1_snplist.chr <- filter(df_gwas1_snplist, CHR == chr)
  df_gwas2_snplist.chr <- filter(df_gwas2_snplist, CHR == chr)
  
  if (nrow(df_gwas1_snplist.chr) != 0){
    # df_ld_per_chr comprises the columns below:
    # SNP_gwas1, 
    # SNP_gwas2_tested, SNP_gwas2_tested_pval, SNP_gwas2_tested_R2, 
    # In_LD ## In_LD means R2 of SNP_gwas2 > R2 threshold, hence SNP_gwas2 is in LD with SNP_gwas1, and SNP_gwas1 is not independent from GWAS2 genetic signal. 
    #       ## (So lower R2 threshold to conduct more strict LD comparison)
    # SNP_gwas2_in_LD, SNP_gwas2_in_LD_pval, SNP_gwas2_in_LD_R2
    df_ld_per_chr <- NULL
    unique_leadsnp <- unique(df_gwas1_snplist.chr$leadSNP)
    for (idx in 1:length(unique_leadsnp)){

      # row.df_gwas1_snplist.chr <- df_gwas1_snplist.chr[idx, ]
      ref_snp <- unique_leadsnp[idx]
      ref_snp.pos <- df_gwas1_leadsnplist[df_gwas1_leadsnplist$leadSNP == ref_snp, ]$POS
      
      df_gwas2_snplist.chr.filter <- filter(df_gwas2_snplist.chr, 
                                                  (POS >= ref_snp.pos - (window * 1000)) &
                                                  (POS < ref_snp.pos + (window * 1000)),
                                                  (P < p_threshold))

      df_gwas2_snplist.chr.filter  <- df_gwas2_snplist.chr.filter  %>%
                              distinct(secondarySNP, .keep_all = TRUE) # Important note 2

      if(nrow(df_gwas2_snplist.chr.filter ) != 0){
        other_snplist <- df_gwas2_snplist.chr.filter$secondarySNP
        
        snp_gwas2_tested <- c()
        snp_gwas2_tested_pval <- c()
        snp_gwas2_tested_r2 <- c()
        is_in_LD <- NA
        snp_gwas2_in_ld <- "NA"
        snp_gwas2_in_ld_pval <- "NA"
        snp_gwas2_in_ld_r2 <- "NA"
      
        # Register a parallel backend
        cl <- makeCluster(n_thread)
        registerDoParallel(cl)

        plink_ld_result <- foreach(other_snp = other_snplist,
                                .export = c("ref_snp", "plink_exe", "file_reference_panel", "chr"),
                                .packages = c("dplyr", "data.table"),
                                .combine = 'bind_rows') %dopar% {
                            compute_snp_ld(snp1 = ref_snp, 
                                          snp2 = other_snp, 
                                          plink_exe = plink_exe, 
                                          file_reference = file_reference_panel, 
                                          chr = chr)
                          }
        stopCluster(cl)
        
        snp_gwas2_tested <- plink_ld_result$other_snp
        df_tmp <- df_gwas2_snplist.chr.filter  %>%
                    dplyr::select(secondarySNP, P) %>%
                    filter(secondarySNP %in% snp_gwas2_tested)
        df_tmp <- merge(df_tmp, plink_ld_result %>% dplyr::select(-ref_snp), by.x = "secondarySNP", by.y = "other_snp", all.x = T)
        # df_tmp[[snp_gwas2_snplist]] <- factor(df_tmp[[snp_gwas2_snplist]], levels = snp_gwas2_tested)
        snp_gwas2_tested_pval <- df_tmp$P
        snp_gwas2_tested_r2 <- as.numeric(df_tmp$R2)

        idx_in_ld <- which(snp_gwas2_tested_r2 > r2_threshold)

        if (identical(idx_in_ld, integer(0))){ # No SNPs from GWAS2 in LD with SNP from GWAS1
          is_in_LD <- FALSE
          snp_gwas2_in_ld_leadsnp <- "NA"
        } else{
          is_in_LD <- TRUE
          snp_gwas2_in_ld <- paste(snp_gwas2_tested[idx_in_ld], collapse = ";")
          snp_gwas2_in_ld_pval <- paste(snp_gwas2_tested_pval[idx_in_ld], collapse = ";")
          snp_gwas2_in_ld_r2 <- paste(snp_gwas2_tested_r2[idx_in_ld], collapse = ";")
          snp_gwas2_in_ld_leadsnp <- paste(unique((df_gwas2_snplist.chr %>% 
                                               filter(secondarySNP %in% snp_gwas2_tested[idx_in_ld]))$leadSNP), collapse = ";")
        }

        snp_gwas2_tested <- paste(snp_gwas2_tested, collapse = ";")
        snp_gwas2_tested_pval <- paste(snp_gwas2_tested_pval, collapse = ";")
        snp_gwas2_tested_r2 <- paste(snp_gwas2_tested_r2, collapse = ";")
        

        df_ld_per_chr <- rbind(df_ld_per_chr,
                        data.frame(SNP_gwas1 = ref_snp,
                            SNP_gwas2_tested = snp_gwas2_tested,
                            SNP_gwas2_tested_pval = snp_gwas2_tested_pval,
                            SNP_gwas2_tested_R2 = snp_gwas2_tested_r2,
                            In_LD = is_in_LD,
                            SNP_gwas2_in_LD = snp_gwas2_in_ld,
                            SNP_gwas2_in_LD_pval = snp_gwas2_in_ld_pval,
                            SNP_gwas2_in_LD_R2 = snp_gwas2_in_ld_r2,
                            snp_gwas2_in_ld_leadsnp = snp_gwas2_in_ld_leadsnp)
                      )
      } else{
        df_ld_per_chr <- rbind(df_ld_per_chr,
                        data.frame(SNP_gwas1 = ref_snp,
                            SNP_gwas2_tested = "NA",
                            SNP_gwas2_tested_pval = "NA",
                            SNP_gwas2_tested_R2 = "NA",
                            In_LD = FALSE,
                            SNP_gwas2_in_LD = "NA",
                            SNP_gwas2_in_LD_pval = "NA",
                            SNP_gwas2_in_LD_R2 = "NA",
                            snp_gwas2_in_ld_leadsnp = "NA")
                      )
      }
    }

    df_ld_result.gwas1 <- rbind(df_ld_result.gwas1, df_ld_per_chr)
  }
}

write.table(df_ld_result.gwas1,
            file = paste0(dir_out, pref_out, ".compare_", trait1, ".to_", trait2, ".", exe_mapper[delim_out]),
            sep = delim_mapper[delim_out],
            row.names = F,
            quote = F)


#######################
# LD - GWAS1 lead SNP
#######################

df_ld_result.gwas2 <- NULL

for (chr in 1:22){
  print(paste0("Processing chromosome ", chr, "..."))
  df_gwas1_snplist.chr <- filter(df_gwas1_snplist, CHR == chr)
  df_gwas2_snplist.chr <- filter(df_gwas2_snplist, CHR == chr)
  
  if (nrow(df_gwas2_snplist.chr) != 0){
    # df_ld_per_chr comprises the columns below:
    # SNP_gwas1, 
    # SNP_gwas2_tested, SNP_gwas2_tested_pval, SNP_gwas2_tested_R2, 
    # In_LD ## In_LD means R2 of SNP_gwas2 > R2 threshold, hence SNP_gwas2 is in LD with SNP_gwas1, and SNP_gwas1 is not independent from GWAS2 genetic signal. 
    #       ## (So lower R2 threshold to conduct more strict LD comparison)
    # SNP_gwas2_in_LD, SNP_gwas2_in_LD_pval, SNP_gwas2_in_LD_R2
    df_ld_per_chr <- NULL
    unique_leadsnp <- unique(df_gwas2_snplist.chr$leadSNP)
    for (idx in 1:length(unique_leadsnp)){

      # row.df_gwas1_snplist.chr <- df_gwas1_snplist.chr[idx, ]
      ref_snp <- unique_leadsnp[idx]
      ref_snp.pos <- df_gwas2_leadsnplist[df_gwas2_leadsnplist$leadSNP == ref_snp, ]$POS
      
      df_gwas1_snplist.chr.filter <- filter(df_gwas1_snplist.chr, 
                                                  (POS >= ref_snp.pos - (window * 1000)) &
                                                  (POS < ref_snp.pos + (window * 1000)),
                                                  (P < p_threshold))

      df_gwas1_snplist.chr.filter  <- df_gwas1_snplist.chr.filter  %>%
                              distinct(secondarySNP, .keep_all = TRUE) # Important note 2

      if(nrow(df_gwas1_snplist.chr.filter ) != 0){
        other_snplist <- df_gwas1_snplist.chr.filter$secondarySNP
        
        snp_gwas1_tested <- c()
        snp_gwas1_tested_pval <- c()
        snp_gwas1_tested_r2 <- c()
        is_in_LD <- NA
        snp_gwas1_in_ld <- "NA"
        snp_gwas1_in_ld_pval <- "NA"
        snp_gwas1_in_ld_r2 <- "NA"
      
        # Register a parallel backend
        cl <- makeCluster(n_thread)
        registerDoParallel(cl)

        plink_ld_result <- foreach(other_snp = other_snplist,
                                .export = c("ref_snp", "plink_exe", "file_reference_panel", "chr"),
                                .packages = c("dplyr", "data.table"),
                                .combine = 'bind_rows') %dopar% {
                            compute_snp_ld(snp1 = ref_snp, 
                                          snp2 = other_snp, 
                                          plink_exe = plink_exe, 
                                          file_reference = file_reference_panel, 
                                          chr = chr)
                          }
        stopCluster(cl)
        
        snp_gwas1_tested <- plink_ld_result$other_snp
        df_tmp <- df_gwas1_snplist.chr.filter  %>%
                    dplyr::select(secondarySNP, P) %>%
                    filter(secondarySNP %in% snp_gwas1_tested)
        df_tmp <- merge(df_tmp, plink_ld_result %>% dplyr::select(-ref_snp), by.x = "secondarySNP", by.y = "other_snp", all.x = T)
        # df_tmp[[snp_gwas1_snplist]] <- factor(df_tmp[[snp_gwas2_snplist]], levels = snp_gwas2_tested)
        snp_gwas1_tested_pval <- df_tmp$P
        snp_gwas1_tested_r2 <- as.numeric(df_tmp$R2)

        idx_in_ld <- which(snp_gwas1_tested_r2 > r2_threshold)

        if (identical(idx_in_ld, integer(0))){ # No SNPs from GWAS1 in LD with SNP from GWAS2
          is_in_LD <- FALSE
        } else{
          is_in_LD <- TRUE
          snp_gwas1_in_ld <- paste(snp_gwas1_tested[idx_in_ld], collapse = ";")
          snp_gwas1_in_ld_pval <- paste(snp_gwas1_tested_pval[idx_in_ld], collapse = ";")
          snp_gwas1_in_ld_r2 <- paste(snp_gwas1_tested_r2[idx_in_ld], collapse = ";")
          snp_gwas1_in_ld_leadsnp <- paste(unique((df_gwas1_snplist.chr %>% 
                                               filter(secondarySNP %in% snp_gwas1_tested[idx_in_ld]))$leadSNP), collapse = ";")
        }

        snp_gwas1_tested <- paste(snp_gwas1_tested, collapse = ";")
        snp_gwas1_tested_pval <- paste(snp_gwas1_tested_pval, collapse = ";")
        snp_gwas1_tested_r2 <- paste(snp_gwas1_tested_r2, collapse = ";")
        

        df_ld_per_chr <- rbind(df_ld_per_chr,
                        data.frame(SNP_gwas1 = ref_snp,
                            SNP_gwas2_tested = snp_gwas1_tested,
                            SNP_gwas2_tested_pval = snp_gwas1_tested_pval,
                            SNP_gwas2_tested_R2 = snp_gwas1_tested_r2,
                            In_LD = is_in_LD,
                            SNP_gwas2_in_LD = snp_gwas1_in_ld,
                            SNP_gwas2_in_LD_pval = snp_gwas1_in_ld_pval,
                            SNP_gwas2_in_LD_R2 = snp_gwas1_in_ld_r2,
                            snp_gwas2_in_ld_leadsnp = snp_gwas1_in_ld_leadsnp)
                      )
      } else{
        df_ld_per_chr <- rbind(df_ld_per_chr,
                        data.frame(SNP_gwas1 = ref_snp,
                            SNP_gwas2_tested = "NA",
                            SNP_gwas2_tested_pval = "NA",
                            SNP_gwas2_tested_R2 = "NA",
                            In_LD = FALSE,
                            SNP_gwas2_in_LD = "NA",
                            SNP_gwas2_in_LD_pval = "NA",
                            SNP_gwas2_in_LD_R2 = "NA",
                            snp_gwas2_in_ld_leadsnp = "NA")
                      )
      }
    }

    df_ld_result.gwas2 <- rbind(df_ld_result.gwas2, df_ld_per_chr)
  }
}

write.table(df_ld_result.gwas2,
            file = paste0(dir_out, pref_out, ".compare_", trait2, ".to_", trait1, ".", exe_mapper[delim_out]),
            sep = delim_mapper[delim_out],
            row.names = F,
            quote = F)

system("rm plink.log plink.nosex",
      wait = TRUE)


