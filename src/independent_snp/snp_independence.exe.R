library(data.table)
library(dplyr)
library(foreach)
library(doParallel)
require(argparse)

#######################
### Notes
# GWAS 1 is the reference SNP
# Compute per chromosome
#######################

#######################
# Parser arguments
#######################
### Define parser arguments
parser <- argparse::ArgumentParser(description=":: SNP independence by P-value and LD ::", 
                                formatter_class="argparse.ArgumentDefaultsHelpFormatter")

# GWAS1 SNP list
parser$add_argument("--gwas1-snplist", dest = "file_gwas1_snplist", required=TRUE, 
                    help="Path to the SNP list. Must have Chr, Pos, and SNP information.")
parser$add_argument("--delim-gwas1-snplist", dest="delim_gwas1_snplist",  required=TRUE,
                    help="Specify delimiter for GWAS summary statistics. Options = [tab, whitespace, comma]")
parser$add_argument("--snp-col-gwas1-snplist", dest="snp_gwas1_snplist",  required=TRUE,
                    help="Specify SNP column name.")
parser$add_argument("--chr-col-gwas1-snplist", dest="chr_gwas1_snplist",  required=TRUE,
                    help="Specify chromosome column name.")
parser$add_argument("--pos-col-gwas1-snplist", dest="pos_gwas1_snplist",  required=TRUE,
                    help="Specify position column name.")

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
file_gwas1_snplist <- args$file_gwas1_snplist
delim_gwas1_snplist <- args$delim_gwas1_snplist
snp_gwas1_snplist <- args$snp_gwas1_snplist
chr_gwas1_snplist <- args$chr_gwas1_snplist
pos_gwas1_snplist <- args$pos_gwas1_snplist
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
                      wait = TRUE)
  ld_result <- grep("R-sq", plink_out, value = TRUE)
  
  if(nchar(ld_result) != 0){
    ld_result_query <- regmatches(ld_result, regexpr("R-sq\\s+=\\s+([0-9.]+([eE][-+]?[0-9]+)?)", ld_result))
    
    r2 <- as.numeric(strsplit(ld_result_query, " ")[[1]][length(strsplit(ld_result_query, " ")[[1]])])
    
  } else{
    r2 <- NA
  }
  
  return (data.frame(ref_snp = snp1,
                     other_snp = snp2,
                     R2 = r2))
}

#######################
# Read data
#######################

df_gwas2 <- fread(file_gwas2, 
                  sep = delim_mapper[delim_gwas2],
                  data.table = F,
                  nThread = n_thread,
                  showProgress = F)
df_gwas1_snplist <- fread(file_gwas1_snplist, 
                          sep = delim_mapper[delim_gwas1_snplist],
                          data.table = F,
                          nThread = n_thread,
                          showProgress = F)


#######################
# Select necessary columns
#######################
df_gwas2.sub <- df_gwas2 %>%
  dplyr::select(!!as.name(chr_gwas2), !!as.name(pos_gwas2), !!as.name(snp_gwas2), !!as.name(p_gwas2))

#######################
# LD
#######################

file_reference_panel <- reference_panel_list[reference_panel]

df_ld_result <- NULL

for (chr in 1:22){
  df_gwas1_snplist.chr <- filter(df_gwas1_snplist, !!as.name(chr_gwas1_snplist) == chr)
  df_gwas2.chr <- filter(df_gwas2.sub, !!as.name(chr_gwas2) == chr)
  
  if (nrow(df_gwas1_snplist.chr) != 0){
    # df_ld_per_chr comprises the columns below:
    # SNP_gwas1, 
    # SNP_gwas2_tested, SNP_gwas2_tested_pval, SNP_gwas2_tested_R2, 
    # In_LD ## In_LD means R2 of SNP_gwas2 > R2 threshold, hence SNP_gwas2 is in LD with SNP_gwas1, and SNP_gwas1 is not independent from GWAS2 genetic signal. 
    #       ## (So lower R2 threshold to conduct more strict LD comparison)
    # SNP_gwas2_in_LD, SNP_gwas2_in_LD_pval, SNP_gwas2_in_LD_R2
    df_ld_per_chr <- NULL
    for (idx in 1:nrow(df_gwas1_snplist.chr)){
      row.df_gwas1_snplist.chr <- df_gwas1_snplist.chr[idx, ]
      ref_snp <- row.df_gwas1_snplist.chr[[snp_gwas1_snplist]]
      ref_snp.pos <- row.df_gwas1_snplist.chr[[pos_gwas1_snplist]]
      
      df_gwas2.chr.filter <- filter(df_gwas2.chr, (!!as.name(pos_gwas2) >= ref_snp.pos - (window * 1000)) &
                                                  (!!as.name(pos_gwas2) < ref_snp.pos + (window * 1000)) &
                                                  (!!as.name(p_gwas2) < p_threshold)
                                    
                                    )
      
      if(nrow(df_gwas2.chr.filter) != 0){
        other_snplist <- df_gwas2.chr.filter[[snp_gwas2]]
        
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
        df_tmp <- df_gwas2.chr.filter %>%
                    dplyr::select(!!as.name(snp_gwas2), !!as.name(p_gwas2)) %>%
                    filter(!!as.name(snp_gwas2) %in% snp_gwas2_tested)
        df_tmp[[snp_gwas2]] <- factor(df_tmp[[snp_gwas2]], levels = snp_gwas2_tested)
        snp_gwas2_tested_pval <- df_tmp[[p_gwas2]]
        snp_gwas2_tested_r2 <- plink_ld_result$R2
      
        idx_in_ld <- which(snp_gwas2_tested_r2 > r2_threshold)
        if (identical(which(idx_in_ld > r2_threshold), integer(0))){ # No SNPs from GWAS2 in LD with SNP from GWAS1
          is_in_LD <- FALSE
        } else{
          is_in_LD <- TRUE
          snp_gwas2_in_ld <- paste(snp_gwas2_tested[idx_in_ld], collapse = ", ")
          snp_gwas2_in_ld_pval <- paste(snp_gwas2_tested_pval[idx_in_ld], collapse = ", ")
          snp_gwas2_in_ld_r2 <- paste(snp_gwas2_tested_r2[idx_in_ld], collapse = ", ")
        }
        
        snp_gwas2_tested <- paste(snp_gwas2_tested, collapse = ", ")
        snp_gwas2_tested_pval <- paste(snp_gwas2_tested_pval, collapse = ", ")
        snp_gwas2_tested_r2 <- paste(snp_gwas2_tested_r2, collapse = ", ")
        

        df_ld_per_chr <- rbind(df_ld_per_chr,
                        data.frame(SNP_gwas1 = ref_snp,
                            SNP_gwas2_tested = snp_gwas2_tested,
                            SNP_gwas2_tested_pval = snp_gwas2_tested_pval,
                            SNP_gwas2_tested_R2 = snp_gwas2_tested_r2,
                            In_LD = is_in_LD,
                            SNP_gwas2_in_LD = snp_gwas2_in_ld,
                            SNP_gwas2_in_LD_pval = snp_gwas2_in_ld_pval,
                            SNP_gwas2_in_LD_R2 = snp_gwas2_in_ld_r2)
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
                            SNP_gwas2_in_LD_R2 = "NA")
                      )
      }
    }

    df_ld_result <- rbind(df_ld_result, df_ld_per_chr)
  }
}

system("rm plink.log plink.nosex",
      wait = TRUE)

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

write.table(df_ld_result,
            file = paste0(dir_out, pref_out, ".", exe_mapper[delim_out]),
            sep = delim_mapper[delim_out],
            row.names = F,
            quote = quote)
