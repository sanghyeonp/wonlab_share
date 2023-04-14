# Load packages from env_R.R
if(!require(rstudioapi)){
    install.packages("rstudioapi", repos = "http://cran.us.r-project.org")
    library(rstudioapi)
}
current <- rstudioapi::getSourceEditorContext()$path
env_file <- paste(dirname(current), "env_R.R", sep="/")
source(env_file)


### Define parser arguments
parser <- argparse::ArgumentParser(description=":: Compute Cohen's kappa coefficient ::", 
                                  formatter_class="argparse.ArgumentDefaultsHelpFormatter")

parser$add_argument("--ref_gwas", required=TRUE,
                    help="Path to the reference GWAS summary statistics.")
parser$add_argument("--ref_gwas_delim", required=FALSE, default="tab",
                    help="Specify the delimiter used for the reference GWAS. Options are [tab, comma, whitespace]. Default='tab'.")
parser$add_argument("--ref_gwas_snp", required=FALSE, default="SNP",
                    help="Specify the SNP column name in the reference GWAS summary statistics. Default='SNP'.")
parser$add_argument("--ref_gwas_beta", required=FALSE, default="Beta",
                    help="Specify the Beta column name in the reference GWAS summary statistics. Default='Beta'.")
parser$add_argument("--ref_name", required=FALSE, default="NA",
                    help="Name of the reference trait.")

parser$add_argument("--alt_gwas", required=TRUE, 
                    help="Path to the reference GWAS summary statistics.")
parser$add_argument("--alt_gwas_delim", required=FALSE, default="tab",
                    help="Specify the delimiter used for the alternative GWAS. Options are [tab, comma, whitespace]. Default='tab'.")
parser$add_argument("--alt_gwas_snp", required=FALSE, default="SNP",
                    help="Specify the SNP column name in the alternative summary statistics. Default='SNP'.")
parser$add_argument("--alt_gwas_beta", required=FALSE, default="Beta",
                    help="Specify the Beta column name in the alternative GWAS summary statistics. Default='Beta'.")
parser$add_argument("--alt_name", required=FALSE, default="NA",
                    help="Name of the alternative trait.")
parser$add_argument("--alt_rev_beta", action='store_true',
                    help="Specify to reverse the effect direction of alternative GWAS. Default = FALSE.")

parser$add_argument("--snplist", required=FALSE, default="NA",
                    help="Path to file with SNP list.")
parser$add_argument("--snplist_FUMA", required=FALSE, default="NA", 
                    help="Path to `leadSNPs.txt` from FUMA to compute Cohen's kappa coefficient on the lead SNPs identified from FUMA.")

parser$add_argument("--outf", required=FALSE, default="cohenkappa",
                    help="Name of the output log.")
parser$add_argument("--outd", required=FALSE, default="NA",
                    help="Path to output directory. Default = current working directory.")

parser$add_argument("--verbose", action='store_true',
                    help="")
parser$add_argument("--rds", action='store_true',
                    help="Save RDS for the Cohen's kappa test result. Default = FALSE.")
parser$add_argument("--table", action='store_true',
                    help="Save table in csv for the Cohen's kappa test result. Default = FALSE.")

### Get parser arguments
args <- parser$parse_args()

ref_gwas <- args$ref_gwas
ref_gwas_delim <- args$ref_gwas_delim
ref_gwas_snp <- args$ref_gwas_snp
ref_gwas_beta <- args$ref_gwas_beta
ref_name <- args$ref_name

alt_gwas <- args$alt_gwas
alt_gwas_delim <- args$alt_gwas_delim
alt_gwas_snp <- args$alt_gwas_snp
alt_gwas_beta <- args$alt_gwas_beta
alt_name <- args$alt_name
alt_rev_beta <- args$alt_rev_beta

snplist <- args$snplist
snplist_fuma <- args$snplist_FUMA

outf <- args$outf
outd <- args$outd
if (outd == "NA"){
  outd <- paste0(getwd(), "/")
}

verbose <- args$verbose
rds <- args$rds
totable <- args$table

### Logs
save_log <- function(log_list, log_, verbose = FALSE){
  log_ <- paste0(log_, "\n")
  if (verbose){
    cat(log_)
  }
  return(c(log_list, log_))
}

### Re-map the delimiters
delim_map <- hash()
delim_map[["tab"]] <- "\t"
delim_map[["comma"]] <- ","
delim_map[["whitespace"]] <- " "

ref_gwas_delim <- delim_map[[ref_gwas_delim]]
alt_gwas_delim <- delim_map[[alt_gwas_delim]]

### Process SNP list read
logs2 <- c()
if (snplist == "NA" & snplist_fuma == "NA"){
  log_ <- paste0("SNP list is NA. All SNPs in the reference GWAS are compared."); logs2 <- save_log(logs2, log_)
  snplist_ <- df_ref$snp
} else{
  if (snplist_fuma == "NA"){
    snplist_ <- fread(snplist, header = F, stringsAsFactors = F)$V1
    log_ <- paste0("SNP list was provided. ", length(snplist_), " SNPs were specified."); logs2 <- save_log(logs2, log_)
  } else{
    snplist_ <- fread(snplist_fuma, sep="\t", header = T, stringsAsFactors = F)$rsID
    log_ <- paste0("SNP list was provided from FUMA `leadSNPs.txt`. ", length(snplist_), " SNPs were specified."); logs2 <- save_log(logs2, log_)
  }
}

# SNP list into a single string
snplist_ <- paste(c(snplist_), collapse=",")

temp <- main(ref_gwas, ref_gwas_delim, ref_gwas_snp, ref_gwas_beta,
                alt_gwas, alt_gwas_delim, alt_gwas_snp, alt_gwas_beta, alt_rev_beta,
                snplist_, ref_name, alt_name)
cohen <- temp[[1]]

### Print the log
logs <- temp[[2]]
logs_final <- c(logs[1:1], logs2, logs[2:length(logs)])
if (verbose){
  for (log_ in logs_final){
    cat(log_)
  }
  cat("\n")
  print(cohen)
}

### Save the log
sink(paste0(outd, outf, ".log"))
for (log_ in logs_final){
  cat(log_)
}
cat("\n")
print(cohen)
sink()

### Save RDS
if (rds){
  saveRDS(cohen, paste0(outd, outf, ".RDS"))
}

### Save table
if (totable){
  result <- temp[[3]]
  write.csv(result, paste0(outd, outf, ".csv"), row.names=FALSE, quote=TRUE)
}
