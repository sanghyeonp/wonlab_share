# :: LocusZoom plot with gene.loc ::
# - Sanghyeon Park
# - 2024.04.20
# - First in-use project: SSc drug repurposing
# - Trait 1이 위에.

# library(BiocManager)
# BiocManager::install("ensembldb")
# BiocManager::install("EnsDb.Hsapeins.v75")
# install.packages("EnsDb.Hsapiens.v75_2.99.0.tar.gz", repos=NULL, type="source")
# install.packages("locuszoomr")
library(data.table)
library(dplyr)
httr::set_config(httr::config(ssl_verifypeer=FALSE))
library(locuszoomr)
library(foreach)
library(doParallel)
# library(cowplot)
options(error=traceback)
require(argparse)


### Pre-defined
exe_plink <- "/data1/sanghyeon/wonlab_contribute/combined/software/plink/plink"
reference <- "/data1/sanghyeon/wonlab_contribute/combined/data_common/reference_panel/1kGp3/EUR/CHR/reference.1kG.EUR.maf_0.005.geno_0.02.chr"

delim_map <- c("tab"="\t", "comma"=",", "whitespace"=" ")

add_trailing_slash <- function(str1){
    if(!grepl("/$", str1)){
        str1 <- paste0(str1, "/")
    }
    return (str1)
}

str_na_to_null <- function(v1){
    if (v1 == "NA"){
        v1 <- NULL
    }
    return (v1)
}

compute_snp_ld <- function(snp1, snp2, plink_exe, file_reference){
    if (snp1 == snp2){
        r2 <- 1
    } else{
        plink_out <- system(paste0(plink_exe, 
                                   " --bfile ", file_reference,
                                   " --ld ", snp1, " ", snp2),
                            intern=TRUE,
                            wait=TRUE,
                            ignore.stderr=TRUE)
        ld_result <- grep("R-sq", plink_out, value=TRUE)
        
        if (identical(ld_result, character(0))){
            r2 <- NA
        } else{
            ld_result <- ld_result[1] # Refer to important note 1
            if (nchar(ld_result) != 0){
                ld_result_query <- regmatches(ld_result, regexpr("R-sq\\s+=\\s+([0-9.]+([eE][-+]?[0-9]+)?)", ld_result))
                
                r2 <- as.numeric(strsplit(ld_result_query, " ")[[1]][length(strsplit(ld_result_query, " ")[[1]])])
                
            } else{
                r2 <- NA
            }
        }
    }
    
    return (data.frame(ref_snp=snp1,
                       other_snp=snp2,
                       R2=r2))
}

####
#######################
# Parser arguments
#######################
### Define parser arguments
parser <- argparse::ArgumentParser(description=":: LocusZoom two plots given snp.loc ::", 
                                   formatter_class="argparse.ArgumentDefaultsHelpFormatter")

# GWAS1
parser$add_argument("--trait1", required=F, default="trait1", type="character",
                    help="Specify the trait name for GWAS 1. Default=trait1.")
parser$add_argument("--gwas1", required=T, type="character",
                    help="Specify the GWAS 1 file path.")
parser$add_argument("--delim-gwas1", dest="delim.gwas1", required=T, type="character",
                    help="Specify delimiter for GWAS summary statistics. Options = [tab, whitespace, comma]")
parser$add_argument("--snp-col-gwas1", dest="snp.gwas1", required=T, type="character",
                    help="Specify SNP column name of GWAS1.")
parser$add_argument("--chr-col-gwas1", dest="chr.gwas1", required=T, type="character",
                    help="Specify CHR column name of GWAS1.")
parser$add_argument("--pos-col-gwas1", dest="pos.gwas1", required=T, type="character",
                    help="Specify POS column name of GWAS1.")
parser$add_argument("--pval-col-gwas1", dest="p.gwas1", required=T, type="character",
                    help="Specify P-value column name of GWAS1.")

# GWAS2
parser$add_argument("--trait2", required=F, default="trait2", type="character",
                    help="Specify the trait name for GWAS 2. Default=trait2.")
parser$add_argument("--gwas2", required=T,  type="character",
                    help="Specify the GWAS 2 file path.")
parser$add_argument("--delim-gwas2", dest="delim.gwas2", required=T, type="character",
                    help="Specify delimiter for GWAS summary statistics. Options = [tab, whitespace, comma]")
parser$add_argument("--snp-col-gwas2", dest="snp.gwas2", required=T, type="character",
                    help="Specify SNP column name of GWAS2.")
parser$add_argument("--chr-col-gwas2", dest="chr.gwas2", required=T, type="character",
                    help="Specify CHR column name of GWAS2.")
parser$add_argument("--pos-col-gwas2", dest="pos.gwas2", required=T, type="character",
                    help="Specify POS column name of GWAS2.")
parser$add_argument("--pval-col-gwas2", dest="p.gwas2", required=T, type="character",
                    help="Specify P-value column name of GWAS2.")

# GENELOC file
parser$add_argument("--geneloc", required=T, type="character",
                    help="Specify the SNP location file. It should contain SNP, CHR, START POS, END POS informations.")
parser$add_argument("--delim-geneloc", dest="delim.geneloc", required=T, type="character",
                    help="Specify delimiter for geneloc file. Options = [tab, whitespace, comma]")
parser$add_argument("--gene-col-geneloc", dest="gene.geneloc", required=T, type="character",
                    help="Specify GENE column name of geneloc file.")
parser$add_argument("--chr-col-geneloc", dest="chr.geneloc", required=T, type="character",
                    help="Specify CHR column name of geneloc file.")
parser$add_argument("--start-col-geneloc", dest="start.geneloc", required=T, type="character",
                    help="Specify START POS column name of geneloc file.")
parser$add_argument("--end-col-geneloc", dest="end.geneloc", required=T, type="character",
                    help="Specify END POS column name of geneloc file.")

# Extract gene
parser$add_argument("--extract-gene-gwas1", dest="extract_gene_gwas1", action="store_true", default=F,
                    help="Specify to extract gene from the given trait1 summary statistics. Default=FALSE.")
parser$add_argument("--gene-extract-col-gwas1", dest="gene.extract.gwas1", default="NA", required=F, type="character",
                    help="Specify the gene column in GWAS1 to extract the specified gene.")
parser$add_argument("--extract-gene-gwas2", dest="extract_gene_gwas2", action="store_true", default=F,
                    help="Specify to extract gene from the given trait2 summary statistics. Default=FALSE.")
parser$add_argument("--gene-extract-col-gwas2", dest="gene.extract.gwas2", default="NA", required=F, type="character",
                    help="Specify the gene column in GWAS2 to extract the specified gene.")
parser$add_argument("--gene-extract-col-geneloc", dest="gene.extract.geneloc", default="NA", required=F, type="character",
                    help="Specify the gene column in GENELOC to extract the specified gene.")
# Parallel processing
parser$add_argument("--n-thread", dest="n_thread",  required=F, default = 1, type = "integer",
                    help="Specify the number of threads.")

# Plot-related
parser$add_argument("--flanking-window", dest="flanking_window", type="integer", required=T,
                    help="Specify the flanking window in kb.")

# Output
parser$add_argument("--out-dir", dest="dir_out", required=F, default=".", type="character",
                    help="Specify output directory. Default=current working directory.")

####
args <- parser$parse_args()

trait.gwas1 <- args$trait1
file.gwas1 <- args$gwas1
delim.gwas1 <- args$delim.gwas1
snp.gwas1 <- args$snp.gwas1
chr.gwas1 <- args$chr.gwas1
pos.gwas1 <- args$pos.gwas1
p.gwas1 <- args$p.gwas1

trait.gwas2 <- args$trait2
file.gwas2 <- args$gwas2
delim.gwas2 <- args$delim.gwas2
snp.gwas2 <- args$snp.gwas2
chr.gwas2 <- args$chr.gwas2
pos.gwas2 <- args$pos.gwas2
p.gwas2 <- args$p.gwas2

file.geneloc <- args$geneloc
delim.geneloc <- args$delim.geneloc
gene.geneloc <- args$gene.geneloc
chr.geneloc <- args$chr.geneloc
start.geneloc <- args$start.geneloc
end.geneloc <- args$end.geneloc

n_thread <- args$n_thread

extract_gene_gwas1 <- args$extract_gene_gwas1
gene.extract.gwas1 <- args$gene.extract.gwas1
if (extract_gene_gwas1 & gene.extract.gwas1 == "NA"){
    stop("--extract-gene-gwas1 is specified. Provide gene column to extract in GWAS1 by --gene-gwas1.")
}
extract_gene_gwas2 <- args$extract_gene_gwas2
gene.extract.gwas2 <- args$gene.extract.gwas2
if (extract_gene_gwas2 & gene.extract.gwas2 == "NA"){
    stop("--extract-gene-gwas2 is specified. Provide gene column to extract in GWAS2 by --gene-gwas1.")
}
gene.extract.geneloc <- args$gene.extract.geneloc
if ((extract_gene_gwas1 | extract_gene_gwas2) & gene.extract.geneloc == "NA"){
    stop("--extract-gene-gwas1 or --extract-gene-gwas2 is specified. Provide gene column to extract in GENELOC by --gene-extract-geneloc.")
}

flanking_window <- args$flanking_window
dir_out <- args$dir_out; dir_out <- add_trailing_slash(dir_out); dir.create(dir_out, showWarnings=F)

#######
if (extract_gene_gwas1){
    df.gwas1 <- fread(file.gwas1, sep=delim_map[delim.gwas1], data.table=F, nThread=n_thread) %>%
        dplyr::select(!!as.name(gene.extract.gwas1), !!as.name(snp.gwas1), !!as.name(chr.gwas1), !!as.name(pos.gwas1), !!as.name(p.gwas1)) %>%
        dplyr::rename(GENE.extract = !!as.name(gene.extract.gwas1),
                      rsid=!!as.name(snp.gwas1),
                      chrom=!!as.name(chr.gwas1),
                      pos=!!as.name(pos.gwas1),
                      p=!!as.name(p.gwas1)) %>%
        dplyr::mutate(chrom = as.character(chrom))
} else{
    df.gwas1 <- fread(file.gwas1, sep=delim_map[delim.gwas1], data.table=F, nThread=n_thread) %>%
        dplyr::select(!!as.name(snp.gwas1), !!as.name(chr.gwas1), !!as.name(pos.gwas1), !!as.name(p.gwas1)) %>%
        dplyr::rename(rsid=!!as.name(snp.gwas1),
                      chrom=!!as.name(chr.gwas1),
                      pos=!!as.name(pos.gwas1),
                      p=!!as.name(p.gwas1)) %>%
        dplyr::mutate(chrom = as.character(chrom))
}

if (extract_gene_gwas2){
    df.gwas2 <- fread(file.gwas2, sep=delim_map[delim.gwas2], data.table=F, nThread=n_thread) %>%
        dplyr::select(!!as.name(gene.extract.gwas2), !!as.name(snp.gwas2), !!as.name(chr.gwas2), !!as.name(pos.gwas2), !!as.name(p.gwas2)) %>%
        rename(GENE.extract = !!as.name(gene.extract.gwas2),
               rsid=!!as.name(snp.gwas2),
               chrom=!!as.name(chr.gwas2),
               pos=!!as.name(pos.gwas2),
               p=!!as.name(p.gwas2)) %>%
        mutate(chrom = as.character(chrom))
} else{
    df.gwas2 <- fread(file.gwas2, sep=delim_map[delim.gwas2], data.table=F, nThread=n_thread) %>%
        dplyr::select(!!as.name(snp.gwas2), !!as.name(chr.gwas2), !!as.name(pos.gwas2), !!as.name(p.gwas2)) %>%
        rename(rsid=!!as.name(snp.gwas2),
               chrom=!!as.name(chr.gwas2),
               pos=!!as.name(pos.gwas2),
               p=!!as.name(p.gwas2)) %>%
        mutate(chrom = as.character(chrom))
}

if (extract_gene_gwas1 | extract_gene_gwas2){
    df.geneloc <- fread(file.geneloc, sep=delim_map[delim.geneloc], data.table=F, nThread=n_thread) %>%
        dplyr::select(!!as.name(gene.extract.geneloc), !!as.name(gene.geneloc), !!as.name(chr.geneloc), !!as.name(start.geneloc), !!as.name(end.geneloc)) %>%
        dplyr::rename(GENE.extract = !!as.name(gene.extract.geneloc),
                GENE=!!as.name(gene.geneloc),
               CHR=!!as.name(chr.geneloc),
               START=!!as.name(start.geneloc),
               END=!!as.name(end.geneloc)) %>%
        dplyr::mutate(CHR = as.character(CHR))
} else{
    df.geneloc <- fread(file.geneloc, sep=delim_map[delim.geneloc], data.table=F, nThread=n_thread) %>%
        dplyr::select(!!as.name(gene.geneloc), !!as.name(chr.geneloc), !!as.name(start.geneloc), !!as.name(end.geneloc)) %>%
        dplyr::rename(GENE=!!as.name(gene.geneloc),
                      CHR=!!as.name(chr.geneloc),
                      START=!!as.name(start.geneloc),
                      END=!!as.name(end.geneloc)) %>%
        dplyr::mutate(CHR = as.character(CHR))
}

######

for (idx in 1:nrow(df.geneloc)){
    gene <- df.geneloc[idx, ]$GENE
    chr <- as.character(df.geneloc[idx, ]$CHR)
    start_pos <- df.geneloc[idx, ]$START
    end_pos <- df.geneloc[idx, ]$END
    
    ## Subset GWAS1 and compute R2 between selected SNP
    if (extract_gene_gwas1){
        gene_to_extract <- df.geneloc[idx, ]$GENE.extract
        df.gwas1.sub <- df.gwas1 %>%
            dplyr::filter(GENE.extract == gene_to_extract) %>%
            dplyr::filter(chrom == chr & pos >= (start_pos - (flanking_window*1000)) & pos <= (end_pos + flanking_window*1000)) %>%
            dplyr::select(chrom, pos, rsid, p)
    } else{
        df.gwas1.sub <- df.gwas1 %>%
            dplyr::filter(chrom == chr & pos >= (start_pos - (flanking_window*1000)) & pos <= (end_pos + flanking_window*1000)) %>%
            dplyr::select(chrom, pos, rsid, p)
    }
    snp1 <- df.gwas1.sub %>%
        arrange(p) %>%
        dplyr::select(rsid) %>%
        dplyr::slice(1) %>%
        dplyr::pull(1)
    
    # Register a parallel backend
    cl <- makeCluster(n_thread)
    registerDoParallel(cl)

    df.gwas1.ld <- foreach(other_snp=df.gwas1.sub$rsid,
                            .export=c("snp1", "exe_plink", "reference", "chr"),
                            .packages=c("dplyr", "data.table"),
                            .combine='bind_rows') %dopar% {
                        compute_snp_ld(snp1=snp1, 
                                        snp2=other_snp, 
                                        plink_exe=exe_plink, 
                                        file_reference=paste0(reference, chr))
                        }
    stopCluster(cl)

    df.gwas1.sub <- merge(df.gwas1.sub, df.gwas1.ld %>% dplyr::select(-ref_snp),
                          by.x="rsid", by.y="other_snp", all.x=T)
    
    ## Subset GWAS2 and compute R2 between selected SNP
    if (extract_gene_gwas2){
        gene_to_extract <- df.geneloc[idx, ]$GENE.extract
        df.gwas2.sub <- df.gwas2 %>%
            dplyr::filter(GENE.extract == gene_to_extract) %>%
            dplyr::filter(chrom == chr & pos >= (start_pos - flanking_window*1000) & pos <= (end_pos + flanking_window*1000)) %>%
            dplyr::select(chrom, pos, rsid, p)
    } else{
        df.gwas2.sub <- df.gwas2 %>%
            dplyr::filter(chrom == chr & pos >= (start_pos - flanking_window*1000) & pos <= (end_pos + flanking_window*1000)) %>%
            dplyr::select(chrom, pos, rsid, p)
    }
    df.gwas2.sub <- df.gwas2 %>%
        dplyr::filter(chrom == chr & pos >= (start_pos - flanking_window*1000) & pos <= (end_pos + flanking_window*1000)) %>%
        dplyr::select(chrom, pos, rsid, p)
    
    
    if (snp1 %in% df.gwas2.sub$rsid){
        snp2 <- snp1
    } else{
        snp2 <- (df.gwas2.sub %>% arrange(p))[1, ]$rsid
    }
    
    # Register a parallel backend
    cl <- makeCluster(n_thread)
    registerDoParallel(cl)

    df.gwas2.ld <- foreach(other_snp=df.gwas2.sub$rsid,
                            .export=c("snp1", "exe_plink", "reference", "chr"),
                            .packages=c("dplyr", "data.table"),
                            .combine='bind_rows') %dopar% {
                        compute_snp_ld(snp1=snp1, 
                                        snp2=other_snp, 
                                        plink_exe=exe_plink, 
                                        file_reference=paste0(reference, chr))
                        }
    stopCluster(cl)
    
    df.gwas2.sub <- merge(df.gwas2.sub, df.gwas2.ld %>% dplyr::select(-ref_snp),
                          by.x="rsid", by.y="other_snp", all.x=T)

    save(df.gwas1.sub, df.gwas2.sub, snp1, snp2, 
        file=paste0(dir_out, "locuszoom.trait1_", trait.gwas1, ".trait2_", trait.gwas2,
                    ".chr", chr, ".snp1_", snp1, ".snp2_", snp2, ".loc_input.RData"))
}