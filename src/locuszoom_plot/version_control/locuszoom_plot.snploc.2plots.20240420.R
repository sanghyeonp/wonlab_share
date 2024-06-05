# :: LocusZoom plot with snp.loc ::
# - Sanghyeon Park
# - 2024.03.31
# - First in-use project: GWAS-by-subtraction of SCZ 
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
parser$add_argument("--trait1", required=F, default="trait1", 
                    help="Specify the trait name for GWAS 1. Default=trait1.")
parser$add_argument("--gwas1", required=T, 
                    help="Specify the GWAS 1 file path.")
parser$add_argument("--delim-gwas1", dest="delim.gwas1", required=T,
                    help="Specify delimiter for GWAS summary statistics. Options = [tab, whitespace, comma]")
parser$add_argument("--snp-col-gwas1", dest="snp.gwas1", required=T,
                    help="Specify SNP column name of GWAS1.")
parser$add_argument("--chr-col-gwas1", dest="chr.gwas1", required=T,
                    help="Specify CHR column name of GWAS1.")
parser$add_argument("--pos-col-gwas1", dest="pos.gwas1", required=T,
                    help="Specify POS column name of GWAS1.")
parser$add_argument("--pval-col-gwas1", dest="p.gwas1", required=T,
                    help="Specify P-value column name of GWAS1.")

# GWAS2
parser$add_argument("--trait2", required=F, default="trait2", 
                    help="Specify the trait name for GWAS 2. Default=trait2.")
parser$add_argument("--gwas2", required=T, 
                    help="Specify the GWAS 2 file path.")
parser$add_argument("--delim-gwas2", dest="delim.gwas2", required=T,
                    help="Specify delimiter for GWAS summary statistics. Options = [tab, whitespace, comma]")
parser$add_argument("--snp-col-gwas2", dest="snp.gwas2", required=T,
                    help="Specify SNP column name of GWAS2.")
parser$add_argument("--chr-col-gwas2", dest="chr.gwas2", required=T,
                    help="Specify CHR column name of GWAS2.")
parser$add_argument("--pos-col-gwas2", dest="pos.gwas2", required=T,
                    help="Specify POS column name of GWAS2.")
parser$add_argument("--pval-col-gwas2", dest="p.gwas2", required=T,
                    help="Specify P-value column name of GWAS2.")

# SNPLOC file
parser$add_argument("--snploc", required=T, 
                    help="Specify the SNP location file. It should contain SNP, CHR, START POS, END POS informations.")
parser$add_argument("--delim-snploc", dest="delim.snploc", required=T,
                    help="Specify delimiter for snploc file. Options = [tab, whitespace, comma]")
parser$add_argument("--snp-col-snploc", dest="snp.snploc", required=T,
                    help="Specify SNP column name of snploc file.")
parser$add_argument("--chr-col-snploc", dest="chr.snploc", required=T,
                    help="Specify CHR column name of snploc file.")
parser$add_argument("--start-col-snploc", dest="start.snploc", required=T,
                    help="Specify START POS column name of snploc file.")
parser$add_argument("--end-col-snploc", dest="end.snploc", required=T,
                    help="Specify END POS column name of snploc file.")

# Parallel processing
parser$add_argument("--n-thread", dest="n_thread",  required=F, default = 1, type = "integer",
                    help="Specify the number of threads.")

# Plot-related
parser$add_argument("--gene-track-off", dest="gene_track_off", action="store_true", default=F,
                    help="Specify to avoid plotting the gene track.")
parser$add_argument("--ylab1", required=F, default="NA",
                    help="Specify the y-axis title for GWAS1. Default=-log10(P).")
parser$add_argument("--ylab2", required=F, default="NA",
                    help="Specify the y-axis title for GWAS2. Default=-log10(P).")

# Output
parser$add_argument("--out-dir", dest="dir_out", required=F, default=".",
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

file.snploc <- args$snploc
delim.snploc <- args$delim.snploc
snp.snploc <- args$snp.snploc
chr.snploc <- args$chr.snploc
start.snploc <- args$start.snploc
end.snploc <- args$end.snploc

n_thread <- args$n_thread

gene_track_off <- args$gene_track_off
ylab1 <- args$ylab1; ylab1 <- str_na_to_null(ylab1) # Default
ylab2 <- args$ylab2; ylab2 <- str_na_to_null(ylab2)


dir_out <- args$dir_out; dir_out <- add_trailing_slash(dir_out); dir.create(dir_out, showWarnings=F)

#######

df.gwas1 <- fread(file.gwas1, sep=delim_map[delim.gwas1], data.table=F, nThread=n_thread) %>%
    dplyr::select(!!as.name(snp.gwas1), !!as.name(chr.gwas1), !!as.name(pos.gwas1), !!as.name(p.gwas1)) %>%
    rename(rsid=!!as.name(snp.gwas1),
           chrom=!!as.name(chr.gwas1),
           pos=!!as.name(pos.gwas1),
           p=!!as.name(p.gwas1)) %>%
    mutate(chrom = as.character(chrom))

df.gwas2 <- fread(file.gwas2, sep=delim_map[delim.gwas2], data.table=F, nThread=n_thread) %>%
    dplyr::select(!!as.name(snp.gwas2), !!as.name(chr.gwas2), !!as.name(pos.gwas2), !!as.name(p.gwas2)) %>%
    rename(rsid=!!as.name(snp.gwas2),
           chrom=!!as.name(chr.gwas2),
           pos=!!as.name(pos.gwas2),
           p=!!as.name(p.gwas2)) %>%
    mutate(chrom = as.character(chrom))

df.snploc <- fread(file.snploc, sep=delim_map[delim.snploc], data.table=F, nThread=n_thread) %>%
    rename(SNP=!!as.name(snp.snploc),
           CHR=!!as.name(chr.snploc),
           START=!!as.name(start.snploc),
           END=!!as.name(end.snploc))


######

for (idx in 1:nrow(df.snploc)){
    snp1 <- df.snploc[idx, ]$SNP
    chr <- as.chracter(df.snploc[idx, ]$CHR)
    start_pos <- df.snploc[idx, ]$START
    end_pos <- df.snploc[idx, ]$END
    
    ## Subset GWAS1 and compute R2 between selected SNP
    df.gwas1.sub <- df.gwas1 %>%
        dplyr::filter(chrom == chr & pos >= start_pos & pos <= end_pos) %>%
        dplyr::select(chrom, pos, rsid, p)
    
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
    df.gwas2.sub <- df.gwas2 %>%
        dplyr::filter(chrom == chr & pos >= start_pos & pos <= end_pos) %>%
        dplyr::select(chrom, pos, rsid, p)
    
    if (snp1 %in% df.gwas2.sub$rsid){
        snp2 <- snp1
    } else{
        snp2 <- (df.gwas2.sub %>% arrange(p))[1, ]$rsid
    }
    
    # Register a parallel backend
    cl <- makeCluster(n_thread)
    registerDoParallel(cl)

    df.gwas2.ld <- foreach(other_snp=df.gwas1.sub$rsid,
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
    
    ## Make loc object for GWAS1 and GWAS2
    if (require(EnsDb.Hsapiens.v75)){
        loc1 <- locus(data=df.gwas1.sub,
                      seqname=chr,
                      xrange=c(start_pos, end_pos),
                      ens_db="EnsDb.Hsapiens.v75",
                      LD="R2",
                      index_snp=snp1)
    }
    if (require(EnsDb.Hsapiens.v75)){
        loc2 <- locus(data=df.gwas2.sub,
                      seqname=chr,
                      xrange=c(start_pos, end_pos),
                      ens_db="EnsDb.Hsapiens.v75",
                      LD="R2",
                      index_snp=snp2)
    }
    
    ## 
    pdf(paste0(dir_out, "locuszoom.trait1_", trait.gwas1, ".trait2_", trait.gwas2,
               ".chr", chr, ".snp1_", snp1, ".snp2_", snp2, ".pdf"), 
        width=6, height=8)
    order_partition <- set_layers(2)
    
    if (snp1 == snp2){
        scatter_plot(loc1, labels=c("index"), 
                     border=FALSE, xtick=FALSE, legend_pos="topright",
                     ylab=ylab1)
        scatter_plot(loc2, labels=c("index"),
                     border=FALSE, xtick=FALSE, legend_pos="topright",
                     ylab=ylab2)
    } else{
        # For Plot 1
        if (snp2 %in% loc1$data$rsid){
            scatter_plot(loc1, labels=c("index", snp2),
                         border=FALSE, xtick=FALSE, legend_pos="topright",
                         ylab=ylab1)
        } else{
            scatter_plot(loc1, labels=c("index"),
                         border=FALSE, xtick=FALSE, legend_pos="topright",
                         ylab=ylab1)
        }
        # For Plot 2
        if (snp1 %in% loc2$data$rsid){
            scatter_plot(loc2, labels=c("index", snp1),
                         border=FALSE, xtick=FALSE, legend_pos="topright",
                         ylab=ylab2)
        } else{
            scatter_plot(loc2, labels=c("index"),
                         border=FALSE, xtick=FALSE, legend_pos="topright",
                         ylab=ylab2)
        }

    }
    
    if (!gene_track_off){
        genetracks(loc1, border=FALSE,
                   filter_gene_biotype="protein_coding",
                   gene_col="grey", exon_col="blue", exon_border="darkgrey")
    }

    par(order_partition)
    dev.off()
    
}
