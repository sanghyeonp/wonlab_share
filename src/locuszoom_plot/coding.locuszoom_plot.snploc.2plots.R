# library(BiocManager)
# BiocManager::install("ensembldb")
# BiocManager::install("EnsDb.Hsapeins.v75")
# install.packages("EnsDb.Hsapiens.v75_2.99.0.tar.gz", repos=NULL, type="source")
# install.packages("locuszoomr")
library(data.table)
library(dplyr)
httr::set_config(httr::config(ssl_verifypeer = FALSE))
library(locuszoomr)
library(cowplot)
options(error = traceback)


### Pre-defined
exe_plink <- "./plink_mac_20231211/plink"
reference <- "./1kg/reference.1kG.EUR.maf_0.005.geno_0.02.chr19"

delim_map <- c("tab" = "\t", "comma" = ",", "space" = " ")

add_trailing_slash <- function(str1){
    if(!grepl("/$", str1)){
        str1 <- paste0(str1, "/")
    }
    return(str1)
}

compute_snp_ld <- function(snp1, snp2, plink_exe, file_reference){
    if (snp1 == snp2){
        r2 <- 1
    } else{
        plink_out <- system(paste0(plink_exe, 
                                   " --bfile ", file_reference,
                                   " --ld ", snp1, " ", snp2),
                            intern = TRUE,
                            wait = TRUE,
                            ignore.stderr = TRUE)
        ld_result <- grep("R-sq", plink_out, value = TRUE)
        
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
    
    return (data.frame(ref_snp = snp1,
                       other_snp = snp2,
                       R2 = r2))
}

####
trait.gwas1 <- "SCZsubBD"
file.gwas1 <- "SCZsubBD.chr19.tsv"
delim.gwas1 <- "tab"
snp.gwas1 <- "SNP"
chr.gwas1 <- "CHR"
pos.gwas1 <- "BP"
p.gwas1 <- "Pval_Estimate"

trait.gwas2 <- "SCZ"
file.gwas2 <- "SCZ.chr19.tsv"
delim.gwas2 <- "tab"
snp.gwas2 <- "SNP"
chr.gwas2 <- "CHR"
pos.gwas2 <- "POS"
p.gwas2 <- "PVAL"

file.snploc <- "snp.loc"
delim.snploc <- "tab"
snp.snploc <- "snp"
chr.snploc <- "chr"
start.snploc <- "start"
end.snploc <- "end"

n_thread <- 1

dir_out <- "./locuszoom.snp_indpendence/"; dir_out <- add_trailing_slash(dir_out); dir.create(dir_out, showWarnings=F)

#######

df.gwas1 <- fread(file.gwas1, sep=delim_map[delim.gwas1], data.table=F, nThread=n_thread) %>%
    rename(rsid = !!as.name(snp.gwas1),
           chrom = !!as.name(chr.gwas1),
           pos = !!as.name(pos.gwas1),
           p = !!as.name(p.gwas1)) %>%
    dplyr::select(rsid, chrom, pos, p)

df.gwas2 <- fread(file.gwas2, sep=delim_map[delim.gwas2], data.table=F, nThread=n_thread) %>%
    rename(rsid = !!as.name(snp.gwas2),
           chrom = !!as.name(chr.gwas2),
           pos = !!as.name(pos.gwas2),
           p = !!as.name(p.gwas2)) %>%
    dplyr::select(rsid, chrom, pos, p)

df.snploc <- fread(file.snploc, sep=delim_map[delim.snploc], data.table=F, nThread=n_thread) %>%
    rename(SNP = !!as.name(snp.snploc),
           CHR = !!as.name(chr.snploc),
           START = !!as.name(start.snploc),
           END = !!as.name(end.snploc))


######

idx <- 2 # Loop this. 1:nrow(df.loc)

snp1 <- df.snploc[idx, ]$SNP
chr <- df.snploc[idx, ]$CHR
start_pos <- df.snploc[idx, ]$START
end_pos <- df.snploc[idx, ]$END

## Subset GWAS1 and compute R2 between selected SNP
df.gwas1.sub <- df.gwas1 %>%
    dplyr::filter(chrom == chr & pos >= start_pos & pos <= end_pos) %>%
    dplyr::select(chrom, pos, rsid, p)

df.gwas1.ld <- data.frame()
for (idx1 in 1:nrow(df.gwas1.sub)){ # Parallelize this
    ld.temp <- compute_snp_ld(snp1 = snp1, 
                              snp2 = df.gwas1.sub[idx1, ]$rsid, 
                              plink_exe = exe_plink, 
                              file_reference = reference)
    df.gwas1.ld <- rbind(df.gwas1.ld, ld.temp)
}

df.gwas1.sub <- merge(df.gwas1.sub, df.gwas1.ld %>% dplyr::select(-ref_snp),
                      by.x = "rsid", by.y = "other_snp", all.x = T)

## Subset GWAS2 and compute R2 between selected SNP
df.gwas2.sub <- df.gwas2 %>%
    dplyr::filter(chrom == chr & pos >= start_pos & pos <= end_pos) %>%
    dplyr::select(chrom, pos, rsid, p)

if (snp1 %in% df.gwas2.sub$rsid){
    snp2 <- snp1
} else{
    snp2 <- (df.gwas2.sub %>% arrange(p))[1, ]$rsid
}

df.gwas2.ld <- data.frame()
for (idx1 in 1:nrow(df.gwas2.sub)){ # Parallelize this
    ld.temp <- compute_snp_ld(snp1 = snp2, 
                              snp2 = df.gwas2.sub[idx1, ]$rsid, 
                              plink_exe = exe_plink, 
                              file_reference = reference)
    df.gwas2.ld <- rbind(df.gwas2.ld, ld.temp)
}

df.gwas2.sub <- merge(df.gwas2.sub, df.gwas2.ld %>% dplyr::select(-ref_snp),
                      by.x = "rsid", by.y = "other_snp", all.x = T)

## Make loc object for GWAS1 and GWAS2
if (require(EnsDb.Hsapiens.v75)){
    loc1 <- locus(data = df.gwas1.sub,
                  seqname = chr,
                  xrange = c(start_pos, end_pos),
                  ens_db = "EnsDb.Hsapiens.v75",
                  LD = "R2",
                  index_snp = snp1)
}
if (require(EnsDb.Hsapiens.v75)){
    loc2 <- locus(data = df.gwas2.sub,
                  seqname = chr,
                  xrange = c(start_pos, end_pos),
                  ens_db = "EnsDb.Hsapiens.v75",
                  LD = "R2",
                  index_snp = snp2)
}

## 
pdf(paste0(dir_out, "locuszoom.trait1_", trait.gwas1, ".trait2_", trait.gwas2,
           ".snp1_", snp1, ".snp2_", snp2, ".pdf"), 
    width = 6, height = 8)
order_partition <- set_layers(2)
if (snp1 == snp2){
    scatter_plot(loc1, labels = c("index"), 
                 border = FALSE, xtick = FALSE, legend_pos = "topright")
    scatter_plot(loc2, labels = c("index"),
                 border = FALSE, xtick = FALSE, legend_pos = "topright")
} else{
    if (snp2 %in% loc1$data$rsid){
        scatter_plot(loc1, labels = c("index", snp2),
                     border = FALSE, xtick = FALSE, legend_pos = "topright")
    } else{
        scatter_plot(loc1, labels = c("index"),
                     border = FALSE, xtick = FALSE, legend_pos = "topright")
    }
    if (snp1 %in% loc2$data$rsid){
        scatter_plot(loc2, labels = c("index", snp1),
                     border = FALSE, xtick = FALSE, legend_pos = "topright")
    } else{
        scatter_plot(loc2, labels = c("index"),
                     border = FALSE, xtick = FALSE, legend_pos = "topright")
    }
    
}

genetracks(loc1, border = FALSE,
           filter_gene_biotype = "protein_coding",
           gene_col = "grey", exon_col = "blue", exon_border = "darkgrey")

par(order_partition)
dev.off()
