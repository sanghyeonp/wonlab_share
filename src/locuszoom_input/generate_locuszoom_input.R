library(data.table)
library(dplyr)
library(argparse)


# Parse command line arguments
parser <- ArgumentParser()
parser$add_argument("--file-in", dest="file_in", type = "character", required = TRUE,
                    help="Input file GWAS")
parser$add_argument("--delim", type = "character", required = TRUE,
                    help="")
parser$add_argument("--snp-col", dest="snp_col", type = "character", required = TRUE,
                    help="")
parser$add_argument("--chr-col", dest="chr_col", type = "character", required = TRUE,
                    help="")
parser$add_argument("--pos-col", dest="pos_col", type = "character", required = TRUE,
                    help="")
parser$add_argument("--a1-col", dest="a1_col", type = "character", required = TRUE,
                    help="")
parser$add_argument("--a2-col", dest="a2_col", type = "character", required = TRUE,
                    help="")
parser$add_argument("--pval-col", dest="pval_col", type = "character", required = TRUE,
                    help="")
parser$add_argument("--out-pref", dest="out_pref", type = "character", required = TRUE,
                    help="")
parser$add_argument("--thread", type = "integer", required = FALSE, default = 1, 
                    help="")


args <- parser$parse_args()
file_in <- args$file_in
delim <- args$delim
snp_col <- args$snp_col
chr_col <- args$chr_col
pos_col <- args$pos_col
a1_col <- args$a1_col
a2_col <- args$a2_col
pval_col <- args$pval_col

file_prefix <- args$out_pref
n_thread <- args$thread

###
delim_map <- c("whitespace" = " ",
                "tab" = "\t", 
                "comma" = ",")

### 1. Read GWAS
df <- fread(file_in, sep=delim_map[delim], data.table = F, nThread = n_thread)

head(df, 5)

### 2. Read 1000 Genome bim file
df.bim <- fread("/data1/sanghyeon/wonlab_contribute/combined/data_common/reference_panel/1kGp3/EUR/reference.1kG.EUR.bim", 
                sep = "\t", data.table = F, header = F, nThread = n_thread, col.names = c("CHR", "SNP", "CM", "BP", "ALT", "REF")) %>%
    dplyr::select(SNP, ALT, REF)

head(df.bim, 5)

### 2. Sort by chromosome, position
df <- df %>%
    # 필요한 column 만 남기고.
    dplyr::select(!!as.name(snp_col), !!as.name(chr_col), !!as.name(pos_col), 
                  !!as.name(a1_col), !!as.name(a2_col), !!as.name(pval_col))
# A1, A2를 1000 Genome을 기반으로 reference, alternative 구분해주고.
df1 <- merge(df, df.bim, by.x = snp_col, by.y = "SNP", all.x = T)

head(df1, 5)

# 추가 manipulation
df1 <- df1 %>%
    # Reference, alternative allele 이 없는 경우, GWAS의 A1이 alternative, A2가 reference로
    dplyr::mutate(REF = ifelse(is.na(REF), !!as.name(a2_col), REF),
                  ALT = ifelse(is.na(ALT), !!as.name(a1_col), ALT)) %>%
    # chr, pos로 sorting하고.
    dplyr::arrange(!!as.name(chr_col), !!as.name(pos_col)) %>%
    # 필요한 column만 남기고 (markername, chr, pos, ref, alt, pval).
    dplyr::select(!!as.name(chr_col), !!as.name(pos_col), REF, ALT, !!as.name(pval_col)) %>%
    # column이름 바꿔주고.
    dplyr::rename("Chr" = !!as.name(chr_col),
                "Pos" = !!as.name(pos_col),
                "P" = !!as.name(pval_col))

head(df1, 5)

### 3. Save sorted summary statistics
write.table(df1,
            paste0(file_prefix, ".tsv"),
            sep="\t", row.names = F, quote = F)
system(paste0("/data1/sanghyeon/wonlab_contribute/combined/software/htslib/bin/bgzip -c ./", file_prefix, ".tsv > ./", file_prefix, ".tsv.bgz"))

### 4. tabix파일 만들기
exe_tabix <- "/data1/sanghyeon/wonlab_contribute/combined/software/htslib/bin/tabix"
system(paste0(exe_tabix, 
            " -s1 -b 2 -e 2 --skip-lines 1 -f ./", file_prefix, ".tsv.bgz"))

### 5. 필요없는 파일 지우기
system(paste0("rm ./", file_prefix, ".tsv"))
