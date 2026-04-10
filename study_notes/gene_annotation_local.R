library(EnsDb.Hsapiens.v75)

ensg_list <- c()
df.gene_annot <- as.data.frame(
    genes(EnsDb.Hsapiens.v75,
          filter = GeneIdFilter(ensg_list),
          columns = c("gene_id", "gene_name", "gene_biotype",
                      "seq_name", "gene_seq_start", "gene_seq_end", "seq_strand"))
)

symbol_list <- c()
df.gene_annot <- as.data.frame(
    genes(EnsDb.Hsapiens.v75,
          filter = GeneNameFilter(symbol_list),
          columns = c("gene_id", "gene_name", "gene_biotype",
                      "seq_name", "gene_seq_start", "gene_seq_end", "seq_strand"))
)
