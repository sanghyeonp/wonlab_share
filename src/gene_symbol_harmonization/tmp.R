suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

harmonize_human_genes <- function(df, gene_col = "gene") {
  x <- df %>%
    mutate(
      gene_input = .data[[gene_col]] %>%
        as.character() %>%
        str_trim() %>%
        str_replace("\\.\\d+$", "") # drop version suffix if present
    )

  is_ensg <- str_detect(x$gene_input, "^ENSG\\d+$")

  # Map ENSG -> SYMBOL, ENTREZ
  map_from_ensg <- tibble(gene_input = unique(x$gene_input[is_ensg])) %>%
    mutate(
      ENSEMBL = gene_input,
      SYMBOL = mapIds(org.Hs.eg.db, ENSEMBL, "SYMBOL", "ENSEMBL", multiVals = "first"),
      ENTREZID = mapIds(org.Hs.eg.db, ENSEMBL, "ENTREZID", "ENSEMBL", multiVals = "first"),
      mapping_source = "ENSEMBL"
    )

  # Map ALIAS/SYMBOL-like -> ENTREZ -> SYMBOL, ENSEMBL
  # First try ALIAS mapping (covers synonyms + old symbols); if NA, try SYMBOL directly.
  keys0 <- unique(x$gene_input[!is_ensg])

  entrez_from_alias <- mapIds(org.Hs.eg.db, keys0, "ENTREZID", "ALIAS", multiVals = "first")
  entrez_from_symbol <- mapIds(org.Hs.eg.db, keys0, "ENTREZID", "SYMBOL", multiVals = "first")

  entrez_final <- ifelse(is.na(entrez_from_alias), entrez_from_symbol, entrez_from_alias)
  source_final <- ifelse(!is.na(entrez_from_alias), "ALIAS", ifelse(!is.na(entrez_from_symbol), "SYMBOL", NA))

  map_from_names <- tibble(gene_input = keys0) %>%
    mutate(
      ENTREZID = unname(entrez_final),
      SYMBOL = mapIds(org.Hs.eg.db, ENTREZID, "SYMBOL", "ENTREZID", multiVals = "first"),
      ENSEMBL = mapIds(org.Hs.eg.db, ENTREZID, "ENSEMBL", "ENTREZID", multiVals = "first"),
      mapping_source = source_final
    )

  map_tbl <- bind_rows(map_from_ensg, map_from_names) %>%
    distinct(gene_input, .keep_all = TRUE) %>%
    mutate(
      mapping_status = case_when(
        is.na(ENSEMBL) & is.na(SYMBOL) ~ "unmapped",
        !is.na(ENSEMBL) ~ "mapped_to_ensembl",
        !is.na(SYMBOL) ~ "mapped_to_symbol_only",
        TRUE ~ "partial"
      ),
      gene_key = coalesce(ENSEMBL, SYMBOL, gene_input) # canonical key
    )

  x %>% left_join(map_tbl, by = "gene_input")
}

# Example:
# res_h <- harmonize_human_genes(res, gene_col="Gene")
# res_h %>% count(mapping_status)
# res_h %>% select(Gene, gene_input, SYMBOL, ENSEMBL, gene_key, mapping_source, mapping_status) %>% head()

#####################
df_example <- tibble::tribble(
    ~method,   ~gene,
    "FUMA",    "NFKB1",
    "MAGMA",   "NF-KB1",          # messy alias-like string (often appears in text)
    "TWAS",    "ENSG00000109320", # NFKB1 (example ENSG, may map depending on annotation release)
    "PWAS",    "TP53",
    "FUMA",    "P53",             # alias
    "MAGMA",   "ENSG00000141510.18", # TP53 with version suffix
    "TWAS",    "PARK2",           # old symbol often used for PRKN
    "MAGMA",   "PRKN",            # approved symbol
    "FUMA",    "C14orf166",       # old symbol for RTRAF (commonly updated)
    "TWAS",    "RTRAF",
    "PWAS",    "  IL6  ",         # leading/trailing whitespace
    "FUMA",    "ENST00000361390", # transcript ID (should remain unmapped in ENSG-only logic)
    "MAGMA",   "XYZNOTAGENE"      # unmapped
)

df_example


#####################
combined <- bind_rows(
  fuma %>% mutate(method = "FUMA"),
  magma %>% mutate(method = "MAGMA"),
  twas %>% mutate(method = "TWAS")
) %>%
  harmonize_human_genes(gene_col = "gene") %>%
  mutate(hit = TRUE) %>%
  distinct(method, gene_key, .keep_all = TRUE)

wide_hits <- combined %>%
  select(gene_key, SYMBOL, ENSEMBL, method, hit) %>%
  tidyr::pivot_wider(names_from = method, values_from = hit, values_fill = FALSE)

wide_hits
