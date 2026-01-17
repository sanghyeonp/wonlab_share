# 패키지
# install.packages(c("dplyr", "stringr"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("biomaRt")

library(dplyr)
library(stringr)
library(biomaRt)

query_gene_annotations_grch37_full <- function(inputs, filter_protein_coding = TRUE) {
    # Ensembl GRCh37 전용
    ens_grch37 <- useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl",
        host = "https://grch37.ensembl.org"
    )

    classify_inputs <- function(x) {
        tibble(input = x) %>%
            mutate(
            input_trim = str_replace(input, "\\.\\d+$", ""),  # ENSG 버전 제거
            type = case_when(
                str_detect(input_trim, "^ENSG\\d+$") ~ "ensembl_gene_id",
                str_detect(input_trim, "^\\d+$")      ~ "entrezgene_id",
                TRUE                                  ~ "symbol_or_synonym"
            )
            )
    }


    attrs <- c(
        "ensembl_gene_id", "hgnc_symbol", "entrezgene_id",
        "external_gene_name", "external_synonym",
        "gene_biotype", "chromosome_name", "start_position", "end_position", "strand"
    )

    # 2-1) Ensembl Gene ID로 바로 매핑
    query_by_ensembl <- function(ids, mart = ens_grch37) {
        if (length(ids) == 0) return(tibble())
        getBM(
            attributes = attrs,
            filters    = "ensembl_gene_id",
            values     = ids,
            mart       = mart
        ) %>% mutate(source_filter = "ensembl_gene_id")
    }

    # 2-2) Entrez Gene ID로 매핑
    # biomaRt에서는 보통 'entrezgene_id' 속성을 사용합니다. :contentReference[oaicite:3]{index=3}
    query_by_entrez <- function(ids, mart = ens_grch37) {
        if (length(ids) == 0) return(tibble())
        getBM(
            attributes = attrs,
            filters    = "entrezgene_id",
            values     = as.character(ids),
            mart       = mart
        ) %>% mutate(source_filter = "entrezgene_id")
    }

    # 2-3) HGNC 공식 심볼로 매핑
    query_by_symbol <- function(symbols, mart = ens_grch37) {
        if (length(symbols) == 0) return(tibble())
        getBM(
            attributes = attrs,
            filters    = "hgnc_symbol",
            values     = symbols,
            mart       = mart
        ) %>% mutate(source_filter = "hgnc_symbol")
    }

    # 2-4) 동의어로 매핑 (external_synonym)
    query_by_synonym <- function(syns, mart = ens_grch37) {
        if (length(syns) == 0) return(tibble())
        getBM(
            attributes = attrs,
            filters    = "external_synonym",
            values     = syns,
            mart       = mart
        ) %>% mutate(source_filter = "external_synonym")
    }

    # 1) 표준 스키마
    gene_table_schema <- function() {
        list(
            ensembl_gene_id     = character(),
            hgnc_symbol         = character(),
            entrezgene_id       = character(),   # 숫자처럼 보여도 캐릭터로 고정 권장
            external_gene_name  = character(),
            external_synonym    = character(),
            gene_biotype        = character(),
            chromosome_name     = character(),
            start_position      = integer(),
            end_position        = integer(),
            strand              = integer()
        )
    }

    # 2) 한 방에 타입/열 스키마 통일하는 함수
    normalize_gene_table <- function(df, schema = gene_table_schema(), drop_extras = TRUE) {
        # df가 NULL 이거나 0행인 경우도 안전하게 처리
        if (is.null(df) || nrow(df) == 0) {
            out <- lapply(schema, function(x) x)  # 0행, 지정 타입 유지
            return(as_tibble(out))
        }

        # 누락 열 채우기
        missing_cols <- setdiff(names(schema), names(df))
        if (length(missing_cols) > 0) {
            for (nm in missing_cols) {
            tgt <- typeof(schema[[nm]])
            df[[nm]] <- switch(
                tgt,
                "character" = as.character(NA),
                "integer"   = as.integer(NA),
                "double"    = as.double(NA),
                "logical"   = as.logical(NA),
                as.character(NA) # 기본
            )
            }
        }

        # 타입 캐스팅
        for (nm in names(schema)) {
            tgt <- typeof(schema[[nm]])
            if (!nm %in% names(df)) next
            # factor 방지: 항상 character로 먼저 펴기
            if (tgt == "character") df[[nm]] <- as.character(df[[nm]])
            if (tgt == "integer")   df[[nm]] <- suppressWarnings(as.integer(df[[nm]]))
            if (tgt == "double")    df[[nm]] <- suppressWarnings(as.double(df[[nm]]))
            if (tgt == "logical")   df[[nm]] <- as.logical(df[[nm]])
        }

        # 열 정렬 및 불필요 열 처리
        if (drop_extras) {
            df <- df[, names(schema), drop = FALSE]
        } else {
            # 표준 열을 앞으로
            std <- df[, names(schema), drop = FALSE]
            extra <- df[, setdiff(names(df), names(schema)), drop = FALSE]
            df <- bind_cols(std, extra)
        }

        as_tibble(df)
    }

    # 심볼로 검색
    res_symbol <- getBM(
        attributes = attrs,
        filters = "hgnc_symbol",
        values = inputs,
        mart = ens_grch37
    ) %>%
        normalize_gene_table()

    # Ensembl Gene ID로 검색
    res_ensg <- getBM(
        attributes = attrs,
        filters = "ensembl_gene_id",
        values = inputs,
        mart = ens_grch37
    ) %>%
        normalize_gene_table()

    # Entrez ID로 검색
    res_entrez <- getBM(
        attributes = attrs,
        filters = "entrezgene_id",
        values = inputs,
        mart = ens_grch37
    ) %>%
        normalize_gene_table()

    # Synonym으로 검색
    res_syn <- getBM(
        attributes = attrs,
        filters = "external_synonym",
        values = inputs,
        mart = ens_grch37
    ) %>%
        normalize_gene_table()

    # 통합
    filter_valid_chromosomes <- function(df) {
        valid_chroms <- c(as.character(1:22), "X", "Y")
        df %>% filter(chromosome_name %in% valid_chroms)
        }

    res_all <- dplyr::bind_rows(res_symbol, res_ensg, res_entrez, res_syn) %>%
        dplyr::distinct()  %>%
        filter_valid_chromosomes()
    
    if (filter_protein_coding) {
        res_all <- res_all %>% filter(gene_biotype == "protein_coding")
    }

    return(res_all)
}

query_gene_annotations_grch37.result_modification <- function(inputs, ref_gene){
    df <- query_gene_annotations_grch37_full(inputs)
    df1 <- df %>%
        dplyr::select(!!as.name(ref_gene), ensembl_gene_id) %>%
        distinct(!!as.name(ref_gene), ensembl_gene_id, .keep_all = TRUE)

    df2 <- df %>%
        dplyr::select(-ensembl_gene_id) %>%
        group_by(!!as.name(ref_gene)) %>%
        # Merge all multiple rows by ;, keeping the order
        reframe(across(everything(), ~ paste(unique(na.omit(.x)), collapse = ";"))) %>%
        ungroup()
}