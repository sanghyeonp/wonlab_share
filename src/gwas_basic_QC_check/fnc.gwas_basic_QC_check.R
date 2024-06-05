# Sanghyeon Park
# 2024.06.04
# Project: Atherosclerosis GWAS

# Workflow
# 1. Subset specified column names
# 2. Generate SNP.tmp column
# 3. Check missing values
# 4. Check invalid values
# 5. Check SNP duplicates
# 6. Check multi-allelic SNPs
# 7. Check INDEL
# 8. Check monomorphic SNPs
# 9. Check non-automsomal SNPs
# 10. Check SNP AF < threshold
# 11. Check SNP INFO < threshold
# 12. Additional criteria given

#### Version control (Current version = 2)
# v1: Initial version
# v2: Allow input of addtional criteria by the users
####

library(dplyr)

### UTILITY
map_delim <- c("txt" = "\t", "tsv" = "\t", "csv" = ",")

###
basic_GWAS_filter_criteria <- function(df_gwas, 
                                    snp_col=NA, chr_col=NA, pos_col=NA, a1_col=NA, a2_col=NA, af_col=NA, 
                                    effect_col=NA, se_col=NA, info_col=NA,
                                    threshold.AF=0.01, threshold.INFO=0.9, prefix="", 
                                    save.snp_count=NA, save.snp_list=NA,
                                    additional_columns=NA, additional_criteria=NA){
    ### Subset specified columns
    cat(":: Subset specified columns ::\n")

    potential_col <- c(snp_col, chr_col, pos_col, a1_col, a2_col, af_col, effect_col, se_col)

    if (!is.na(additional_columns)) {
        split_cols <- strsplit(additional_columns, ";")[[1]]
        # Assign each part to a variable dynamically
        for (i in seq_along(split_cols)) {
            assign(paste0("additional_col", i), split_cols[i])
            potential_col <- c(potential_col, split_cols[i])
        }
    }

    valid_col <- potential_col[which(!is.na(potential_col))]

    df <- df_gwas %>%
        dplyr::select(any_of(valid_col)) %>%
        mutate(row_index = 1:nrow(df_gwas))

    cat(paste0("\t[N SNP] Initial: ", prettyNum(nrow(df), big.mark = ",", scientific = FALSE), "\n"))
    # Make SNP column if not present
    if (is.na(snp_col)) {
        df <- mutate(df, SNP.tmp = paste(!!as.name(chr_col), !!as.name(pos_col), 
                                    pmin(!!as.name(a1_col), !!as.name(a2_col)),
                                    pmax(!!as.name(a1_col), !!as.name(a2_col)), sep=":"))
    } else{
        df <- df %>% 
            mutate(SNP.tmp = ifelse(!is.na(!!as.name(snp_col)), !!as.name(snp_col),
                                                paste(!!as.name(chr_col), !!as.name(pos_col), 
                                                    pmin(!!as.name(a1_col), !!as.name(a2_col)),
                                                    pmax(!!as.name(a1_col), !!as.name(a2_col)), sep=":")))
    }

    ### CHECK: missing values
    cat("\n:: CHECK missing value ::\n")
    nsnp.missing.chr <- NA; row.missing.chr <- c(); df.missing.chr <- data.frame()
    nsnp.missing.pos <- NA; row.missing.pos <- c(); df.missing.pos <- data.frame()
    nsnp.missing.a1 <- NA; row.missing.a1 <- c(); df.missing.a1 <- data.frame()
    nsnp.missing.a2 <- NA; row.missing.a2 <- c(); df.missing.a2 <- data.frame()
    nsnp.missing.af <- NA; row.missing.af <- c(); df.missing.af <- data.frame()
    nsnp.missing.effect <- NA; row.missing.effect <- c(); df.missing.effect <- data.frame()
    nsnp.missing.se <- NA; row.missing.se <- c(); df.missing.se <- data.frame()
    nsnp.missing.info <- NA; row.missing.info <- c(); df.missing.info <- data.frame()
    nsnp.missing.combined <- NA; row.missing.combined <- c(); df.missing.combined <- data.frame()
    # missing Chromosome
    if(!is.na(chr_col)){
        df.tmp <- dplyr::filter(df, is.na(!!as.name(chr_col)))
        nsnp.missing.chr <- nrow(df.tmp); row.missing.chr <- df.tmp$row_index
        df.missing.chr <- data.frame(row_index = row.missing.chr,
                                    criteria = rep("Missing CHR", nsnp.missing.chr))
    }
    # missing Position
    if(!is.na(pos_col)){
        df.tmp <- dplyr::filter(df, is.na(!!as.name(pos_col)))
        nsnp.missing.pos <- nrow(df.tmp); row.missing.pos <- df.tmp$row_index
        df.missing.pos <- data.frame(row_index = row.missing.pos,
                                    criteria = rep("Missing POS", nsnp.missing.pos))
    }
    # missing Allele1
    if(!is.na(a1_col)){
        df.tmp <- dplyr::filter(df, is.na(!!as.name(a1_col)))
        nsnp.missing.a1 <- nrow(df.tmp); row.missing.a1 <- df.tmp$row_index
        df.missing.a1 <- data.frame(row_index = row.missing.a1,
                                    criteria = rep("Missing A1", nsnp.missing.a1))
    }
    # missing Allele1
    if(!is.na(a2_col)){
        df.tmp <- dplyr::filter(df, is.na(!!as.name(a2_col)))
        nsnp.missing.a2 <- nrow(df.tmp); row.missing.a2 <- df.tmp$row_index
        df.missing.a2 <- data.frame(row_index = row.missing.a2,
                                    criteria = rep("Missing A2", nsnp.missing.a2))
    }
    # missing Allele frequency
    if(!is.na(af_col)){
        df.tmp <- dplyr::filter(df, is.na(!!as.name(af_col)))
        nsnp.missing.af <- nrow(df.tmp); row.missing.af <- df.tmp$row_index
        df.missing.af <- data.frame(row_index = row.missing.af,
                                    criteria = rep("Missing AF", nsnp.missing.af))
    }
    # missing Effect
    if(!is.na(effect_col)){
        df.tmp <- dplyr::filter(df, is.na(!!as.name(effect_col)))
        nsnp.missing.effect <- nrow(df.tmp); row.missing.effect <- df.tmp$row_index
        df.missing.effect <- data.frame(row_index = row.missing.effect,
                                    criteria = rep("Missing Effect", nsnp.missing.effect))
    }
    # missing SE
    if(!is.na(se_col)){
        df.tmp <- dplyr::filter(df, is.na(!!as.name(se_col)))
        nsnp.missing.se <- nrow(df.tmp); row.missing.se <- df.tmp$row_index
        df.missing.se <- data.frame(row_index = row.missing.se,
                                        criteria = rep("Missing SE", nsnp.missing.se))
    }
    # missing INFO
    if(!is.na(info_col)){
        df.tmp <- dplyr::filter(df, is.na(!!as.name(info_col)))
        nsnp.missing.info <- nrow(df.tmp); row.missing.info <- df.tmp$row_index
        df.missing.info <- data.frame(row_index = row.missing.info,
                                        criteria = rep("Missing INFO", nsnp.missing.info))
    }

    nsnp.missing.combined <- sum(c(nsnp.missing.chr, nsnp.missing.pos, nsnp.missing.a1, nsnp.missing.a2,
                                nsnp.missing.af, nsnp.missing.effect, nsnp.missing.se, nsnp.missing.info), na.rm=T)
    row.missing.combined <- unique(c(row.missing.chr, row.missing.pos, row.missing.a1, row.missing.a2, row.missing.af,
                                    row.missing.effect, row.missing.se, row.missing.info))
    df.missing.combined <- as.data.frame(bind_rows(df.missing.chr, df.missing.pos, 
                                    df.missing.a1, df.missing.a2, df.missing.af,
                                    df.missing.effect, df.missing.se, df.missing.info) %>%
        group_by(row_index) %>%
        summarise(criteria = paste(criteria, collapse = ";")) %>%
        arrange(row_index))

    cat(paste0("\t[N SNP] Missing CHR: ", prettyNum(nsnp.missing.chr, big.mark = ",", scientific = FALSE), "\n"))
    cat(paste0("\t[N SNP] Missing POS: ", prettyNum(nsnp.missing.pos, big.mark = ",", scientific = FALSE), "\n"))
    cat(paste0("\t[N SNP] Missing A1: ", prettyNum(nsnp.missing.a1, big.mark = ",", scientific = FALSE), "\n"))
    cat(paste0("\t[N SNP] Missing A2: ", prettyNum(nsnp.missing.a2, big.mark = ",", scientific = FALSE), "\n"))
    cat(paste0("\t[N SNP] Missing AF: ", prettyNum(nsnp.missing.af, big.mark = ",", scientific = FALSE), "\n"))
    cat(paste0("\t[N SNP] Missing Effect: ", prettyNum(nsnp.missing.effect, big.mark = ",", scientific = FALSE), "\n"))
    cat(paste0("\t[N SNP] Missing SE: ", prettyNum(nsnp.missing.se, big.mark = ",", scientific = FALSE), "\n"))
    cat(paste0("\t[N SNP] Missing INFO: ", prettyNum(nsnp.missing.info, big.mark = ",", scientific = FALSE), "\n"))
    cat(paste0("\t[N SNP] Missing Combined: ", prettyNum(nsnp.missing.combined, big.mark = ",", scientific = FALSE), "\n"))

    ### CHECK: invalid values
    cat("\n:: CHECK invalid value ::\n")
    nsnp.invalid.a1 <- NA; row.invalid.a1 <- c(); df.invalid.a1 <- data.frame()
    nsnp.invalid.a2 <- NA; row.invalid.a2 <- c(); df.invalid.a2 <- data.frame()
    nsnp.invalid.af <- NA; row.invalid.af <- c(); df.invalid.af <- data.frame()
    nsnp.invalid.effect <- NA; row.invalid.effect <- c(); df.invalid.effect <- data.frame()
    nsnp.invalid.se <- NA; row.invalid.se <- c(); df.invalid.se <- data.frame()
    nsnp.invalid.info <- NA; row.invalid.info <- c(); df.invalid.info <- data.frame()
    nsnp.invalid.combined <- NA; row.invalid.combined <- c(); df.invalid.combined <- data.frame()

    if(!is.na(a1_col)){
        df.tmp <- filter(df, grepl("[^ATCGatcg]", !!as.name(a1_col)))
        nsnp.invalid.a1 <- nrow(df.tmp); row.invalid.a1 <- df.tmp$row_index
        df.invalid.a1 <- data.frame(row_index = row.invalid.a1,
                                    criteria = rep("Invalid A1", nsnp.invalid.a1))
    }
    if(!is.na(a2_col)){
        df.tmp <- filter(df, grepl("[^ATCGatcg]", !!as.name(a2_col)))
        nsnp.invalid.a2 <- nrow(df.tmp); row.invalid.a2 <- df.tmp$row_index
        df.invalid.a2 <- data.frame(row_index = row.invalid.a2,
                                    criteria = rep("Invalid A2", nsnp.invalid.a2))
    }
    if(!is.na(af_col)){
        df.tmp <- filter(df, !!as.name(af_col) <0 | !!as.name(af_col) >1)
        nsnp.invalid.af <- nrow(df.tmp); row.invalid.af <- df.tmp$row_index
        df.invalid.af <- data.frame(row_index = row.invalid.af,
                                    criteria = rep("Invalid AF", nsnp.invalid.af))
    }
    if(!is.na(effect_col)){
        df.tmp <- filter(df, abs(!!as.name(effect_col)) == Inf)
        nsnp.invalid.effect <- nrow(df.tmp); row.invalid.effect <- df.tmp$row_index
        df.invalid.effect <- data.frame(row_index = row.invalid.effect,
                                    criteria = rep("Invalid Effect", nsnp.invalid.effect))
    }
    if(!is.na(se_col)){
        df.tmp <- filter(df, !!as.name(se_col) <=0 | !!as.name(se_col) == Inf)
        nsnp.invalid.se <- nrow(df.tmp); row.invalid.se <- df.tmp$row_index
        df.invalid.se <- data.frame(row_index = row.invalid.se,
                                    criteria = rep("Invalid SE", nsnp.invalid.se))
    }
    if(!is.na(info_col)){
        df.tmp <- filter(df, !!as.name(info_col) <0 | !!as.name(info_col) >1)
        nsnp.invalid.info <- nrow(df.tmp); row.invalid.info <- df.tmp$row_index
        df.invalid.info <- data.frame(row_index = row.invalid.info,
                                    criteria = rep("Invalid INFO", nsnp.invalid.info))
    }

    nsnp.invalid.combined <- sum(c(nsnp.invalid.a1, nsnp.invalid.a2, nsnp.invalid.af,
                                nsnp.invalid.effect, nsnp.invalid.se, nsnp.invalid.info), na.rm=T)
    row.invalid.combined <- unique(c(row.invalid.a1, row.invalid.a2, row.invalid.af, row.invalid.effect,
                                    row.invalid.se, row.invalid.info))
    df.invalid.combined <- as.data.frame(bind_rows(df.invalid.a1, df.invalid.a2, df.invalid.af,
                                                df.invalid.effect, df.invalid.se, df.invalid.info) %>%
                                            group_by(row_index) %>%
                                            summarise(criteria = paste(criteria, collapse = ";")) %>%
                                            arrange(row_index))

    cat(paste0("\t[N SNP] Invalid A1: ", prettyNum(nsnp.invalid.a1, big.mark = ",", scientific = FALSE), "\n"))
    cat(paste0("\t[N SNP] Invalid A2: ", prettyNum(nsnp.invalid.a2, big.mark = ",", scientific = FALSE), "\n"))
    cat(paste0("\t[N SNP] Invalid AF: ", prettyNum(nsnp.invalid.af, big.mark = ",", scientific = FALSE), "\n"))
    cat(paste0("\t[N SNP] Invalid Effect: ", prettyNum(nsnp.invalid.effect, big.mark = ",", scientific = FALSE), "\n"))
    cat(paste0("\t[N SNP] Invalid SE: ", prettyNum(nsnp.invalid.se, big.mark = ",", scientific = FALSE), "\n"))
    cat(paste0("\t[N SNP] Invalid INFO: ", prettyNum(nsnp.invalid.info, big.mark = ",", scientific = FALSE), "\n"))
    cat(paste0("\t[N SNP] Invalid Combined: ", prettyNum(nsnp.invalid.combined, big.mark = ",", scientific = FALSE), "\n"))

    ### CHECK: SNP duplicates
    cat("\n:: CHECK duplicates ::\n")
    nsnp.duplicate <- NA; row.duplicate <- c(); df.duplicate <- data.frame()
    df.tmp <- df[duplicated(df$SNP.tmp) | duplicated(df$SNP.tmp, fromLast = TRUE), ]
    nsnp.duplicate <- nrow(df.tmp); row.duplicate <- df.tmp$row_index
    df.duplicate <- data.frame(row_index = row.duplicate,
                            criteria = rep("Duplicates", nsnp.duplicate))

    cat(paste0("\t[N SNP] Duplicates: ", prettyNum(nsnp.duplicate, big.mark = ",", scientific = FALSE), "\n"))

    ### CHECK: multi-allelic SNP
    cat("\n:: CHECK multi-allelic ::\n")
    nsnp.multi_allelic <- NA; row.multi_allelic <- c(); df.multi_allelic <- data.frame()
    if (!is.na(a1_col) & !is.na(a2_col)){
        df.tmp <- filter(df, grepl(",", !!as.name(a1_col)) | grepl(",", !!as.name(a2_col)))
        nsnp.multi_allelic <- nrow(df.tmp); row.multi_allelic <- df.tmp$row_index
        df.multi_allelic <- data.frame(row_index = row.multi_allelic,
                                    criteria = rep("Multi-allelic", nsnp.multi_allelic))
    }
    cat(paste0("\t[N SNP] Multi-allelic: ", prettyNum(nsnp.multi_allelic, big.mark = ",", scientific = FALSE), "\n"))

    ### CHECK: INDEL
    cat("\n:: CHECK INDEL ::\n")
    nsnp.indel <- NA; row.indel <- c(); df.indel <- data.frame()
    if (!is.na(a1_col) & !is.na(a2_col)){
        df.tmp <- filter(df, nchar(!!as.name(a1_col)) >1 | nchar(!!as.name(a2_col)) >1)
        nsnp.indel <- nrow(df.tmp); row.indel <- df.tmp$row_index
        df.indel <- data.frame(row_index = row.indel,
                            criteria = rep("INDEL", nsnp.indel))
    }
    cat(paste0("\t[N SNP] INDEL: ", prettyNum(nsnp.indel, big.mark = ",", scientific = FALSE), "\n"))


    ### CHECK: Monomorphic SNP
    cat("\n:: CHECK monomorphic SNP ::\n")
    nsnp.monomorphic <- NA; row.monomorphic <- c(); df.monomorphic <- data.frame()
    if (!is.na(af_col)){
        df.tmp <- filter(df, !!as.name(af_col) == 0 | !!as.name(af_col) == 1)
        nsnp.monomorphic <- nrow(df.tmp); row.monomorphic <- df.tmp$row_index
        df.monomorphic <- data.frame(row_index = row.monomorphic,
                                    criteria = rep("Monomorphic", nsnp.monomorphic))
    }
    cat(paste0("\t[N SNP] Monomorphic SNP: ", prettyNum(nsnp.monomorphic, big.mark = ",", scientific = FALSE), "\n"))


    ### CHECK: Non-autosome
    cat("\n:: CHECK non-autosome ::\n")
    nsnp.non_autosome <- NA; row.non_autosome <- c(); df.non_autosome <- data.frame()
    if (!is.na(chr_col)){
        df.tmp <- df %>%
            mutate(!!as.name(chr_col) := as.character(!!as.name(chr_col))) %>%
            filter(!(!!as.name(chr_col) %in% as.character(1:22)))
        nsnp.non_autosome <- nrow(df.tmp); row.non_autosome <- df.tmp$row_index
        df.non_autosome <- data.frame(row_index = row.non_autosome,
                                    criteria = rep("Non-autosome", nsnp.non_autosome))
    }
    cat(paste0("\t[N SNP] Non-autosome: ", prettyNum(nsnp.non_autosome, big.mark = ",", scientific = FALSE), "\n"))


    ### CHECK: AF < threshold
    cat("\n:: CHECK AF < AF_threshold ::\n")
    cat(paste0("\t[AF threshold] ", threshold.AF, "\n"))
    nsnp.af_below_thres <- NA; row.af_below_thres <- c(); df.af_below_thres <- data.frame()
    if (!is.na(af_col) & !is.na(threshold.AF)){
        df.tmp <- filter(df, !!as.name(af_col) < threshold.AF)
        nsnp.af_below_thres <- nrow(df.tmp); row.af_below_thres <- df.tmp$row_index
        df.af_below_thres <- data.frame(row_index = row.af_below_thres,
                                        criteria = rep(paste0("AF<", threshold.AF), nsnp.af_below_thres))
    }
    cat(paste0("\t[N SNP] AF <", threshold.AF, ": ", prettyNum(nsnp.af_below_thres, big.mark = ",", scientific = FALSE), "\n"))

    ### CHECK: Imputation quality
    cat("\n:: CHECK INFO < INFO_threshold ::\n")
    cat(paste0("\t[INFO threshold] ", threshold.INFO, "\n"))
    nsnp.INFO_below_thres <- NA; row.INFO_below_thres <- c(); df.INFO_below_thres <- data.frame()
    if (!is.na(info_col) & !is.na(threshold.INFO)){
        df.tmp <- filter(df, !!as.name(info_col) < threhsold.INFO)
        nsnp.INFO_below_thres <- nrow(df.tmp); row.INFO_below_thres <- df.tmp$row_index
        df.INFO_below_thres <- data.frame(row_index = row.INFO_below_thres,
                                        criteria = rep(paste0("INFO<", threshold.INFO), nsnp.INFO_below_thres))
    }
    cat(paste0("\t[N SNP] INFO <", threshold.INFO, ": ", prettyNum(nsnp.INFO_below_thres, big.mark = ",", scientific = FALSE), "\n"))

    ################################################################################
    ### CHECK: Additional user-specified criteria
    cat("\n:: CHECK Additional user-specified criteria ::\n")

    split_criteria <- strsplit(additional_criteria, ";")[[1]]
    cat(paste0("\tNumber of user-specified criteria: ", length(split_criteria), "\n\n"))

    list.nsnp.addtional_criteria <- list(); row.additional_criteria <- list(); list.df.additional_criteria <- list()
    for (i in 1:length(split_criteria)){
        cat(paste0("\t[Additional criteria ", i, "]: ", split_criteria[i], "\n"))
        cat(paste0("\t[Column being used]: ", get(paste0("additional_col", i)), "\n"))

        df.tmp <- filter(df, eval(parse(text = split_criteria[i])))
        list.nsnp.addtional_criteria[[i]] <- nrow(df.tmp); row.additional_criteria[[i]] <- df.tmp$row_index
        list.df.additional_criteria[[i]] <- data.frame(row_index = row.additional_criteria[[i]],
                                                    criteria = rep(split_criteria[i], nrow(df.tmp)))

        cat(paste0("\t[N SNP] Additional criteria ", i, ": ", prettyNum(nrow(df.tmp), big.mark = ",", scientific = FALSE), "\n\n"))
    }

    ################################################################################
    df.filter.All <- bind_rows(df.missing.combined, df.invalid.combined, df.duplicate, df.multi_allelic, df.indel,
                            df.monomorphic, df.non_autosome, df.af_below_thres, df.INFO_below_thres) %>%
        group_by(row_index) %>%
        summarise(criteria = paste(criteria, collapse = ";")) %>%
        arrange(row_index)

    if (!is.na(additional_criteria)){
        for (i in 1:length(split_criteria)){
            df.filter.All <- bind_rows(df.filter.All, list.df.additional_criteria[[i]]) %>%
                group_by(row_index) %>%
                summarise(criteria = paste(criteria, collapse = ";")) %>%
                arrange(row_index)
        }
    }

    if (!is.na(save.snp_list)){
        if (save.snp_list == "rds"){
            saveRDS(df.filter.All, paste0("snp_list.", prefix, ".", save.snp_list))
        } else{
            write.table(df.filter.All,
                paste0("snp_list.", prefix, ".", save.snp_list),
                sep=map_delim[save.snp_list], row.names=F, quote=F)
        }
    }

    values <- c(nrow(df),
                nsnp.missing.chr, nsnp.missing.pos, nsnp.missing.a1, nsnp.missing.a2, nsnp.missing.af,
                nsnp.missing.effect, nsnp.missing.se, nsnp.missing.info, nsnp.missing.combined,
                nsnp.invalid.a1, nsnp.invalid.a2, nsnp.invalid.af, nsnp.invalid.effect, nsnp.invalid.se,
                nsnp.invalid.info, nsnp.invalid.combined, nsnp.duplicate, nsnp.multi_allelic, nsnp.indel,
                nsnp.monomorphic, nsnp.non_autosome, nsnp.af_below_thres, nsnp.INFO_below_thres)

    if (!is.na(additional_criteria)){
        for (i in 1:length(split_criteria)){
            values <- c(values, list.nsnp.addtional_criteria[[i]])
        }
    }

    values <- c(values, nrow(df.filter.All), nrow(df) - nrow(df.filter.All))

    df.nsnp_count <- data.frame(matrix(values, ncol = length(values), nrow = 1))

    column_name <- c("Initial",
                    "Missing CHR", "Missing POS", "Missing A1", "Missing A2", "Missing AF",
                    "Missing Effect", "Missing SE", "Missing INFO", "Missing Combined",
                    "Invalid A1", "Invalid A2", "Invalid AF", "Invalid Effect", "Invalid SE",
                    "Invalid INFO", "Invalid Combined", "Duplicates", "Multi-allelic", "INDEL",
                    "Monomorphic", "Non-autosome", paste0("AF<", threshold.AF), paste0("INFO<", threshold.INFO))
    if (!is.na(additional_criteria)){
        for (i in 1:length(split_criteria)){
            column_name <- c(column_name, split_criteria[i])
        }
    }
    column_name <- c(column_name, "All (unique)", "Final")

    colnames(df.nsnp_count) <- column_name

    if (!is.na(save.snp_count)){
        write.table(df.nsnp_count,
                    paste0("nsnp_count.", prefix, ".", save.snp_count),
                    sep=map_delim[save.snp_count], row.names=F, quote=F)
    }
}
