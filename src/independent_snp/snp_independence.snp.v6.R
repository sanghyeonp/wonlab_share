############################################################
# GWAS signal comparison by LD
#
# For each GWAS1 lead SNP:
#
# 1. Find GWAS2 variants satisfying:
#       P < p_threshold
#       within +/- window kb
#
# 2. Calculate LD using reference panel
#
# 3. Determine:
#       In_LD = TRUE
#           if >=1 candidate has R2 > r2_threshold
#
#       need_check = yes
#           if In_LD == FALSE
#           AND >=1 candidate pair has unresolved R2
#
#
# Two output files:
#
# [1] Summary
#     One row per GWAS1 lead SNP
#
# [2] Pairwise detail
#     One row per
#       GWAS1 lead SNP x GWAS2 GWS candidate SNP
############################################################

############################################################
# snp_independence.snp.v4.R의 업데이트 버전
#   - 속도 매우 향상
############################################################


############################################################
# Libraries
############################################################

library(data.table)
library(argparse)


############################################################
# Stop immediately on error
############################################################

options(
    error = function() {
        traceback(2)
        quit(
            save = "no",
            status = 1,
            runLast = FALSE
        )
    }
)


############################################################
# Parser
############################################################

parser <- argparse::ArgumentParser(
    description = ":: SNP signal comparison by P-value and LD ::",
    formatter_class = "argparse.ArgumentDefaultsHelpFormatter"
)


############################################################
# Trait names
############################################################

parser$add_argument(
    "--trait1",
    required = TRUE
)

parser$add_argument(
    "--trait2",
    required = TRUE
)


############################################################
# GWAS1 SNP list
############################################################

parser$add_argument(
    "--gwas1-snplist",
    dest = "file_gwas1_snplist",
    required = TRUE
)

parser$add_argument(
    "--delim-gwas1-snplist",
    dest = "delim_gwas1_snplist",
    required = TRUE
)

parser$add_argument(
    "--snp-col-gwas1-snplist",
    dest = "snp_gwas1_snplist",
    required = TRUE
)


############################################################
# GWAS1 summary statistics
############################################################

parser$add_argument(
    "--gwas1",
    dest = "file_gwas1",
    required = TRUE
)

parser$add_argument(
    "--delim-gwas1",
    dest = "delim_gwas1",
    required = TRUE
)

parser$add_argument(
    "--snp-col-gwas1",
    dest = "snp_gwas1",
    required = TRUE
)

parser$add_argument(
    "--chr-col-gwas1",
    dest = "chr_gwas1",
    required = TRUE
)

parser$add_argument(
    "--pos-col-gwas1",
    dest = "pos_gwas1",
    required = TRUE
)

parser$add_argument(
    "--pval-col-gwas1",
    dest = "p_gwas1",
    required = TRUE
)


############################################################
# GWAS2 summary statistics
############################################################

parser$add_argument(
    "--gwas2",
    dest = "file_gwas2",
    required = TRUE
)

parser$add_argument(
    "--delim-gwas2",
    dest = "delim_gwas2",
    required = TRUE
)

parser$add_argument(
    "--snp-col-gwas2",
    dest = "snp_gwas2",
    required = TRUE
)

parser$add_argument(
    "--chr-col-gwas2",
    dest = "chr_gwas2",
    required = TRUE
)

parser$add_argument(
    "--pos-col-gwas2",
    dest = "pos_gwas2",
    required = TRUE
)

parser$add_argument(
    "--pval-col-gwas2",
    dest = "p_gwas2",
    required = TRUE
)


############################################################
# Reference panel
############################################################

parser$add_argument(
    "--reference-panel",
    dest = "reference_panel",
    required = FALSE,
    default = "1kG"
)

parser$add_argument(
    "--genome-build",
    dest = "genome_build",
    required = FALSE,
    default = "GRCh37"
)


############################################################
# Thresholds
############################################################

parser$add_argument(
    "--r2-threshold",
    dest = "r2_threshold",
    required = TRUE,
    type = "numeric"
)

parser$add_argument(
    "--window",
    required = TRUE,
    type = "integer"
)

parser$add_argument(
    "--pval-threshold",
    dest = "p_threshold",
    required = TRUE,
    type = "numeric"
)


############################################################
# Thread
############################################################

parser$add_argument(
    "--thread",
    required = FALSE,
    default = 4,
    type = "integer"
)


############################################################
# Output
############################################################

parser$add_argument(
    "--dir-out",
    dest = "dir_out",
    required = FALSE,
    default = "NA"
)

parser$add_argument(
    "--prefix-out",
    dest = "pref_out",
    required = FALSE,
    default = "NA"
)

parser$add_argument(
    "--delim-out",
    dest = "delim_out",
    required = FALSE,
    default = "comma"
)


############################################################
# Parse
############################################################

args <- parser$parse_args()


trait1 <- args$trait1
trait2 <- args$trait2

file_gwas1_snplist <- args$file_gwas1_snplist
delim_gwas1_snplist <- args$delim_gwas1_snplist
snp_gwas1_snplist <- args$snp_gwas1_snplist

file_gwas1 <- args$file_gwas1
delim_gwas1 <- args$delim_gwas1
snp_gwas1 <- args$snp_gwas1
chr_gwas1 <- args$chr_gwas1
pos_gwas1 <- args$pos_gwas1
p_gwas1 <- args$p_gwas1

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

n_thread <- args$thread

dir_out <- args$dir_out
pref_out <- args$pref_out
delim_out <- args$delim_out


############################################################
# PLINK
############################################################

plink_exe <- paste0(
    "/data1/sanghyeon/wonlab_contribute/combined/",
    "software/plink/plink"
)


############################################################
# Reference panels
############################################################

reference_panel_list <- c(

    "1kG" = paste0(
        "/data1/sanghyeon/wonlab_contribute/combined/",
        "data_common/reference_panel/1kGp3/EUR/CHR/",
        "reference.1kG.EUR.maf_0.005.geno_0.02.chr"
    ),

    "UKB_random_10K" = paste0(
        "/data1/sanghyeon/wonlab_contribute/combined/",
        "data_common/reference_panel/UKB_random10k_noQC/CHR/",
        "UKB_random10k_noqc_chr"
    ),

    "KoGES_random_10k" = paste0(
        "/data1/sanghyeon/wonlab_contribute/combined/",
        "data_common/reference_panel/KoGES/CHR/",
        "KoGES_random10k.rsid.b37.drop_dup_chr"
    )
)


############################################################
# Delimiters
############################################################

delim_mapper <- c(
    "tab" = "\t",
    "whitespace" = " ",
    "comma" = ","
)

exe_mapper <- c(
    "tab" = "tsv",
    "whitespace" = "txt",
    "comma" = "csv"
)


############################################################
# Basic checks
############################################################

if (!(reference_panel %in% names(reference_panel_list))) {
    stop(
        "Unknown reference panel: ",
        reference_panel
    )
}


if (!(delim_out %in% names(delim_mapper))) {
    stop(
        "Invalid output delimiter."
    )
}


if (r2_threshold < 0 || r2_threshold > 1) {
    stop(
        "R2 threshold must be between 0 and 1."
    )
}


if (window <= 0) {
    stop(
        "Window must be > 0 kb."
    )
}


if (n_thread < 1) {
    stop(
        "Thread must be >= 1."
    )
}


if (!file.exists(plink_exe)) {
    stop(
        "PLINK not found: ",
        plink_exe
    )
}


############################################################
# Output directory
############################################################

if (dir_out == "NA") {
    dir_out <- "."
}


dir.create(
    dir_out,
    recursive = TRUE,
    showWarnings = FALSE
)


if (pref_out == "NA") {

    pref_out <- paste0(
        "snp_independence.",
        tools::file_path_sans_ext(
            basename(file_gwas1_snplist)
        )
    )
}


############################################################
# Utility
############################################################

normalize_chr <- function(x) {

    x <- as.character(x)

    x <- gsub(
        "^chr",
        "",
        x,
        ignore.case = TRUE
    )

    suppressWarnings(
        as.integer(x)
    )
}


############################################################
# GWAS1 lead SNP list
############################################################

message(
    "Reading GWAS1 lead SNP list..."
)


df_gwas1_snplist <- fread(

    file_gwas1_snplist,

    sep = delim_mapper[[delim_gwas1_snplist]],

    select = snp_gwas1_snplist,

    data.table = TRUE,

    nThread = n_thread,

    showProgress = FALSE
)


setnames(
    df_gwas1_snplist,
    snp_gwas1_snplist,
    "leadSNP"
)


df_gwas1_snplist <- unique(
    df_gwas1_snplist,
    by = "leadSNP"
)


############################################################
# Preserve original lead SNP order
############################################################

df_gwas1_snplist[
    ,
    input_order := .I
]


############################################################
# GWAS1 summary statistics
############################################################

message(
    "Reading GWAS1 summary statistics..."
)


df_gwas1 <- fread(

    file_gwas1,

    sep = delim_mapper[[delim_gwas1]],

    select = c(
        snp_gwas1,
        chr_gwas1,
        pos_gwas1,
        p_gwas1
    ),

    data.table = TRUE,

    nThread = n_thread,

    showProgress = FALSE
)


setnames(

    df_gwas1,

    old = c(
        snp_gwas1,
        chr_gwas1,
        pos_gwas1,
        p_gwas1
    ),

    new = c(
        "SNP",
        "CHR",
        "POS",
        "P"
    )
)


df_gwas1 <- df_gwas1[
    grepl("^rs", SNP)
]


df_gwas1[
    ,
    CHR := normalize_chr(CHR)
]

df_gwas1[
    ,
    POS := as.numeric(POS)
]

df_gwas1[
    ,
    P := as.numeric(P)
]


############################################################
# Avoid duplicated rsIDs
############################################################

setorder(
    df_gwas1,
    SNP,
    P
)


df_gwas1 <- unique(
    df_gwas1,
    by = "SNP"
)


############################################################
# Map lead SNP position
############################################################

df_gwas1_snplist <- merge(

    df_gwas1_snplist,

    df_gwas1[
        ,
        .(
            leadSNP = SNP,
            CHR,
            POS,
            P
        )
    ],

    by = "leadSNP",

    all.x = TRUE,

    sort = FALSE
)


############################################################
# Check missing lead SNP
############################################################

missing_leads <- df_gwas1_snplist[
    is.na(CHR) |
    is.na(POS),
    unique(leadSNP)
]


if (length(missing_leads) > 0) {

    stop(
        paste0(
            "Lead SNP(s) could not be mapped in GWAS1:\n",
            paste(
                missing_leads,
                collapse = ", "
            )
        )
    )
}


############################################################
# Restore input ordering
############################################################

setorder(
    df_gwas1_snplist,
    input_order
)


############################################################
# GWAS2
############################################################

message(
    "Reading GWAS2 summary statistics..."
)


df_gwas2 <- fread(

    file_gwas2,

    sep = delim_mapper[[delim_gwas2]],

    select = c(
        snp_gwas2,
        chr_gwas2,
        pos_gwas2,
        p_gwas2
    ),

    data.table = TRUE,

    nThread = n_thread,

    showProgress = FALSE
)


setnames(

    df_gwas2,

    old = c(
        snp_gwas2,
        chr_gwas2,
        pos_gwas2,
        p_gwas2
    ),

    new = c(
        "SNP",
        "CHR",
        "POS",
        "P"
    )
)


df_gwas2 <- df_gwas2[
    grepl("^rs", SNP)
]


df_gwas2[
    ,
    CHR := normalize_chr(CHR)
]

df_gwas2[
    ,
    POS := as.numeric(POS)
]

df_gwas2[
    ,
    P := as.numeric(P)
]


############################################################
# Keep GWS variants only
############################################################

df_gwas2 <- df_gwas2[
    !is.na(CHR) &
    !is.na(POS) &
    !is.na(P) &
    P < p_threshold
]


############################################################
# If duplicate rsID:
# retain smallest P
############################################################

setorder(
    df_gwas2,
    SNP,
    P
)


df_gwas2 <- unique(
    df_gwas2,
    by = "SNP"
)


message(
    "GWAS2 SNPs with P < ",
    format(
        p_threshold,
        scientific = TRUE
    ),
    ": ",
    nrow(df_gwas2)
)


############################################################
# Reference panel
############################################################

file_reference_panel <-
    reference_panel_list[[reference_panel]]


message(
    "Using reference panel: ",
    reference_panel
)


message(
    "Reference prefix: ",
    file_reference_panel
)


############################################################
# Result containers
############################################################

summary_list <- list()

detail_list <- list()

summary_index <- 1L

detail_index <- 1L


############################################################
# Main chromosome loop
############################################################

for (chr in 1:22) {


    message("")
    message(
        "========================================"
    )

    message(
        "Processing chromosome ",
        chr
    )

    message(
        "========================================"
    )


    ########################################################
    # GWAS1 lead SNPs
    ########################################################

    df_lead_chr <- df_gwas1_snplist[
        CHR == chr
    ]


    if (nrow(df_lead_chr) == 0) {

        message(
            "No GWAS1 lead SNPs. Skipping."
        )

        next
    }


    query_snps <- unique(
        df_lead_chr$leadSNP
    )


    message(
        "GWAS1 query SNPs: ",
        length(query_snps)
    )


    ########################################################
    # GWAS2 significant SNPs
    ########################################################

    df_gwas2_chr <- df_gwas2[
        CHR == chr
    ]


    message(
        "GWAS2 significant SNPs: ",
        nrow(df_gwas2_chr)
    )


    ########################################################
    # Default summary
    ########################################################

    df_summary_chr <- data.table(

        SNP_gwas1 =
            query_snps,

        In_LD =
            FALSE,

        need_check =
            "no",

        n_candidate =
            0L,

        n_r2_computed =
            0L,

        n_r2_not_computed =
            0L,

        SNP_gwas2_in_LD =
            "NA",

        SNP_gwas2_in_LD_pval =
            "NA",

        SNP_gwas2_in_LD_R2 =
            "NA",

        SNP_gwas2_not_computed =
            "NA",

        SNP_gwas2_not_computed_pval =
            "NA"
    )


    ########################################################
    # No significant SNP on chromosome
    ########################################################

    if (nrow(df_gwas2_chr) == 0) {

        summary_list[[summary_index]] <-
            df_summary_chr

        summary_index <-
            summary_index + 1L

        next
    }


    ########################################################
    # Generate candidate pairs
    ########################################################

    lead_tmp <- df_lead_chr[
        ,
        .(
            SNP_A = leadSNP,
            POS_A = POS,
            tmp_join = 1L
        )
    ]


    gwas2_tmp <- df_gwas2_chr[
        ,
        .(
            SNP_B = SNP,
            POS_B = POS,
            P = P,
            tmp_join = 1L
        )
    ]


    df_candidate_pairs <- merge(

        lead_tmp,

        gwas2_tmp,

        by = "tmp_join",

        allow.cartesian = TRUE,

        sort = FALSE
    )


    df_candidate_pairs[
        ,
        tmp_join := NULL
    ]


    ########################################################
    # Restrict to +/- window
    ########################################################

    window_bp <-
        window * 1000


    df_candidate_pairs <-
        df_candidate_pairs[
            abs(POS_B - POS_A) <= window_bp
        ]


    ########################################################
    # No candidate within window
    ########################################################

    if (nrow(df_candidate_pairs) == 0) {

        message(
            "No GWAS2 GWS SNPs within +/- ",
            window,
            " kb."
        )


        summary_list[[summary_index]] <-
            df_summary_chr

        summary_index <-
            summary_index + 1L

        next
    }


    ########################################################
    # Distance
    ########################################################

    df_candidate_pairs[
        ,
        distance_bp :=
            abs(POS_B - POS_A)
    ]


    message(
        "Candidate pairs within +/- ",
        window,
        " kb: ",
        nrow(df_candidate_pairs)
    )


    ########################################################
    # Initialize LD information
    ########################################################

    df_candidate_pairs[
        ,
        `:=`(
            R2 = as.numeric(NA),
            r2_computed = FALSE
        )
    ]


    ########################################################
    # Same SNP:
    #
    # R2 = 1
    ########################################################

    df_candidate_pairs[
        SNP_A == SNP_B,
        `:=`(
            R2 = 1,
            r2_computed = TRUE
        )
    ]


    ########################################################
    # Reference files
    ########################################################

    ref_prefix <- paste0(
        file_reference_panel,
        chr
    )


    bim_file <- paste0(
        ref_prefix,
        ".bim"
    )

    bed_file <- paste0(
        ref_prefix,
        ".bed"
    )

    fam_file <- paste0(
        ref_prefix,
        ".fam"
    )


    if (
        !file.exists(bim_file) ||
        !file.exists(bed_file) ||
        !file.exists(fam_file)
    ) {

        stop(
            "Reference panel files missing: ",
            ref_prefix
        )
    }


    ########################################################
    # Reference SNP IDs
    ########################################################

    ref_snps <- fread(

        bim_file,

        header = FALSE,

        select = 2,

        data.table = TRUE,

        nThread = n_thread,

        showProgress = FALSE
    )


    setnames(
        ref_snps,
        names(ref_snps)[1],
        "SNP"
    )


    ref_snp_vector <-
        unique(ref_snps$SNP)


    ########################################################
    # Is each SNP in reference?
    ########################################################

    df_candidate_pairs[
        ,
        lead_in_reference :=
            SNP_A %in% ref_snp_vector
    ]


    df_candidate_pairs[
        ,
        gwas2_in_reference :=
            SNP_B %in% ref_snp_vector
    ]


    ########################################################
    # Pairs eligible for PLINK
    ########################################################

    df_pairs_to_plink <- df_candidate_pairs[
        SNP_A != SNP_B &
        lead_in_reference == TRUE &
        gwas2_in_reference == TRUE
    ]


    message(
        "Pairs eligible for PLINK: ",
        nrow(df_pairs_to_plink)
    )


    message(
        "Pairs missing >=1 SNP in reference: ",
        df_candidate_pairs[
            SNP_A != SNP_B &
            (
                lead_in_reference == FALSE |
                gwas2_in_reference == FALSE
            ),
            .N
        ]
    )


    ########################################################
    # Batch PLINK
    ########################################################

    if (nrow(df_pairs_to_plink) > 0) {


        ####################################################
        # Query SNPs
        ####################################################

        query_snps_for_plink <-
            unique(
                df_pairs_to_plink$SNP_A
            )


        ####################################################
        # Only variants needed for calculation
        ####################################################

        extract_snps <- unique(
            c(
                df_pairs_to_plink$SNP_A,
                df_pairs_to_plink$SNP_B
            )
        )


        ####################################################
        # Temporary files
        ####################################################

        query_file <- tempfile(
            pattern = paste0(
                "query_chr",
                chr,
                "_"
            ),
            fileext = ".txt"
        )


        extract_file <- tempfile(
            pattern = paste0(
                "extract_chr",
                chr,
                "_"
            ),
            fileext = ".txt"
        )


        plink_out_prefix <- tempfile(
            pattern = paste0(
                "batch_ld_chr",
                chr,
                "_"
            )
        )


        ####################################################
        # Write SNP files
        ####################################################

        fwrite(

            data.table(
                SNP = query_snps_for_plink
            ),

            query_file,

            col.names = FALSE,

            quote = FALSE
        )


        fwrite(

            data.table(
                SNP = extract_snps
            ),

            extract_file,

            col.names = FALSE,

            quote = FALSE
        )


        ####################################################
        # IMPORTANT
        #
        # --ld-window-r2 0
        #
        # We want ALL calculable R2 values.
        #
        # R2 threshold is applied later in R.
        ####################################################

        plink_args <- c(

            "--bfile",
            ref_prefix,

            "--extract",
            extract_file,

            "--r2",

            "--ld-snp-list",
            query_file,

            "--ld-window",
            "99999999",

            "--ld-window-kb",
            as.character(window),

            "--ld-window-r2",
            "0",

            "--threads",
            as.character(n_thread),

            "--out",
            plink_out_prefix
        )


        message(
            "Running batch PLINK LD calculation..."
        )


        plink_status <- system2(
            command = plink_exe,
            args = plink_args
        )


        if (plink_status != 0) {

            stop(
                paste0(
                    "PLINK failed on chromosome ",
                    chr,
                    ". Check ",
                    plink_out_prefix,
                    ".log"
                )
            )
        }


        ####################################################
        # Read PLINK LD output
        ####################################################

        ld_file <- paste0(
            plink_out_prefix,
            ".ld"
        )


        if (file.exists(ld_file)) {


            df_ld <- fread(
                ld_file,
                data.table = TRUE,
                showProgress = FALSE
            )


            required_cols <- c(
                "SNP_A",
                "SNP_B",
                "R2"
            )


            if (
                nrow(df_ld) > 0 &&
                all(required_cols %in% names(df_ld))
            ) {


                ################################################
                # Keep finite R2 only
                ################################################

                df_ld <- df_ld[
                    is.finite(R2)
                ]


                ################################################
                # Keep unique SNP pairs
                ################################################

                df_ld <- unique(

                    df_ld[
                        ,
                        .(
                            SNP_A,
                            SNP_B,
                            R2
                        )
                    ],

                    by = c(
                        "SNP_A",
                        "SNP_B"
                    )
                )


                ################################################
                # Rename before merge
                ################################################

                setnames(
                    df_ld,
                    "R2",
                    "R2_plink"
                )


                ################################################
                # Merge into complete expected pair table
                #
                # If a pair is absent from PLINK output,
                # R2_plink remains NA.
                ################################################

                df_candidate_pairs <- merge(

                    df_candidate_pairs,

                    df_ld,

                    by = c(
                        "SNP_A",
                        "SNP_B"
                    ),

                    all.x = TRUE,

                    sort = FALSE
                )


                ################################################
                # Populate calculated R2
                ################################################

                df_candidate_pairs[
                    SNP_A != SNP_B &
                    !is.na(R2_plink),
                    `:=`(
                        R2 = R2_plink,
                        r2_computed = TRUE
                    )
                ]


                df_candidate_pairs[
                    ,
                    R2_plink := NULL
                ]
            }
        }


        ####################################################
        # Delete temporary files
        ####################################################

        temp_files <- unique(
            c(
                query_file,
                extract_file,
                Sys.glob(
                    paste0(
                        plink_out_prefix,
                        ".*"
                    )
                )
            )
        )


        unlink(
            temp_files,
            force = TRUE
        )
    }


    ########################################################
    # Pair-level LD classification
    ########################################################

    df_candidate_pairs[
        ,
        In_LD_pair := as.logical(NA)
    ]


    df_candidate_pairs[
        r2_computed == TRUE &
        !is.na(R2),
        In_LD_pair :=
            R2 > r2_threshold
    ]


    ########################################################
    # Reason / status
    #
    # Useful for pairwise-detail file
    ########################################################

    df_candidate_pairs[
        ,
        r2_status := "not_computed"
    ]


    df_candidate_pairs[
        SNP_A == SNP_B,
        r2_status := "same_SNP_R2_1"
    ]


    df_candidate_pairs[
        SNP_A != SNP_B &
        (
            lead_in_reference == FALSE |
            gwas2_in_reference == FALSE
        ),
        r2_status := "not_computed_missing_reference"
    ]


    df_candidate_pairs[
        SNP_A != SNP_B &
        lead_in_reference == TRUE &
        gwas2_in_reference == TRUE &
        r2_computed == FALSE,
        r2_status := "not_computed_by_PLINK"
    ]


    df_candidate_pairs[
        r2_computed == TRUE &
        SNP_A != SNP_B,
        r2_status := "computed"
    ]


    ########################################################
    # Sort candidates
    #
    # Important because semicolon-separated fields
    # should have matching order.
    ########################################################

    setorder(
        df_candidate_pairs,
        SNP_A,
        P,
        SNP_B
    )


    ########################################################
    # Save pairwise detail
    ########################################################

    df_detail_chr <- df_candidate_pairs[
        ,
        .(
            SNP_gwas1 =
                SNP_A,

            POS_gwas1 =
                POS_A,

            SNP_gwas2 =
                SNP_B,

            POS_gwas2 =
                POS_B,

            P_gwas2 =
                P,

            distance_bp =
                distance_bp,

            lead_in_reference =
                lead_in_reference,

            gwas2_in_reference =
                gwas2_in_reference,

            R2 =
                R2,

            r2_computed =
                r2_computed,

            In_LD_pair =
                In_LD_pair,

            r2_status =
                r2_status
        )
    ]


    detail_list[[detail_index]] <-
        df_detail_chr


    detail_index <-
        detail_index + 1L


    ########################################################
    # Summary per GWAS1 SNP
    ########################################################

    df_pair_summary <- df_candidate_pairs[
        ,
        {

            ################################################
            # Confirmed LD hits
            ################################################

            idx_ld <- which(
                r2_computed == TRUE &
                !is.na(R2) &
                R2 > r2_threshold
            )


            ################################################
            # Unresolved pairs
            ################################################

            idx_unresolved <- which(
                r2_computed == FALSE
            )


            ################################################
            # Is any confirmed LD hit present?
            ################################################

            confirmed_ld <-
                length(idx_ld) > 0


            ################################################
            # Need manual check?
            ################################################

            need_manual_check <-
                !confirmed_ld &&
                length(idx_unresolved) > 0


            ################################################
            # Output
            ################################################

            list(

                In_LD =
                    confirmed_ld,

                need_check =
                    ifelse(
                        need_manual_check,
                        "yes",
                        "no"
                    ),

                n_candidate =
                    .N,

                n_r2_computed =
                    sum(
                        r2_computed
                    ),

                n_r2_not_computed =
                    sum(
                        !r2_computed
                    ),


                ################################################
                # Confirmed LD variants
                ################################################

                SNP_gwas2_in_LD =
                    if (
                        length(idx_ld) > 0
                    ) {
                        paste(
                            SNP_B[idx_ld],
                            collapse = ";"
                        )
                    } else {
                        "NA"
                    },


                SNP_gwas2_in_LD_pval =
                    if (
                        length(idx_ld) > 0
                    ) {
                        paste(
                            P[idx_ld],
                            collapse = ";"
                        )
                    } else {
                        "NA"
                    },


                SNP_gwas2_in_LD_R2 =
                    if (
                        length(idx_ld) > 0
                    ) {
                        paste(
                            R2[idx_ld],
                            collapse = ";"
                        )
                    } else {
                        "NA"
                    },


                ################################################
                # R2 not computed variants
                ################################################

                SNP_gwas2_not_computed =
                    if (
                        length(idx_unresolved) > 0
                    ) {
                        paste(
                            SNP_B[idx_unresolved],
                            collapse = ";"
                        )
                    } else {
                        "NA"
                    },


                SNP_gwas2_not_computed_pval =
                    if (
                        length(idx_unresolved) > 0
                    ) {
                        paste(
                            P[idx_unresolved],
                            collapse = ";"
                        )
                    } else {
                        "NA"
                    }
            )
        },

        by = .(
            SNP_gwas1 = SNP_A
        )
    ]


    ########################################################
    # Merge with default results
    #
    # SNPs with no GWS candidate remain:
    #
    # In_LD = FALSE
    # need_check = no
    # n_candidate = 0
    ########################################################

    idx <- match(
        df_pair_summary$SNP_gwas1,
        df_summary_chr$SNP_gwas1
    )


    df_summary_chr[
        idx,
        `:=`(

            In_LD =
                df_pair_summary$In_LD,

            need_check =
                df_pair_summary$need_check,

            n_candidate =
                df_pair_summary$n_candidate,

            n_r2_computed =
                df_pair_summary$n_r2_computed,

            n_r2_not_computed =
                df_pair_summary$n_r2_not_computed,

            SNP_gwas2_in_LD =
                df_pair_summary$SNP_gwas2_in_LD,

            SNP_gwas2_in_LD_pval =
                df_pair_summary$SNP_gwas2_in_LD_pval,

            SNP_gwas2_in_LD_R2 =
                df_pair_summary$SNP_gwas2_in_LD_R2,

            SNP_gwas2_not_computed =
                df_pair_summary$SNP_gwas2_not_computed,

            SNP_gwas2_not_computed_pval =
                df_pair_summary$SNP_gwas2_not_computed_pval
        )
    ]


    ########################################################
    # Console summary
    ########################################################

    message(
        "In_LD = TRUE: ",
        sum(
            df_summary_chr$In_LD
        ),
        "/",
        nrow(df_summary_chr)
    )


    message(
        "need_check = yes: ",
        sum(
            df_summary_chr$need_check == "yes"
        ),
        "/",
        nrow(df_summary_chr)
    )


    ########################################################
    # Store summary
    ########################################################

    summary_list[[summary_index]] <-
        df_summary_chr


    summary_index <-
        summary_index + 1L
}


############################################################
# Combine summary
############################################################

if (length(summary_list) == 0) {

    stop(
        "No chromosome-level results were generated."
    )
}


df_summary <- rbindlist(

    summary_list,

    use.names = TRUE,

    fill = TRUE
)


############################################################
# Restore original lead SNP order
############################################################

df_order <- df_gwas1_snplist[
    ,
    .(
        SNP_gwas1 = leadSNP,
        input_order
    )
]


df_summary <- merge(

    df_order,

    df_summary,

    by = "SNP_gwas1",

    all.x = TRUE,

    sort = FALSE
)


setorder(
    df_summary,
    input_order
)


df_summary[
    ,
    input_order := NULL
]


############################################################
# Summary column order
############################################################

setcolorder(

    df_summary,

    c(
        "SNP_gwas1",

        "In_LD",

        "need_check",

        "n_candidate",

        "n_r2_computed",

        "n_r2_not_computed",

        "SNP_gwas2_in_LD",

        "SNP_gwas2_in_LD_pval",

        "SNP_gwas2_in_LD_R2",

        "SNP_gwas2_not_computed",

        "SNP_gwas2_not_computed_pval"
    )
)


############################################################
# Combine pairwise detail
############################################################

if (length(detail_list) > 0) {

    df_detail <- rbindlist(

        detail_list,

        use.names = TRUE,

        fill = TRUE
    )

} else {

    ########################################################
    # Empty detail table with expected columns
    ########################################################

    df_detail <- data.table(

        SNP_gwas1 = character(),

        POS_gwas1 = numeric(),

        SNP_gwas2 = character(),

        POS_gwas2 = numeric(),

        P_gwas2 = numeric(),

        distance_bp = numeric(),

        lead_in_reference = logical(),

        gwas2_in_reference = logical(),

        R2 = numeric(),

        r2_computed = logical(),

        In_LD_pair = logical(),

        r2_status = character()
    )
}


############################################################
# Output paths
############################################################

output_summary <- file.path(

    dir_out,

    paste0(
        pref_out,
        ".compare_",
        trait1,
        ".to_",
        trait2,
        ".",
        exe_mapper[[delim_out]]
    )
)


output_detail <- file.path(

    dir_out,

    paste0(
        pref_out,
        ".compare_",
        trait1,
        ".to_",
        trait2,
        ".pairwise_detail.",
        exe_mapper[[delim_out]]
    )
)


############################################################
# Save summary
############################################################

fwrite(

    df_summary,

    file = output_summary,

    sep = delim_mapper[[delim_out]],

    quote = FALSE,

    na = "NA"
)


############################################################
# Save pairwise detail
############################################################

fwrite(

    df_detail,

    file = output_detail,

    sep = delim_mapper[[delim_out]],

    quote = FALSE,

    na = "NA"
)


############################################################
# Final report
############################################################

message("")
message(
    "========================================"
)

message(
    "Analysis completed."
)

message("")
message(
    "Summary:"
)

message(
    output_summary
)

message("")
message(
    "Pairwise detail:"
)

message(
    output_detail
)

message("")
message(
    "Total lead SNPs: ",
    nrow(df_summary)
)

message(
    "In_LD = TRUE: ",
    sum(
        df_summary$In_LD,
        na.rm = TRUE
    )
)

message(
    "need_check = yes: ",
    sum(
        df_summary$need_check == "yes",
        na.rm = TRUE
    )
)

message(
    "Total candidate pairs: ",
    nrow(df_detail)
)

message(
    "R2 not computed pairs: ",
    sum(
        !df_detail$r2_computed,
        na.rm = TRUE
    )
)

message(
    "========================================"
)
