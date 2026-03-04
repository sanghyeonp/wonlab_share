library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(scales)

set.seed(1)

## -------------------------
## Inputs
## -------------------------
file.gwas <- "/data1/sanghyeon/Projects/CMR_GWAS/SAIGE.Junbean/v15/LVEF/SAIGE.step2.fullGRMforNull.LOCO.LVEF.rsID.tsv"
col_snp <- "SNP"; col_chr <- "CHR"; col_pos <- "POS"; col_p <- "p.value"

# Point 관련.
size.point <- 1.2
alpha.point <- 1

# Annotation 하는 SNP 관련.
size.label <- 6 # SNP annotation text label 크기.
color.label <- "black"
size.annot <- 3 # Annotation 할 SNP 크기

# Manhattan plot color
color1 <- "#BFBFBF"
color2 <- "#D9D9D9"

color.hline.GWS <- "red"
color.hline.suggestive <- "blue"

size.axis_text_x <- 16; size.axis_text_y <- 16
size.axis_title <- 20

# Plot axis title
y_axis_title <- "--log<sub>10</sub>(*P*<sub>LVEF</sub>)"

# Output file prefix
outf_pref <- "figure_1__panel_a"

## -------------------------
## Annotation sets (edit only this)
## - name: category id
## - snps: rsIDs
## - locus_color: stripe color for locus points
## - marker_fill: fill color for the SNP marker
## - marker_shape: ggplot shape (21-25 recommended because they support fill)
## - priority: larger number wins if locus windows overlap
## -------------------------

### 1. Annotate 할 SNP 불러와서 vector로 구성하기.
lead_snp_list <- fread("/data1/sanghyeon/Projects/CMR_GWAS/FUMA/v15/leadSNP_table.all_strain.csv", data.table = FALSE) %>%
    filter(Strain == "LVEF") %>%
    pull(rsID)


unreported_snp_list <- c("rs1400771", "rs12521291", "rs11786896")

# Exclude replicated SNPs from lead SNP
lead_snp_list <- setdiff(lead_snp_list, unreported_snp_list)

### 2. Annotate 할 SNP vector 마다, annotation set 구성하기.
# annot_sets: 순서는 가장 순선적으로 annotate 되어야할 (priority가 높은) 것 먼저.
# locus_color: locus의 color. 
# marker_fill: annotation 할 SNP을 어떤 색상으로 표현할 것인지.
annot_sets <- list(
    unreported = list(
        snps = unreported_snp_list,
        locus_color = "purple",
        marker_fill = "purple",
        marker_shape = 23,
        priority = 1
    ),
    lead = list(
        snps = lead_snp_list,
        locus_color = "#6de722",
        marker_fill = "#6de722",
        marker_shape = 21,
        priority = 2
    )
)

## -------------------------
## Other options
## -------------------------
locus_window_bp <- 250000
keep_p_thresh <- 0.05
sample_frac_nonsig <- 0.1

## -------------------------
## Read GWAS
## -------------------------
df <- fread(file.gwas, data.table = FALSE)

df <- df %>%
    distinct(across(all_of(c(col_snp))), .keep_all = TRUE)

## -------------------------
## Helper: build SNP map (SNP -> type)
## -------------------------
snp_map <- bind_rows(lapply(names(annot_sets), function(nm) {
    tibble(
        snp = annot_sets[[nm]]$snps,
        snp_type = nm
    )})) %>%
    distinct() %>%
    filter(!is.na(snp) & snp != "")

all_snps_to_annotate <- unique(snp_map$snp)

## -------------------------
## Helper: build locus intervals for each set
## -------------------------
get_intervals <- function(df, snps, set_name, window_bp, col_snp, col_chr, col_pos) {
    df %>%
        dplyr::filter(.data[[col_snp]] %in% snps) %>%
        transmute(
            locus_type = set_name,
            chr = .data[[col_chr]],
            start_locus = pmax(.data[[col_pos]] - window_bp, 0L),
            end_locus = .data[[col_pos]] + window_bp
        ) %>%
        distinct()
}

intervals <- bind_rows(lapply(names(annot_sets), function(nm) {
    get_intervals(df, annot_sets[[nm]]$snps, nm, locus_window_bp, col_snp, col_chr, col_pos)
}))

## -------------------------
## Keep dense points inside any locus; reduce elsewhere
## (fast overlap using data.table)
## -------------------------
df_dt <- as.data.table(df)
intervals_dt <- as.data.table(intervals)
setkey(intervals_dt, chr, start_locus, end_locus)

# Create point-interval representation for variants
points_dt <- df_dt[, .(
    chr = get(col_chr),
    start = get(col_pos),
    end = get(col_pos),
    SNP = get(col_snp),
    P = get(col_p)
)]

# Overlap join: variants within any locus interval
in_locus_dt <- foverlaps(points_dt, intervals_dt, 
                         by.x = c("chr", "start", "end"),
                         by.y = c("chr", "start_locus", "end_locus"),
                         type = "within", nomatch = 0L)

# Convert to data frame, keep only original columns
df_in_locus <- as.data.frame(in_locus_dt) %>%
    transmute(
        !!col_snp := SNP,
        !!col_chr := chr,
        !!col_pos := start,
        !!col_p   := P
    ) %>%
    distinct()

df_outside <- anti_join(df, df_in_locus, by = c(col_snp, col_chr, col_pos))

sig_outside <- df_outside %>% filter(.data[[col_p]] < keep_p_thresh)

nonsig_outside <- df_outside %>%
    filter(.data[[col_p]] >= keep_p_thresh) %>%
    group_by(.data[[col_chr]]) %>%
    sample_frac(sample_frac_nonsig)

gwas_data <- bind_rows(df_in_locus, sig_outside, nonsig_outside) %>%
    distinct()

## -------------------------
## Add cumulative position (bp_cum)
## -------------------------
df.cum <- gwas_data %>%
    group_by(.data[[col_chr]]) %>%
    summarise(chr_len = max(.data[[col_pos]]), .groups = "drop") %>%
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(.data[[col_chr]], tot)

gwas_data <- gwas_data %>%
    left_join(df.cum, by = col_chr) %>%
    mutate(bp_cum = .data[[col_pos]] + tot)

axis_set <- gwas_data %>%
    group_by(.data[[col_chr]]) %>%
    summarise(center = mean(bp_cum), .groups = "drop")

ylim <- -log10(min(gwas_data[[col_p]])) * 1.1
xlim <- range(gwas_data$bp_cum)

## -------------------------
## Assign locus_type per variant (with priority)
## -------------------------
# Build a per-point locus hit table with priority resolution
hits_dt <- foverlaps(
    as.data.table(gwas_data)[, .(chr = get(col_chr), start = get(col_pos), end = get(col_pos),
                                 bp_cum = bp_cum, SNP = get(col_snp), P = get(col_p), tot = tot)],
    intervals_dt,
    by.x = c("chr", "start", "end"),
    by.y = c("chr", "start_locus", "end_locus"),
    type = "within",
    nomatch = NA
)

# attach priority, pick best locus_type if multiple overlaps
priority_tbl <- tibble(
    locus_type = names(annot_sets),
    priority = sapply(names(annot_sets), function(nm) annot_sets[[nm]]$priority)
)

hits_df <- hits_dt %>%
    as.data.frame() %>%
    transmute(
        !!col_snp := SNP,
        !!col_chr := chr,
        !!col_pos := start,
        !!col_p   := P,
        bp_cum = bp_cum,
        tot = tot,
        locus_type = locus_type
    ) %>%
    left_join(priority_tbl, by = "locus_type") %>%
    arrange(desc(priority)) %>%
    distinct(across(all_of(c(col_snp, col_chr, col_pos, col_p, "bp_cum", "tot"))),
             .keep_all = TRUE) %>%
    select(-priority)

gwas_data <- gwas_data %>%
    left_join(hits_df, by = c(col_snp, col_chr, col_pos, col_p, "bp_cum", "tot"))

## SNP marker info
gwas_data <- gwas_data %>%
    left_join(snp_map, by = setNames("snp", col_snp)) %>%
    mutate(snp_type = ifelse(is.na(snp_type), "other", snp_type))

## -------------------------
## Base Manhattan plot
## -------------------------
p_base <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(.data[[col_p]]), color = factor(.data[[col_chr]]))) +
    scale_x_continuous(
        labels = axis_set[[col_chr]],
        breaks = axis_set$center,
        expand = expansion(mult = 0.005, add = 0),
        limits = xlim
    ) +
    scale_y_continuous(
        expand = expansion(mult = c(0, 0.1)),
        limits = c(0, ylim),
        breaks = pretty_breaks()
    ) +
    scale_color_manual(values = rep(c(color1, color2), length.out = nrow(axis_set))) +
    labs(x = "Chromosome", y = y_axis_title) +
    theme_bw(base_family = "Helvetica") +
    theme(
        axis.title.y = ggtext::element_markdown(size = size.axis_title, face = "bold", color = "black"),
        axis.title.x = element_text(size = size.axis_title, face = "bold", color = "black"),
        axis.text.x  = element_text(size = size.axis_text_x, color = "black"),
        axis.text.y  = element_text(size = size.axis_text_y, color = "black"),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt")
    )

## -------------------------
## Add locus stripes (overlay points) for each locus_type
## -------------------------
add_locus_overlays <- function(p, gwas_data, annot_sets, col_p) {
    
    # Draw in reverse order so first annot_set appears on top
    for (nm in rev(names(annot_sets))) {
        p <- p + geom_point(
            data = gwas_data %>% filter(locus_type == nm),
            aes(x = bp_cum, y = -log10(.data[[col_p]])),
            inherit.aes = FALSE,
            color = annot_sets[[nm]]$locus_color,
            size = size.point,
            alpha = 1
        )
    }
    
    p
}

## -------------------------
## Background panel
## -------------------------
p_bg <- p_base +
    geom_point(size = size.point, alpha = alpha.point)

p_bg <- add_locus_overlays(p_bg, gwas_data, annot_sets, col_p) +
    theme(
        axis.title.y = ggtext::element_markdown(size = size.axis_title, face = "bold", color = "white"),
        axis.title.x = element_text(size = size.axis_title, face = "bold", color = "white"),
        axis.text.x  = element_text(size = size.axis_text_x, color = "white"),
        axis.text.y  = element_text(size = size.axis_text_y, color = "white"),
        axis.line = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.1, color = "black"),
        axis.ticks.y = element_blank()
    )

ggsave(
    paste0(outf_pref, ".background.png"),
    p_bg,
    device = "png",
    type = "cairo",
    dpi = 600,
    width = 14, height = 8, units = "in",
    bg = "transparent"
)

## -------------------------
## Annotation panel: add loci + SNP markers + labels
## -------------------------
annot_points <- gwas_data %>%
    filter(.data[[col_snp]] %in% all_snps_to_annotate) %>%
    distinct()


# Marker aesthetics per type (derived from annot_sets)
shape_map <- setNames(sapply(names(annot_sets), function(nm) annot_sets[[nm]]$marker_shape), names(annot_sets))
fill_map  <- setNames(sapply(names(annot_sets), function(nm) annot_sets[[nm]]$marker_fill),  names(annot_sets))

# p_annot <- p_base +
    # geom_point(size = size.point, alpha = alpha.point)

p_annot <- p_base +
    geom_point(
        data = annot_points,
        aes(shape = snp_type, fill = snp_type),
        color = "black",
        size = size.annot,
        stroke = 1
    ) +
    ggrepel::geom_text_repel(
        data = annot_points,
        aes(label = .data[[col_snp]]),
        size = size.label,
        color = color.label,
        nudge_y = 1,
        max.overlaps = Inf,
        min.segment.length = 0,
        segment.color = "black",
        force = 1,
        show.legend = FALSE
    ) +
    scale_shape_manual(values = shape_map) +
    scale_fill_manual(values = fill_map) +
    geom_hline(yintercept = -log10(5e-8), color = color.hline.GWS, linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = -log10(1e-5), color = color.hline.suggestive, linetype = "dashed", alpha = 0.7)

ggsave(
    paste0(outf_pref, ".annot.pdf"),
    p_annot,
    device = cairo_pdf,
    width = 14, height = 8, units = "in",
    bg = "transparent"
)


