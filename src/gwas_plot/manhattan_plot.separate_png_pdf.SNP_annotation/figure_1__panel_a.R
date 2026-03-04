library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(scales)
set.seed(1)


file.gwas <- "/data1/sanghyeon/Projects/CMR_GWAS/SAIGE.Junbean/v15/LVEF/SAIGE.step2.fullGRMforNull.LOCO.LVEF.rsID.tsv"
df <- fread(file.gwas, data.table=F, nThread=10)

# df.copy <- df

lead_snp <- fread("/data1/sanghyeon/Projects/CMR_GWAS/FUMA/v15/leadSNP_table.all_strain.csv", data.table=F) %>%
    filter(Strain == "LVEF") %>%
    pull(rsID)

replicated_snp <- c("rs1400771", "rs12521291", "rs11786896")

col_snp<-"SNP"; col_chr <- "CHR"; col_pos <- "POS"; col_p <- "p.value"

## Defaults
annot.snplist <- lead_snp; annot_scale <- 1.5
color.annot <- "#F4A261"
annot.label <- T
size.label <- 6; color.label <- "black"

size.point <- 2; alpha.point <- 1
# color1 <- "#2A9D8F"; color2 <- "#1E6460"
color1 <- "#BB5566"; color2 <- "#883E47"
hline.GWS <- T; color.hline.GWS <- "red"
hline.suggestive <- T; color.hline.suggestive <- "blue"
size.axis_text_x <- 16; size.axis_text_y <- 16
size.axis_title <- 20

####
## Reduce computation by reducing non-significant pints
sig_data <- df %>% filter(!!as.name(col_p) < 0.05)

notsig_data <- df %>%
    filter(!!as.name(col_p) >= 0.05) %>%
    group_by(!!as.name(col_chr)) %>%
    sample_frac(0.1)

gwas_data <- bind_rows(sig_data, notsig_data)

## 1. column-wise data
# df.cum <- gwas_data %>%
#     # Compute chromosome size
#     group_by(!!as.name(col_chr)) %>%
#     summarise(max_bp = max(!!as.name(col_pos))) %>%
#     # Calculate cumulative position of each chromosome
#     mutate(bp_add = cumsum(as.numeric(max_bp)), default=0) %>%
#     dplyr::select(!!as.name(col_chr), bp_add)

df.cum <- gwas_data %>%
    # Compute chromosome size
    group_by(!!as.name(col_chr)) %>%
    summarise(chr_len = max(!!as.name(col_pos))) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
    dplyr::select(-chr_len)

gwas_data <- gwas_data %>%
    inner_join(df.cum, by=col_chr) %>%
    mutate(bp_cum = !!as.name(col_pos) + tot)

axis_set <- gwas_data %>%
    group_by(!!as.name(col_chr)) %>%
    summarize(center = mean(bp_cum))

ylim <- -log10(min(gwas_data[[col_p]]))*1.1
xlim <- c(min(gwas_data$bp_cum), max(gwas_data$bp_cum))

################
gwas_data$snp_type <- ifelse(
    gwas_data[[col_snp]] %in% replicated_snp,
    "replicated",
    ifelse(gwas_data[[col_snp]] %in% annot.snplist, "annot", "other")
)


#################
p_base <- ggplot(gwas_data, aes(x = bp_cum, y = -log10(.data[[col_p]]),
                                color = factor(.data[[col_chr]]))) +
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
    labs(x = "Chromosome", y = "--log<sub>10</sub>(*P*<sub>LVEF</sub>)") +
    theme_bw(base_family = "Helvetica") +  # use a font that exists on Linux
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
        plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt")  # lock margins explicitly
    )


p_bg <- p_base + geom_point(size = size.point, alpha = alpha.point) +
    theme(
        axis.title.y = ggtext::element_markdown(size = size.axis_title, face = "bold", color = "white"),
        axis.title.x = element_text(size = size.axis_title, face = "bold", color = "white"),
        axis.text.x  = element_text(size = size.axis_text_x, color = "white"),
        axis.text.y  = element_text(size = size.axis_text_y, color = "white"),
        axis.line = element_blank(), 
        axis.ticks.x = element_line(linewidth=0.1, color="black"),
        axis.ticks.y = element_blank()
    )

p_annot <- p_base +
    geom_point(
        data = gwas_data %>% filter(.data[[col_snp]] %in% annot.snplist | .data[[col_snp]] %in% replicated_snp),
        aes(shape = snp_type, fill = snp_type),
        color = "black",
        size = size.point * annot_scale, stroke = 1.2
    ) +
    ggrepel::geom_text_repel(
        data = gwas_data %>% filter(.data[[col_snp]] %in% annot.snplist | .data[[col_snp]] %in% replicated_snp),
        aes(label = .data[[col_snp]]),
        size = size.label, color = color.label,
        nudge_y = 1,
        segment.curvature = -0.1,
        min.segment.length = 0,
        segment.color = "black",
        force = 1,
        show.legend = FALSE
    ) +
    scale_shape_manual(values = c("annot" = 21, "replicated" = 23)) +
    scale_fill_manual(values = c("annot" = color.annot, "replicated" = "purple")) +
    geom_hline(yintercept = -log10(5e-8), color = color.hline.GWS, linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = -log10(1e-5), color = color.hline.suggestive, linetype = "dashed", alpha = 0.7)

# 2) Save PNG using cairo to match cairo_pdf more closely
ggsave(
    "figure_1__panel_a.background.png",
    p_bg,
    device = "png",
    type = "cairo",
    dpi = 600,
    width = 14, height = 8, units = "in",
    bg = "transparent"
)

# 3) Save PDF with cairo_pdf
ggsave(
    "figure_1__panel_a.annot.pdf",
    p_annot,
    device = cairo_pdf,
    width = 14, height = 8, units = "in",
    bg = "transparent"
)


