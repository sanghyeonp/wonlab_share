library(dplyr)
library(ggplot2)
library(ggrepel)

plot_manhattan <- function(gwas, 
                        snp_col, chr_col, pos_col, p_col, 
                        beta_col='NA', se_col='NA',
                        snps_to_annotate=c('NA'), color_annotate='darkorange1',
                        color1='grey50', color2='grey', 
                        chr_select=c('NA'),
                        img_type='png', dpi=300,
                        outf='NA', outd='NA',
                        width=180, height=100, units='mm',
                        scale=1){

    if (outf == "NA"){
        outf <- paste0("manhattan.", basename(gwas))
    }
    if (outd == "NA"){
        outd <- getwd()
    }

    dpi <- as.numeric(dpi)
    width <- as.numeric(width)
    height <- as.numeric(height)

    ### Read input file
    df <- fread(gwas)

    ### Rename columns
    if((beta_col != "NA") & (se_col != "NA")){
        df <- df %>%    
            select(all_of(c(snp_col, chr_col, pos_col, beta_col, se_col))) %>%
            rename("SNP" := snp_col,
                    "CHR" := chr_col,
                    "POS" := pos_col,
                    "BETA" := beta_col,
                    "SE" := se_col
            ) %>%
            mutate(PVAL = 2 * pnorm(abs(BETA / SE), lower.tail=F))
    } else{
        df <- df %>% 
            select(all_of(c(snp_col, chr_col, pos_col, p_col))) %>%
            rename("SNP" := snp_col,
                    "CHR" := chr_col,
                    "POS" := pos_col,
                    "PVAL" := p_col
            )
    }

    ### Filter chromosomes
    if (chr_select[1] != "NA"){
        chr_select <- as.numeric(chr_select)
        df <- df %>% 
                filter(CHR %in% chr_select)
    }

    n_chr <- length(unique(df$CHR))

    #####################
    # # 여기 P-value를 새로 계산하고, given p-value랑 차이가 많이 나는지 확인하는 코드 작성하기.
    # # 차이가 많이 안 난다면, 새로 계산한 P-value 이용. 
    # if(beta_col != "NA" & se_col != "NA"){
    #   df <- df %>% 
    #     rename("BETA" := beta_col,
    #            "SE" := se_col
    #     )
    #   df$PVAL_cal <- 2 * pnorm(-abs(df$BETA / df$SE))
    #   
    # }


    #####################

    ### Select columns
    df <- df %>%
        select(all_of(c("SNP", "CHR", "POS", "PVAL")))

    ### Check rows with NA
    df_na <- df[!complete.cases(df), ]

    if(nrow(df_na) != 0){
        cat(paste0("Number of removed SNPs with missing value: ", nrow(df_na)))

        write.table(df_na, "manhattan_SNPs_with_missing_value.txt",
                    sep="\t", row.names=FALSE, col.names=TRUE
        )
    }

    ### Remove rows with NA
    df <- df[complete.cases(df), ]

    ### Prepare the dataset
    df_plot <- df %>% 
        # Compute chromosome size
        group_by(CHR) %>%
        summarise(chr_len = max(POS)) %>%

        # Calculate cumulative position of each chromosome
        mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
        select(-chr_len) %>%

        # Add this info to the initial dataset
        left_join(df, ., by=c("CHR"="CHR")) %>%

        # Add a cumulative position of each SNP
        arrange(CHR, POS) %>%
        mutate(BPcum = POS + tot) %>%

        # Add highlight and annotation information
        mutate(is_highlight = ifelse(SNP %in% snps_to_annotate, "yes", "no"))


    ### Prepare X-axis
    axisdf <- df_plot %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


    ### Plot
    p1 <- ggplot(df_plot, aes(x=BPcum, y=-log10(PVAL))) +

        # Show all points
        geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
        scale_color_manual(values = rep(c(color1, color2), n_chr)) +
        
        # Add horizontal line (Genome-wide significant)
        geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red", linewidth = 0.8) +
        
        # Add horizontal line (Suggestive)
        geom_hline(yintercept = -log10(1e-5), linetype = "dashed", color = "blue", linewidth = 0.8) +
        
        # custom X axis:
        scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
        # scale_y_continuous(limits = c(0, -log10(min(df_plot$PVAL))), expand = c(0, 0) ) +     # remove space between plot area and x axis
        
        # Add highlighted points
        geom_point(data=subset(df_plot, is_highlight=="yes"), color=color_annotate, size=1.5) +

        # Label name
        labs(x = "Chromosome",
            y = expression(bold("-log"[10](italic(P))))
            ) +
        
        # Custom the theme:
        theme_bw() +
        theme( 
            legend.position="none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.title = element_text(size = 14, face="bold"),
            axis.text = element_text(size = 14)
        )

    if(img_type == "pdf"){
        ggsave(paste0(outf, ".", img_type), p1, 
        scale = scale, device = "pdf",
        width = width, height = height,
        units = units)
    } else{
        ggsave(paste0(outf, ".", img_type), p1, 
        scale = scale, dpi = dpi,
        width = width, height = height,
        units = units)
    }

}


################################################################################################

plot_manhattan_df_input <- function(df, 
                        snp_col, chr_col, pos_col, p_col, 
                        beta_col='NA', se_col='NA',
                        snps_to_annotate=c('NA'), color_annotate='darkorange1',
                        color1='grey50', color2='grey', point_size=1, point_size_annot=1.5,
                        chr_select=c('NA'),
                        img_type='png', dpi=300,
                        outf='NA', outd='NA',
                        width=180, height=100, units='mm',
                        scale=1, save_img=T){

    if (outf == "NA"){
        outf <- "manhattan"
    }
    if (outd == "NA"){
        outd <- getwd()
    }

    dpi <- as.numeric(dpi)
    width <- as.numeric(width)
    height <- as.numeric(height)

    ### Read input file
    # df <- fread(gwas)

    ### Rename columns
    if((beta_col != "NA") & (se_col != "NA")){
        df <- df %>%    
            select(all_of(c(snp_col, chr_col, pos_col, beta_col, se_col))) %>%
            rename("SNP" := snp_col,
                    "CHR" := chr_col,
                    "POS" := pos_col,
                    "BETA" := beta_col,
                    "SE" := se_col
            ) %>%
            mutate(PVAL = 2 * pnorm(-abs(BETA / SE)))
    } else{
        df <- df %>% 
            select(all_of(c(snp_col, chr_col, pos_col, p_col))) %>%
            rename("SNP" := snp_col,
                    "CHR" := chr_col,
                    "POS" := pos_col,
                    "PVAL" := p_col
            )
    }

    ### Filter chromosomes
    if (chr_select[1] != "NA"){
        chr_select <- as.numeric(chr_select)
        df <- df %>% 
                filter(CHR %in% chr_select)
    }

    n_chr <- length(unique(df$CHR))

    #####################
    # # 여기 P-value를 새로 계산하고, given p-value랑 차이가 많이 나는지 확인하는 코드 작성하기.
    # # 차이가 많이 안 난다면, 새로 계산한 P-value 이용. 
    # if(beta_col != "NA" & se_col != "NA"){
    #   df <- df %>% 
    #     rename("BETA" := beta_col,
    #            "SE" := se_col
    #     )
    #   df$PVAL_cal <- 2 * pnorm(-abs(df$BETA / df$SE))
    #   
    # }


    #####################

    ### Select columns
    df <- df %>%
        select(all_of(c("SNP", "CHR", "POS", "PVAL")))

    ### Check rows with NA
    df_na <- df[!complete.cases(df), ]

    if(nrow(df_na) != 0){
        cat(paste0("Number of removed SNPs with missing value: ", nrow(df_na)))

        write.table(df_na, "manhattan_SNPs_with_missing_value.txt",
                    sep="\t", row.names=FALSE, col.names=TRUE
        )
    }

    ### Remove rows with NA
    df <- df[complete.cases(df), ]

    ### Prepare the dataset
    df_plot <- df %>% 
        # Compute chromosome size
        group_by(CHR) %>%
        summarise(chr_len = max(POS)) %>%

        # Calculate cumulative position of each chromosome
        mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
        select(-chr_len) %>%

        # Add this info to the initial dataset
        left_join(df, ., by=c("CHR"="CHR")) %>%

        # Add a cumulative position of each SNP
        arrange(CHR, POS) %>%
        mutate(BPcum = POS + tot) %>%

        # Add highlight and annotation information
        mutate(is_highlight = ifelse(SNP %in% snps_to_annotate, "yes", "no"))


    ### Prepare X-axis
    axisdf <- df_plot %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


    ### Plot
    p1 <- ggplot(df_plot, aes(x=BPcum, y=-log10(PVAL))) +

        # Show all points
        geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=point_size) +
        scale_color_manual(values = rep(c(color1, color2), n_chr)) +
        
        # Add horizontal line (Genome-wide significant)
        geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red", linewidth = 0.8) +
        
        # Add horizontal line (Suggestive)
        geom_hline(yintercept = -log10(1e-5), linetype = "dashed", color = "blue", linewidth = 0.8) +
        
        # custom X axis:
        scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center ) +
        scale_y_continuous(expand=expansion(mult = c(0, 0.1))) +     # remove space between plot area and x axis
        
        # Add highlighted points
        geom_point(data=subset(df_plot, is_highlight=="yes"), color=color_annotate, size=point_size_annot) +

        # Label name
        labs(x = "Chromosome",
            y = expression(bold("-log"[10](italic(P))))
            ) +
        
        # Custom the theme:
        theme_bw() +
        theme( 
            legend.position="none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.title = element_text(size = 14, face="bold"),
            axis.text = element_text(size = 14)
        )

    if (save_img){
        if(img_type == "pdf"){
            ggsave(paste0(outf, ".", img_type), p1, 
            scale = scale, device = "pdf",
            width = width, height = height,
            units = units)
        } else{
            ggsave(paste0(outf, ".", img_type), p1, 
            scale = scale, dpi = dpi,
            width = width, height = height,
            units = units)
        }
    }


    return (p1)
}








