

make_annovar_input <- function(
    gwas, delim_in,
    chr_col, pos_col, ref_col, alt_col,
    outf,
    nthread
){
    df <- as.data.frame(fread(gwas, 
                            sep=delim_in, 
                            header = TRUE,
                            nThread = nthread,
                            showProgress = FALSE))
    
    pos_col2 <- paste0(pos_col, "_2")

    df_temp <- df %>%
        dplyr::select(!!as.name(chr_col), !!as.name(pos_col), !!as.name(ref_col), !!as.name(alt_col)) %>%
        mutate(!!as.name(pos_col2) := !!as.name(pos_col)) %>%
        dplyr::select(!!as.name(chr_col), !!as.name(pos_col), !!as.name(pos_col2), !!as.name(ref_col), !!as.name(alt_col))

    filename_annovin <- paste0(outf, ".annovin")
    write.table(df_temp,
                filename_annovin,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE)
    
    filename_annovin_flip <- paste0(outf, ".flip.annovin")
    write.table(dplyr::select(df_temp, !!as.name(chr_col), !!as.name(pos_col2), !!as.name(pos_col), !!as.name(alt_col), !!as.name(ref_col)),
                filename_annovin_flip,
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE)

    return (list(filename_annovin, filename_annovin_flip, df))
}
