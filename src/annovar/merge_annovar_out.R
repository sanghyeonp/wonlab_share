
merge_annovar_out <- function(
    filename_annovin_annot,
    filename_annovin_flip_annot,
    nthread
){
    # Result without flipping
    df <- as.data.frame(fread(filename_annovin_annot, 
                            header = FALSE,
                            sep = "\t",
                            nThread = nthread,
                            showProgress = FALSE))

    # Get column values to query
    idx_col_list <- c(3, 4, 5)

    df <- mutate(df, 
                flipped = FALSE)  
    df$temp_var <- apply(df[, idx_col_list], 1, function(x) paste(x, collapse = ":"))

    unflipped_index_value_list <- df$temp_var

    # Result with flipping
    df_flip <- as.data.frame(fread(filename_annovin_flip_annot, 
                                    header = FALSE,
                                    sep = "\t",
                                    nThread = nthread,
                                    showProgress = FALSE))
    df_flip <- mutate(df_flip, 
                        flipped = TRUE)
    df_flip$temp_var <- apply(df_flip[, idx_col_list], 1, function(x) paste(x, collapse = ":"))

    flipped_index_value_list <- df_flip$temp_var

    ###
    cat(paste0("\nMERGE: Merging unflipped and flipped ANNOVAR output files...\n\t", 
        "Unflipped ANNOVAR output: ", filename_annovin_annot, " (Rows: ", format(nrow(df), big.mark = ",", scientific = FALSE), ")", "\n\t", 
        "Flipped ANNOVAR ouput: ", filename_annovin_flip_annot, " (Rows: ", format(nrow(df_flip), big.mark = ",", scientific = FALSE), ")", "\n"))
    ###

    # Query rows to keep from flipped result
    rows_to_keep <- !(flipped_index_value_list %in% unflipped_index_value_list)

    df_flip_to_keep <- df_flip[rows_to_keep, ]

    # Merge
    df_merged <- rbind(df, df_flip_to_keep)

    df_merged <- df_merged %>% 
        select(-temp_var)

    if (args$rm_flip_col){
        df_merged <- df_merged %>% 
            select(-flipped)
    }

    # Write
    file_path <- paste0(filename_annovin_annot, ".merged")

    ###
    cat(paste0("\tMerged ANNOVAR output: ", file_path, " (Rows: ", format(nrow(df_merged), big.mark = ",", scientific = FALSE), ")", "\n"))
    ###

    write.table(df_merged, 
                file = file_path, 
                quote = FALSE, 
                sep = "\t", 
                row.names = FALSE, 
                col.names = FALSE)

    return (file_path)
}

