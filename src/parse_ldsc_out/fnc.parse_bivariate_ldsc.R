# Sanghyeon Park
# 2024.06.05

library(stringr)

parse_bivariate_ldsc <- function(log_file, trait1=NA, trait2=NA){
    log_content <- readLines(log_file)

    genetic_correlation <- NA; se <- NA; z_score <- NA; p_value <- NA
    for (line in log_content) {
        if (grepl("Genetic Correlation:", line)) {
            # Extract Genetic Correlation
            query1 <- gsub("Genetic Correlation: ", "", line)
            genetic_correlation <- str_split(query1, " ")[[1]][1]
            se <- gsub("\\)", "", gsub("\\(", "", paste(str_split(query1, " ")[[1]][-1], collapse=" ")))
        }
        if (grepl("Z-score:", line)) {
            # Extract Z-score
            z_score <- gsub("Z-score: ", "", line)
        }
        if (grepl("P:", line)) {
            # Extract P-value
            p_value <- gsub("P: ", "", line)
        }
    }
    if (is.na(trait1) | is.na(trait2)) {
        return(data.frame(rg=genetic_correlation, se=se, Z=z_score, P=p_value))
    } 
    return(data.frame(trait1=trait1, trait2=trait2, rg=genetic_correlation, se=se, Z=z_score, P=p_value))
}
