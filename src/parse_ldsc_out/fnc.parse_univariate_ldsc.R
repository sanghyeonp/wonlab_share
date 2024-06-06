# Sanghyeon Park
# 2024.06.05

library(stringr)

parse_univariate_ldsc <- function(log_file, trait=NA){
    log_content <- readLines(log_file)

    observed_h2 <- NA; observed_h2_se <- NA; liability_h2 <- NA; liability_h2_se <- NA
    lambda_gc <- NA; mean_chi_sq <- NA; intercept <- NA; intercept_se <- NA
    ratio <- NA; ratio_se <- NA

    for (line in log_content) {
        if (grepl("Total Observed scale h2:", line)) {
            # Extract Genetic Correlation
            query1 <- gsub("Total Observed scale h2: ", "", line)
            observed_h2 <- str_split(query1, " ")[[1]][1]
            observed_h2_se <- gsub("\\)", "", gsub("\\(", "", paste(str_split(query1, " ")[[1]][-1], collapse=" ")))
        }
        if (grepl("Total Liability scale h2:", line)) {
            # Extract Genetic Correlation
            query1 <- gsub("Total Liability scale h2: ", "", line)
            liability_h2 <- str_split(query1, " ")[[1]][1]
            liability_h2_se <- gsub("\\)", "", gsub("\\(", "", paste(str_split(query1, " ")[[1]][-1], collapse=" ")))
        }
        if (grepl("Lambda GC:", line)) {
            lambda_gc <- gsub("Lambda GC: ", "", line)
        }
        if (grepl("Mean Chi^2:", line)) {
            mean_chi_sq <- gsub("Mean Chi^2: ", "", line)
        }
        if (grepl("Intercept:", line)) {
            query1 <- gsub("Intercept: ", "", line)
            intercept <- str_split(query1, " ")[[1]][1]
            intercept_se <- gsub("\\)", "", gsub("\\(", "", paste(str_split(query1, " ")[[1]][-1], collapse=" ")))
        }
        if (grepl("Ratio", line)) {
            if (grepl("Ratio <", line)){
                ratio <- gsub("\\.", "", gsub("Ratio ", "", line))
            } else{
                query1 <- gsub("Ratio: ", "", line)
                ratio <- str_split(query1, " ")[[1]][1]
                ratio_se <- gsub("\\)", "", gsub("\\(", "", paste(str_split(query1, " ")[[1]][-1], collapse=" ")))
            }
        }
    }
    if (is.na(trait)) {
        return(data.frame(`Observed-scale h2`=observed_h2, `Observed-scale h2 SE`=observed_h2_se,
                        `Liability-scale h2`=liability_h2, `Liability-scale h2 SE`=liability_h2_se,
                        `Lambda GC`=lambda_gc, `Mean Chi^2`=mean_chi_sq, 
                        `Intercept`=intercept, `Intercept SE`=intercept_se,
                        `Ratio`=ratio, `Ratio SE`=ratio_se))
    } 
    return(data.frame(trait=trait, `Observed-scale h2`=observed_h2, `Observed-scale h2 SE`=observed_h2_se,
                `Liability-scale h2`=liability_h2, `Liability-scale h2 SE`=liability_h2_se,
                `Lambda GC`=lambda_gc, `Mean Chi^2`=mean_chi_sq, 
                `Intercept`=intercept, `Intercept SE`=intercept_se,
                `Ratio`=ratio, `Ratio SE`=ratio_se))
}