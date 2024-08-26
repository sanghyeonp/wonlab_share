library(stringr)

parse_popcorn <- function(log_file, trait1=NA, trait2=NA){
    log_content <- readLines(log_file)

    keep_count <- 0
    df.out <- data.frame()
    for (line in log_content) {
        if (keep_count > 0){
            line1 <- strsplit(line, " ")[[1]]; line2 <- line1[which(line1 != "")]; line2 <- gsub("\n", "", line2)
            df.out <- rbind(df.out,
                            data.frame("Variable" = line2[1],
                                    "Val(obs)" = line2[2],
                                    "SE" = line2[3],
                                    "Z" = line2[4],
                                    "P(Z)" = line2[5], check.names=FALSE))
            keep_count <- keep_count - 1
        }
        if (grepl("*P (Z)*", line)){
            keep_count <- 3
        }
    }
    return(df.out)
}
