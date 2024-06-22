library(dplyr)

reg_result_to_df <- function(reg.obj, rename_col=T, odds_ratio=T){
    df <- as.data.frame(summary(reg.obj)$coefficients)
    variable_name <- rownames(df)
    df$Variable <- variable_name; rownames(df) <- NULL
    df <- df[,c(5, 1, 2, 3, 4)]
    if (odds_ratio){
        df <- df %>%
            mutate(Estimate.lo = Estimate - (1.96*`Std. Error`),
                   Estimate.hi = Estimate + (1.96*`Std. Error`),
                   OR = exp(Estimate), OR.lo = exp(Estimate.lo), OR.hi = exp(Estimate.hi))
    }
    if (rename_col){
        rename_lookup <- c("SE"="Std. Error", "Z"="z value", "P"="Pr(>|z|)")
        df <- df %>%
            rename(any_of(rename_lookup))
    }
    return(df)
}
