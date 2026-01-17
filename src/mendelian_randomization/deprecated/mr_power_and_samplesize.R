compute_power_binary <- function(N, K, OR, R2, alpha=0.05){
    # https://sb452.shinyapps.io/power/
    # Supplementary data page 8 (Burgess, S. Sample size and power calculations in Mendelian randomization with a single instrumental variable and a binary outcome. Int J Epidemiol 43, 922–929 (2014).)
    
    # N: total sample size of outcome (ncase + ncontrol)
    # K: sample prevalence in outcome
    # alpha: Type I error rate
    # OR: True odds ratio for exposure -> outcome
    # R2: Proportion of variance of exposure explained by all of the IVs
    
    ratio <- 1/(K/(1-K))
    power <- pnorm(sqrt(N*R2*(ratio/(1+ratio))*(1/(1+ratio)))*log(OR)-qnorm(1-alpha/2))
    return(power)
}

compute_samplesize_binary <- function(power, K, OR, R2, alpha=0.05){
    # https://sb452.shinyapps.io/power/
    # Supplementary data page 8 (Burgess, S. Sample size and power calculations in Mendelian randomization with a single instrumental variable and a binary outcome. Int J Epidemiol 43, 922–929 (2014).)
    
    # power: power
    # K: sample prevalence in outcome
    # OR: True odds ratio for exposure -> outcome
    # alpha: Type I error rate
    # R2: Proportion of variance of exposure explained by all of the IVs

    ratio <- 1/(K/(1-K))
    N <- (qnorm(1-alpha/2)+qnorm(power))^2/log(OR)^2/R2/(ratio/(1+ratio))/(1/(1+ratio))
    return(N)
}

compute_power_continuous <- function(N, K, b, R2, alpha=0.05){
    # https://sb452.shinyapps.io/power/
    # Supplementary data page 8 (Burgess, S. Sample size and power calculations in Mendelian randomization with a single instrumental variable and a binary outcome. Int J Epidemiol 43, 922–929 (2014).)
    
    # N: total sample size of outcome (ncase + ncontrol)
    # K: sample prevalence in outcome
    # b: True causal effect for exposure -> outcome
    # R2: Proportion of variance of exposure explained by all of the IVs
    # alpha: Type I error rate
    
    ratio <- 1/(K/(1-K))
    power <- pnorm(sqrt(N*R2*(ratio/(1+ratio))*(1/(1+ratio)))*b-qnorm(1-alpha/2))
    return(power)
}

compute_samplesize_continuous <- function(power, K, b, R2, alpha=0.05){
    # https://sb452.shinyapps.io/power/
    # Supplementary data page 8 (Burgess, S. Sample size and power calculations in Mendelian randomization with a single instrumental variable and a binary outcome. Int J Epidemiol 43, 922–929 (2014).)
    
    # power: power
    # K: sample prevalence in outcome
    # alpha: Type I error rate
    # b: True causal effect for exposure -> outcome
    # R2: Proportion of variance of exposure explained by all of the IVs

    ratio <- 1/(K/(1-K))
    N <- (qnorm(1-alpha/2)+qnorm(power))^2/b^2/R2/(ratio/(1+ratio))/(1/(1+ratio))
    return(N)
}