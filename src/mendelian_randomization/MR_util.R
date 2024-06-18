compute_R2 <- function(b, SE, MAF, N){
    # b: effect size of the genetic variant on the exposure
    # SE: standard error of the effect size
    # MAF: minor allele frequency of the genetic variant
    # N: total sample size of exposure
    if (MAF > 0.5) MAF <- 1-MAF
    
    R2 <- (2*b^2*MAF*(1-MAF)) / ((2*b^2*MAF*(1-MAF)) + (SE^2*2*N*MAF*(1-MAF)))
    return(R2)
}

compute_f_stat <- function(R2, N, k=1){
    # R2: Proportion of variance of exposure explained by all of the IVs
    # N: total sample size of exposure
    # k: number of IVs (in most of cases, we compute F for a single IV)
    f_stat <- (R2*(N-k-1))/(k*(1-R2))
    return(f_stat)
}

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
