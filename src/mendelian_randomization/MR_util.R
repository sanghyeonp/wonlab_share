
compute_R2 <- function(b, maf){
    ## Reference: Burgess, S., Dudbridge, F. & Thompson, S.G. Combining information on multiple instrumental variables in Mendelian randomization: comparison of allele score and summarized data methods. Statistics in Medicine 35, 1880-1906 (2016).
    # b: marginal effect
    # maf: minor allele frequency
    b <- as.numeric(b); maf <- as.numeric(maf)
    if(maf > 0.5){maf <- 1-maf}
    r2 <- 2*b^2 * maf * (1-maf)
    return(r2)
}

compute_F <- function(r2, n, k){
    ## Reference: Burgess, S., Dudbridge, F. & Thompson, S.G. Combining information on multiple instrumental variables in Mendelian randomization: comparison of allele score and summarized data methods. Statistics in Medicine 35, 1880-1906 (2016).
    # r2: Proportion of variance of exposure explained by the SNP
    # n: Sample size of exposure
    # k: Number of IVs
    r2 <- as.numeric(r2); n <- as.numeric(n); k <- as.numeric(k)
    F <- ((n-k-1)/k) * (r2/(1-r2))
    return(F)
}

# compute_f_stat <- function(b, se){
#     b <- as.numeric(b); se <- as.numeric(se)
#     return (b^2 / se^2)
# }



########################################################################################
### Power calculation
# ShinyApp (https://cnsgenomics.shinyapps.io/mRnd/)의 function을 그대로 활용
# https://github.com/kn3in/mRnd/blob/master/functions.R
# Citation: Brion, Marie-Jo A., Konstantin Shakhbazov, and Peter M. Visscher. "Calculating statistical power in Mendelian randomization studies." International journal of epidemiology 42.5 (2013): 1497-1501.
########################################################################################

mRnd.binary_outcome <- function(N, R2xz, K, OR, alpha=0.05, epower=NA) {
    ### Parameters
    # N: Total sample size of outcome
    # R2xz: Proportion of variance of exposure explained by the IV (use compute_R2())
    # K: Proportion of cases in the outcome
    # OR: Hypothetical true OR
    # alpha: Type I error rate
    # epower: Hypothetical power
    N <- as.numeric(N); R2xz <- as.numeric(R2xz); K <- as.numeric(K); OR <- as.numeric(OR); alpha <- as.numeric(alpha); epower <- as.numeric(epower)

    threschi <- qchisq(1 - alpha, 1) # threshold chi(1) scale
    # f.value <- 1 + N * R2xz / (1 - R2xz)

    if (is.na(epower)) {
        b_MR <- K * ( OR/ (1 + K * (OR - 1)) -1)
        
        v_MR <- (K * (1-K) - b_MR^2) / (N*R2xz)
        NCP <- b_MR^2 / v_MR
        
        # 2-sided test
        power <- 1 - pchisq(threschi, 1, NCP)
        return(power)
    } else {
        # Calculation of sample size given power
        z1 <- qnorm(1 - alpha / 2)
        z2 <- qnorm(epower)
        Z  <- (z1 + z2)^2

        b_01 <- K * ( OR/ (1 + K * (OR - 1)) -1)
        f <- K * (1-K) - b_01^2
        n1 <- Z * f / (b_01^2 * R2xz)
        n1 <- ceiling(n1)
        return(n1)
    }
}

mRnd.continuous <- function(N, byx, bOLS, R2xz, varx, vary, alpha=0.05, epower=NA) {
    ### Parameters
    # N: Total sample size of outcome
    # byx: Hypothetical true causal effect
    # bOLS: 
    # R2xz: Proportion of variance of exposure explained by the IV (use compute_R2())
    # varx: Variance of exposure (if GWAS was done on the standardized values, varx = 1)
    # vary: Variance of outcome (if GWAS was done on the standardized values, vary = 1)
    # alpha: Type I error rate
    # epower: Hypothetical power

    threschi <- qchisq(1 - alpha, 1) # threshold chi(1) scale
    # f.value <- 1 + N * R2xz / (1 - R2xz)
    con <- (bOLS - byx) * varx # covariance due to YX confounding
    vey <- vary - byx * varx * (2 * bOLS - byx)
    
    if (vey < 0) {
        cat("\nError: Invalid input. The provided parameters result in a negative estimate for variance of the error term in the two-stage least squares model.\n")
        return(NA)
    } else {
        if (is.na(epower)) {
            b2sls <- byx + con / (N * R2xz)
            v2sls <- vey / (N * R2xz * varx)
            NCP <- b2sls^2 / v2sls
            # 2-sided test
            power <- 1 - pchisq(threschi, 1, NCP)
            return(power)
        } else {
            # Calculation of sample size given power
            z1 <- qnorm(1 - alpha / 2)
            z2 <- qnorm(epower)
            Z  <- (z1 + z2)^2
            # Solve quadratic equation in N
            a <- (byx * R2xz)^2
            b <- R2xz * (2 * byx * con - Z * vey / varx)
            c <- con^2
            N1 <- ceiling((-b + sqrt(b^2 - 4 * a * c)) / (2 * a))
            return(N1)
        }
    }
}



########################################################################################
### Power curve
# OR = 1 을 기준으로 symmetric 하기 때문에 그냥 OR = 1 부터 보여주면 됨.
########################################################################################

or_list <- seq(1, 2, by=0.01)
power_list <- sapply(or_list, function(or){mRnd.binary_outcome(N=79429, R2xz=0.07775461, K=7495/79429, OR=or)})
df.plot <- data.frame(OR=or_list, Power=power_list)
library(ggplot2)
ggplot(data=df.plot, aes(x=OR, y=Power)) +
    geom_point() +
    geom_hline(yintercept=0.8)














### If R2 is known
# compute_f_stat <- function(R2, N, k=1){
#     # R2: Proportion of variance of exposure explained by all of the IVs
#     # N: total sample size of exposure
#     # k: number of IVs (in most of cases, we compute F for a single IV)
#     f_stat <- (R2*(N-k-1))/(k*(1-R2))
#     return(f_stat)
# }


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
