compute_power_binary <- function(N, K, OR, R2, alpha=0.05){
    # N: total sample size of outcome (ncase + ncontrol)
    # alpha: Type I error rate
    # K: sample prevalence in outcome
    # OR: True odds ratio for exposure -> outcome
    # R2: Proportion of variance of exposure explained by all of the IVs
    
    chisquare_threshold <- qchisq(1-alpha, 1)
    
    b_mr <- K * ((OR/(1 + (K*(OR-1)))) - 1)
    v_mr <- ((K*(1-K)) - b_mr^2)/(N*R2)
    NCP <- b_mr^2 / v_mr
    
    power <- 1 - pchisq(chisquare_threshold, 1, NCP)

    return(power)  

    # threschi <- qchisq(1 - alpha, 1) # threshold chi(1) scale
    # # f.value <- 1 + N * R2xz / (1 - R2xz)
    # b_MR <- K * ( OR/ (1 + K * (OR - 1)) -1)
    # v_MR <- (K * (1-K) - b_MR^2) / (N*R2)
    # NCP <- b_MR^2 / v_MR
    
    # # 2-sided test
    # power <- 1 - pchisq(threschi, 1, NCP) 

    # return(power) 
}

compute_samplesize_binary <- function(K, OR, R2, power, alpha=0.05){
    # Calculation of sample size given power
    z1 <- qnorm(1 - alpha / 2)
    z2 <- qnorm(power)
    Z  <- (z1 + z2)^2

    b_01 <- K * ( OR/(1 + K * (OR - 1)) - 1)
    f <- K * (1-K) - b_01^2
    N1 <- Z * f / (b_01^2 * R2)
    N1 <- ceiling(N1)
    
    return(N1)
}