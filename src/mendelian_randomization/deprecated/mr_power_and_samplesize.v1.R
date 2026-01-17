"Version 1
Brion, M.-J. A., Shakhbazov, K. & Visscher, P. M. Calculating statistical power in Mendelian randomization studies. International Journal of Epidemiology 42, 1497–1501 (2013).
https://github.com/kn3in/mRnd/blob/master/functions.R
위 논문과 mRnd web tool을 참고하여 작성하였지만, 아래 여러 단점들로 Burgess 2014 논문의 방법을 사용하도록 변경할 예정.
    - Replication이 어렵다는 단점
    - Burgess 2014에 비해서 citation이 적음 (약 1/3)
    - Burgess 2014보다 조금 더 inflated power를 계산해주는 것 같음
"


compute_power_binary <- function(N, K, OR, R2, alpha=0.05){
    # https://github.com/kn3in/mRnd/blob/master/functions.R
    
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
    # https://github.com/kn3in/mRnd/blob/master/functions.R
    
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