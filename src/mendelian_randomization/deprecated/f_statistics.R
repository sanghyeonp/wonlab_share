compute_R2 <- function(b, se, maf, N){
    # b: effect size of the genetic variant on the exposure
    # se: standard error of the effect size
    # maf: minor allele frequency of the genetic variant
    # N: total sample size of exposure
    if (maf > 0.5){
        maf <- 1-maf
    }
    R2 <- (2*b^2*maf*(1-maf)) / ((2*b^2*maf*(1-maf)) + (se^2*2*N*maf*(1-maf)))
    return(R2)
}

compute_f_stat <- function(R2, N, k=1){
    # R2: Proportion of variance of exposure explained by all of the IVs
    # N: total sample size of exposure
    # k: number of IVs (in most of cases, we compute F for a single IV)
    f_stat <- (R2*(N-k-1))/(k*(1-R2))
    return(f_stat)
}