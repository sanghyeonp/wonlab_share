# Compute Beta
## from Odds ratio
Beta = ln(OR)
```
# In R
beta <- log(OR)

# References:
#   - https://www.biostars.org/p/276869/
#   - https://www.biostars.org/p/336405/
```

## from Z-score
```
Beta = z / sqrt(2p(1− p)(n + z^2))

## Notation
# p: minor allele frequency
# n: sample size

# References:
#   - https://www.biostars.org/p/319584/
#   - https://www.nature.com/articles/ng.3538 (SMR 논문 supplementary text에 설명 있음)
```

# Compute Odds ratio
## from Beta
```
or <- exp(beta)
```

# Compute 95% confidence interval
## 95% confidence interval of Beta
```
# In R
beta_lo_ci <- beta - 1.96 * se
beta_up_ci <- beta + 1.96 * se
```

# Compute Z-score
## from Beta and SE
```
# In R
Z <- beta/se
```
## from Beta and P-value
```
# In R
Z <- sign(beta) * sqrt(qchisq(P, df=1, lower=F))

# Reference:
#   - Genomic SEM (https://github.com/GenomicSEM/GenomicSEM/blob/master/R/sumstats_main.R)
```

# Compute P-value
## from Z-score
```
# In R
p <- pchisq(Z^2, df=1, lower=F)
p <- 2*pnorm(abs(Z), lower.tail=F) # From Genomic SEM Wiki

# Reference:
#   - https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS2.html
```

# Compute SE of Beta
## from Beta and P-value
```
# In R
tail <- 2 # 1 for one tailed and 2 for two tailed
se <- abs(beta/ qnorm(P/tail))

# Reference:
#   - https://www.biostars.org/p/431875/
```
## from Z-score
```
SE = 1 / sqrt(2p(1− p)(n + Z^2))

## Notation
# p: minor allele frequency
# n: sample size

# References:
#   - https://www.biostars.org/p/319584/
#   - https://www.nature.com/articles/ng.3538 (SMR 논문 supplementary text에 설명 있음)
```

## from Odds ratio and SE of Odds ratio (이 방법은 추천하지 않음)
"OR_up나 OR_lo가 동일한 final value를 돌려줘야 할 텐데 SE가 다름."  
Obtain 95% CI of OR and use it to compute SE of beta.
```
# In R
beta <- log(OR)

OR_up <- OR + (OR.se * 1.96)
OR_lo <- OR - (OR.se * 1.96)

beta_up <- log(OR_up)

# Since, beta_lo = beta - (beta.se * 1.96)
beta.se <- (beta - beta_lo) / 1.96

# References:
# - https://www.researchgate.net/post/How_can_I_calculate_beta_coefficient_and_its_error_from_Odds_Ratio_from_GWAS_summary_Statisitcs
# - https://www.biostars.org/p/276869/
```

# OR과 OR에 대한 SE가 있을 때, Beta와 Beta에 대한 SE 계산하기
```
# 방법 1

beta <- log(OR)
tail <- 2
se <- abs(beta/qnorm(P/tail))
Z <- beta/se

# Reference:
#   - https://www.biostars.org/p/431875/
```

```
# 방법 2

beta <- log(OR)
Z <- sign(beta) * sqrt(qchisq(P, df=1, lower=F)),
se <- beta/Z

# Reference:
#   - Genomic SEM (https://github.com/GenomicSEM/GenomicSEM/blob/master/R/sumstats_main.R)
```
