
# GWAS p-value calculation from Z-score
https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS2.html

```
file2 <- "/data1/inshik/Glaucoma_TLCPD_MR/open_gwas/input/continuous-5255-both_sexes-irnt.tsv.gz"

dff <- fread(file2, nThread=2)

dff_sub <- dff[, c('chr', 'pos', 'beta_EUR', 'se_EUR', 'neglog10_pval_EUR')]


dff_sub$z <- dff_sub$beta_EUR / dff_sub$se_EUR

dff_sub$p_cal <- pchisq(dff_sub$z^2, df = 1, lower = F)

dff_sub$neglog10_p_cal <- -log10(dff_sub$p_cal)
dff_sub$ln_p_cal <- log(dff_sub$p_cal)

head(dff_sub, 20)

```

```
# P-value from beta and se
2*pnorm(abs(beta/se), lower.tail=F)

pchisq((beta / se)^2, df = 1, lower = F)
```

```
# New SE from beta and P-value
abs(beta)/qnorm(`P-value`/2, lower.tail=F)
```

```
# beta and se calculation from z-score
Beta = z / sqrt(2p(1− p)(n + z^2)) and
SE =1 / sqrt(2p(1− p)(n + z^2))

## Notation
# p: minor allele frequency
# n: sample size

## Reference
# https://www.biostars.org/p/319584/
# https://www.nature.com/articles/ng.3538 (SMR 논문 supplementary text에 설명 있음)
```
