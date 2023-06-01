
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