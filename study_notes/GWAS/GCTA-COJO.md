
- GCTA-COJO is used to identify SNPs that are "truly" independent.
- COJO can be ran by 
    1. feeding the whole GWAS summary statistics and identify independent SNPS
    2. feeding the clumped SNPs (e.g., using PLINK) and identify independent SNPs among the lead SNPs.

- Second approach was used in the below references:
    - https://www.nature.com/articles/s41467-023-36013-1
    - https://www.nature.com/articles/s41593-021-00908-3

- Second approach can be done by using `--extract` parameter.

- Anyway, both approach will result in three output files:
    1. `*.jma.cojo`: file with independent associations identified in the stepwise model selection
    2. `*.ldr.cojo`: file showing the LD correlations between the SNPs
    3. `*.cma.cojo`: file with conditional analysis results for all other SNPs that were not selected as independent associations
    - reference: https://cnsgenomics.com/data/teaching/GNGWS23/module1/9_independentLociPrac.html

- Basically, all SNPs in `*.cma.cojo` will be the SNPs that did not PASS conditional analysis, hence considered as SNPs that are not independent. We have to proceed with further analysis with the SNPs in `*.jma.cojo` as these SNPs are significant after conditional analysis, thus can be considered as independent SNPs.
- It can be observed by looking at the `pC` and `pJ` columns in `*.cma.cojo` and `*.jma.cojo` files respectively, which represents the P-value of SNP after conditional analysis.

# Reference
- Result from Clumping vs COJO  
https://cnsgenomics.com/data/teaching/GNGWS23/module1/9_independentLociPrac.html