

# SNP

## Structural variant

https://gatk.broadinstitute.org/hc/en-us/articles/9022476791323-Structural-Variants

## Bi-allelic SNP vs Multi-allelic SNP vs Variant 

https://gatk.broadinstitute.org/hc/en-us/articles/360035890771-Biallelic-vs-Multiallelic-sites

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3198575/

```
# Example: multi-allelic
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
22	16051453	.	A	C,G	100	PASS	AC=478,17;AF=0.0954473,0.00339457;AN=5008;NS=2504;DP=22548;EAS_AF=0.0744,0;AMR_AF=0.1239,0;AFR_AF=0.003,0;EUR_AF=0.0746,0.003;SAS_AF=0.2434,0.0143;AA=.|||;VT=SNP;MULTI_ALLELIC

# Example: bi-allelic
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
0	22	16050075	.	A	G	100	PASS	AC=1;AF=0.000199681;AN=5008;NS=2504;DP=8012;EA

# Example: variant
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
22	16052167	.	A	AAAAC	100	PASS	AC=2367;AF=0.472644;AN=5008;NS=2504;DP=22414;EAS_AF=0.6687;AMR_AF=0.5548;AFR_AF=0.357;EUR_AF=0.3698;SAS_AF=0.4744;VT=INDEL	True

```