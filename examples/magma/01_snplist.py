import pandas as pd

file = "./sumstat/pgcAN2.2019-07.vcf.minimal.tsv"

df = pd.read_csv(file, sep="\t", index_col=False)

df[['SNP', 'CHR', 'POS']].to_csv("AN.snp.loc", sep="\t", index=False, header=False)
