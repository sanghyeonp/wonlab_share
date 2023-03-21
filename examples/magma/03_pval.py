import pandas as pd

file = "./sumstat/pgcAN2.2019-07.vcf.minimal.tsv"

df = pd.read_csv(file, sep="\t", index_col=False)

# Effective N calculation based on https://www.nature.com/articles/nprot.2014.071
# Neff = 2 / (1/Ncase + 1/Ncontrol)

def effective_n(ncase, ncontrol):
    return 2 / ((1/ncase) + (1/ncontrol))

df['N'] = df.apply(lambda row: effective_n(row['NCASE'], row['NCONTROL']), axis=1)

df.fillna("NA")

df[['SNP', 'P', 'N']].to_csv("AN.pval", sep="\t", index=False, header=True)
