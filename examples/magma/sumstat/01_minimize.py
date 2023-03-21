import pandas as pd
import code
# code.interact(local=dict(globals(), **locals()))

file = "pgcAN2.2019-07.vcf.tsv"

with open(file, 'r') as f:
    rows = f.readlines()
    for row_idx, row in enumerate(rows):
        if '#' not in row:
            break

df = pd.read_csv(file, sep="\t", index_col=False, skiprows=row_idx)

# print(df.columns)

df = df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'BETA', 'SE', 'PVAL', 'NCAS', 'NCON']]
df.columns = ['CHR', 'POS', 'SNP', 'REF', 'ALT', 'BETA', 'SE', 'P', 'NCASE', 'NCONTROL']
df.to_csv("pgcAN2.2019-07.vcf.minimal.tsv", sep="\t", index=False)

# code.interact(local=dict(globals(), **locals()))