import pandas as pd
import os
import code
# code.interact(local=dict(globals(), **locals()))

dir_1kGp3 = "/data/public/1kG/phase3/plink"

abc_rsid_map_file = "1000g_rsid_chrbp"

df_map = pd.read_csv(os.path.join(dir_1kGp3, abc_rsid_map_file), sep=" ", index_col=False, names=['chrpos', 'rsid'])

eur_bim = "1kg.phase3.auto.snp.qc.chrposabc.EUR.clean.bim"

df_eur_bim = pd.read_csv(os.path.join(dir_1kGp3, eur_bim), sep="\t", index_col=False, names=['CHR', 'abc', 'CM', 'POS', 'A1', 'A2'])
df_eur_bim['chrpos'] = df_eur_bim['abc'].apply(lambda x: ':'.join(x.split(sep=":")[:2]))

df_eur_bim = pd.merge(df_eur_bim, df_map, how="left", on="chrpos")

"""
# Mapping 파일에 같은 chr:pos에 여러 rsID가 있는 것 같음.
# dbSNP에 확인해본 결과 
#   7:21455345 => rs39290
#   7:75706634 => rs10274404
>>> df_eur_bim.loc[df_eur_bim.duplicated(subset=['chrpos'], keep=False)]
         CHR             abc  CM       POS A1 A2      chrpos        rsid
3504337    7  7:21455345:A:G   0  21455345  A  G  7:21455345     rs39290
3504338    7  7:21455345:A:G   0  21455345  A  G  7:21455345  rs58407843
3659461    7  7:75706634:A:G   0  75706634  A  G  7:75706634  rs10274404
3659462    7  7:75706634:A:G   0  75706634  A  G  7:75706634  rs67988737
"""

df_eur_bim.drop([3504338, 3659462], axis=0, inplace=True)

df_eur_bim.drop(columns=['abc', 'chrpos'], inplace=True)
df_eur_bim1 = df_eur_bim[['CHR', 'rsid', 'CM', 'POS', 'A1', 'A2']]

df_eur_bim1['rsid'] = df_eur_bim1['rsid'].fillna('.')

df_eur_bim1.to_csv("./EUR/1kg.phase3.auto.snp.qc.chrposabc.EUR.clean.rsid.new.bim", sep="\t", index=False, header=False)
del df_eur_bim1

df_eur_bim['chrpos'] = df_eur_bim.apply(lambda row: "{}:{}".format(row['CHR'], row['POS']), axis=1)
df_eur_bim2 = df_eur_bim[['CHR', 'chrpos', 'CM', 'POS', 'A1', 'A2']]
df_eur_bim2.to_csv("./EUR/1kg.phase3.auto.snp.qc.chrposabc.EUR.clean.chrpos.new.bim", sep="\t", index=False, header=False)
del df_eur_bim; del df_eur_bim2

#################################################################
eas_bim = "1kg.phase3.auto.snp.qc.chrposabc.EAS.clean.bim"

df_eas_bim = pd.read_csv(os.path.join(dir_1kGp3, eas_bim), sep="\t", index_col=False, names=['CHR', 'abc', 'CM', 'POS', 'A1', 'A2'])
df_eas_bim['chrpos'] = df_eas_bim['abc'].apply(lambda x: ':'.join(x.split(sep=":")[:2]))
df_eas_bim = pd.merge(df_eas_bim, df_map, how="left", on="chrpos")

"""
# Mapping 파일에 같은 chr:pos에 여러 rsID가 있는 것 같음.
# dbSNP에 확인해본 결과 
#   7:21455345 => rs39290
#   7:75706634 => rs10274404
>>> df_eas_bim.loc[df_eas_bim.duplicated(subset=['chrpos'], keep=False)]
         CHR             abc  CM       POS A1 A2      chrpos        rsid
3087446    7  7:21455345:A:G   0  21455345  A  G  7:21455345     rs39290
3087447    7  7:21455345:A:G   0  21455345  A  G  7:21455345  rs58407843
3226529    7  7:75706634:A:G   0  75706634  G  A  7:75706634  rs10274404
3226530    7  7:75706634:A:G   0  75706634  G  A  7:75706634  rs67988737
"""

df_eas_bim.drop([3087447, 3226530], axis=0, inplace=True)

df_eas_bim.drop(columns=['chrpos', 'abc'], inplace=True)
df_eas_bim1 = df_eas_bim[['CHR', 'rsid', 'CM', 'POS', 'A1', 'A2']]

df_eas_bim1['rsid'] = df_eas_bim1['rsid'].fillna('.')

df_eas_bim1.to_csv("./EAS/1kg.phase3.auto.snp.qc.chrposabc.EAS.clean.rsid.new.bim", sep="\t", index=False, header=False)

df_eas_bim['chrpos'] = df_eas_bim.apply(lambda row: "{}:{}".format(row['CHR'], row['POS']), axis=1)
df_eas_bim2 = df_eas_bim[['CHR', 'chrpos', 'CM', 'POS', 'A1', 'A2']]
df_eas_bim2.to_csv("./EAS/1kg.phase3.auto.snp.qc.chrposabc.EAS.clean.chrpos.new.bim", sep="\t", index=False, header=False)

