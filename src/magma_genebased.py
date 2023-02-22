from src.util import *

import pandas as pd
import os
import csv
from statsmodels.sandbox.stats.multicomp import multipletests


def MAGMA_genebased2Table(result_dir, snp2gene_table=None, verbose=False):
    genebased_file = os.path.join(result_dir, "magma.genes.out")
    assert FileExists(genebased_file), "Place `magma.genes.out` in the result directory."

    logs = []
    logs.append("Running MAGMA gene-based analysis to table...")

    df = pd.read_csv(genebased_file, sep="\t", index_col=False)
    
    logs.append("Number of genes tested: {:,}".format(len(df)))
    logs.append("Bonferroni correction threshold: {}".format(0.05 / len(df)))

    def pval_sig(p):
        return 'Yes' if p < 0.05 else 'No'

    def bonferroni_pval_cal(p):
        return multipletests(p, alpha=0.05, method='bonferroni')[1]
    
    bonferroni_pval = bonferroni_pval_cal(df['P'].tolist())
    df['Bonferroni P-value'] = bonferroni_pval
    df['Bonferroni significant'] = df['Bonferroni P-value'].apply(lambda x: pval_sig(x))

    logs.append("Number of genes passing Bonferroni: {:,}".format(len(df[df['Bonferroni significant'] == "Yes"])))

    df.rename(columns={"GENE":"Ensembl Gene ID",
                    "SYMBOL":"Gene name", 
                    "CHR":"Chromosome",
                    "START":"Start position",
                    "STOP":"End position",
                    "NSNPS": "N SNPs",
                    "NPARAM":"N Param",
                    "ZSTAT":"Z-score",
                    "P":"P-value",
                    }, inplace=True)

    df = df[["Ensembl Gene ID", "Gene name", "Chromosome", "Start position", "End position", "N SNPs", "N Param", "Z-score", "P-value", "Bonferroni P-value", "Bonferroni significant"]]

    df.sort_values(by=["Bonferroni P-value"], inplace=True)

    df.to_csv("MAGMA_genebased_table.csv", sep=",", index=False,
            quoting=csv.QUOTE_ALL)

    df[df["Bonferroni significant"] == "Yes"].to_csv("MAGMA_genebased_Bonferroni_table.csv", sep=",", index=False,
                                                        quoting=csv.QUOTE_ALL)

    if snp2gene_table is None and os.path.exists(os.path.join(os.getcwd(), "leadSNP_to_gene_mapping_table.csv")):
        snp2gene_table = os.path.join(os.getcwd(), "leadSNP_to_gene_mapping_table.csv")

    if snp2gene_table is not None:
        df2 = pd.read_csv(snp2gene_table, sep=",", index_col=False)
        magma_genebased_genes = df[df["Bonferroni significant"] == "Yes"]['Ensembl Gene ID'].tolist()
        snp2gene_genes = df2[(df2['posMap'] == 'Yes') & (df2['eqtlMap'] == 'Yes') & (df2['ciMap'] == 'Yes')]["Ensembl ID"].tolist()

        logs.append("Number of genes identified from 3 stragies of SNP2GENE and MAGMA gene-based analysis: {:,}".format(len(set(magma_genebased_genes).intersection(set(snp2gene_genes)))))

    logs.append("")

    if verbose:
        [print(l) for l in logs]
    
    with open("MAGMA_genebased_table.log", 'w') as f:
        f.writelines([l+"\n" for l in logs])
