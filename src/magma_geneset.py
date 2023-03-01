from src.util import *

import pandas as pd
import os
import csv
from statsmodels.sandbox.stats.multicomp import multipletests
import code
# code.interact(local=dict(globals(), **locals()))

def MAGMA_geneset2Table(result_dir, outd, verbose=False):
    geneset_file = os.path.join(result_dir, "magma.gsa.out")
    assert FileExists(geneset_file), "Place `magma.gsa.out` in the result directory."

    logs = []
    logs.append("Running MAGMA gene-set analysis to table...")
    logs.append("Column description: Beta and Beta.Std is the non-standardized and semi-standardized regression coefficient from generalized linear regression described in [1]. SE is the standard error of the non-standardized regression coefficient. P value is derived with a one-tailed Z-tests of Beta divided by SE (not adjusted for multiple testing).")
    logs.append("Source: Karlsson Linnér, R. et al. Multivariate analysis of 1.5 million people identifies genetic associations with traits related to self-regulation and addiction. Nat Neurosci 24, 1367–1376 (2021).")

    

    with open(geneset_file, 'r') as f:
        lines = f.readlines()
        idx2start = len([idx for idx, l in enumerate(lines) if '#' in l]) + 1
        logs_ = ["\t" + l.strip() for l in lines[:idx2start]]
        lines = lines[idx2start:]
        lines = [row.strip() for row in lines]
        lines = [row.split(sep=" ") for row in lines]
        lines = [[l for l in row if l] for row in lines]
    
    logs.append("First {} lines in magma.gsa.out:".format(idx2start))
    logs += logs_

    df = pd.DataFrame(lines, columns=["X1", "Type", "N genes", "Beta", "Beta STD", "SE", "P-value", "Gene set"])

    df['P-value'] = df['P-value'].astype(float)

    logs.append("Number of gene-sets tested: {:,}".format(len(df)))
    logs.append("Bonferroni correction threshold: {}".format(0.05 / len(df)))

    gene_set_list = df['Gene set'].tolist()
    def counter(list1, string1):
        count1 = 0
        for v in list1:
            if string1 in v:
                count1 += 1
        return count1

    logs.append("Number of curated gene-sets: {:,}".format(counter(gene_set_list, "Curated_gene_sets:")))
    logs.append("Number of GO gene-sets: {:,}".format(counter(gene_set_list, "GO_")))

    def pval_sig(p):
        return 'Yes' if p < 0.05 else 'No'

    def bonferroni_pval_cal(p):
        return multipletests(p, alpha=0.05, method='bonferroni')[1]
    
    bonferroni_pval = bonferroni_pval_cal(df['P-value'].tolist())
    df['Bonferroni P-value'] = bonferroni_pval
    df['Bonferroni significant'] = df['Bonferroni P-value'].apply(lambda x: pval_sig(x))

    logs.append("Number of gene-sets passing Bonferroni: {:,}".format(len(df[df['Bonferroni significant'] == "Yes"])))

    df = df[["Gene set", "N genes", "Beta", "Beta STD", "SE", "P-value", "Bonferroni P-value", "Bonferroni significant"]]

    df.sort_values(by=["Bonferroni P-value"], inplace=True)

    df['Gene set'] = df['Gene set'].apply(lambda x: x.replace("_", " "))

    df.to_csv(os.path.join(outd, "MAGMA_geneset_table.csv"), sep=",", index=False,
            quoting=csv.QUOTE_ALL)

    df[df["Bonferroni significant"] == "Yes"].to_csv(os.path.join(outd, "MAGMA_geneset_Bonferroni_table.csv"), sep=",", index=False,
            quoting=csv.QUOTE_ALL)
    
    logs.append("")
    
    if verbose:
        [print(l) for l in logs]

    with open(os.path.join(outd, "MAGMA_geneset_table.log"), 'w') as f:
        f.writelines([l+"\n" for l in logs])