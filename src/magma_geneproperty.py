from src.util import *

import pandas as pd
import os
import csv
from statsmodels.sandbox.stats.multicomp import multipletests


def MAGMA_geneproperty_GTEx_SpecificTissue(result_dir, verbose=False):
    file = os.path.join(result_dir, "magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out")
    assert FileExists(file), "Place `magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out` in the result directory."

    logs = []
    logs.append("Running MAGMA property analysis of GTEx v8 specific-tissue to table...")

    logs.append("Column description: Beta and Beta.Std is the non-standardized and semi-standardized regression coefficient from generalized linear regression described in [1]. SE is the standard error of the non-standardized regression coefficient. P value is derived with a one-tailed Z-tests of Beta divided by SE (not adjusted for multiple testing).")
    logs.append("Source: Karlsson Linnér, R. et al. Multivariate analysis of 1.5 million people identifies genetic associations with traits related to self-regulation and addiction. Nat Neurosci 24, 1367–1376 (2021).")

    logs.append("First 6 lines in magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out:")

    with open(file, 'r') as f:
        lines = f.readlines()
        logs += ["\t" + l.strip() for l in lines[:6]]
        lines = lines[6:]
        lines = [row.strip() for row in lines]
        lines = [row.split(sep=" ") for row in lines]
        lines = [[l for l in row if l] for row in lines]
    
    df = pd.DataFrame(lines, columns=["X1", "Type", "N genes", "Beta", "Beta STD", "SE", "P-value", "Tissue"])
    
    df['P-value'] = df['P-value'].astype(float)

    logs.append("Number of tissues tested: {:,}".format(len(df)))
    logs.append("Bonferroni correction threshold: {}".format(0.05 / len(df)))

    def pval_sig(p):
        return 'Yes' if p < 0.05 else 'No'

    def bonferroni_pval_cal(p):
        return multipletests(p, alpha=0.05, method='bonferroni')[1]
    
    bonferroni_pval = bonferroni_pval_cal(df['P-value'].tolist())
    df['Bonferroni P-value'] = bonferroni_pval
    df['Bonferroni significant'] = df['Bonferroni P-value'].apply(lambda x: pval_sig(x))

    logs.append("Number of tissue passing Bonferroni: {:,}".format(len(df[df['Bonferroni significant'] == "Yes"])))

    def fdr_pval_cal(p):
        return multipletests(p, alpha=0.05, method='fdr_bh')[1]
    
    fdr_pval = fdr_pval_cal(df['P-value'].tolist())
    df['FDR P-value'] = fdr_pval
    df['FDR significant'] = df['FDR P-value'].apply(lambda x: pval_sig(x))

    df = df[["Tissue", "N genes", "Beta", "Beta STD", "SE", "P-value", "Bonferroni P-value", "Bonferroni significant", "FDR P-value", "FDR significant"]]

    logs.append("Number of tissue passing FDR: {:,}".format(len(df[df['FDR significant'] == "Yes"])))

    df.sort_values(by=["Bonferroni P-value"], inplace=True)

    df['Tissue'] = df['Tissue'].apply(lambda x: x.replace("_", " "))

    df.to_csv("MAGMA_geneproperty_GTExV8_specifictissue_table.csv", sep=",", index=False,
            quoting=csv.QUOTE_ALL)

    logs.append("")
    
    if verbose:
        [print(l) for l in logs]

    with open("MAGMA_geneproperty_GTExV8_specifictissue_table.log", 'w') as f:
        f.writelines([l+"\n" for l in logs])


def MAGMA_geneproperty_GTEx_GeneralTissue(result_dir, verbose=False):
    file = os.path.join(result_dir, "magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out")
    assert FileExists(file), "Place `magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out` in the result directory."

    logs = []
    logs.append("Running MAGMA property analysis of GTEx v8 general tissue to table...")

    logs.append("Column description: Beta and Beta.Std is the non-standardized and semi-standardized regression coefficient from generalized linear regression described in [1]. SE is the standard error of the non-standardized regression coefficient. P value is derived with a one-tailed Z-tests of Beta divided by SE (not adjusted for multiple testing).")
    logs.append("Source: Karlsson Linnér, R. et al. Multivariate analysis of 1.5 million people identifies genetic associations with traits related to self-regulation and addiction. Nat Neurosci 24, 1367–1376 (2021).")

    logs.append("First 6 lines in magma_exp_gtex_v8_ts_general_avg_log2TPM.gsa.out:")

    with open(file, 'r') as f:
        lines = f.readlines()
        logs += ["\t" + l.strip() for l in lines[:6]]
        lines = lines[6:]
        lines = [row.strip() for row in lines]
        lines = [row.split(sep=" ") for row in lines]
        lines = [[l for l in row if l] for row in lines]
    
    
    df = pd.DataFrame(lines, columns=["Tissue", "Type", "N genes", "Beta", "Beta STD", "SE", "P-value"])
    df['P-value'] = df['P-value'].astype(float)

    logs.append("Number of tissues tested: {:,}".format(len(df)))
    logs.append("Bonferroni correction threshold: {}".format(0.05 / len(df)))

    def pval_sig(p):
        return 'Yes' if p < 0.05 else 'No'

    def bonferroni_pval_cal(p):
        return multipletests(p, alpha=0.05, method='bonferroni')[1]
    
    bonferroni_pval = bonferroni_pval_cal(df['P-value'].tolist())
    df['Bonferroni P-value'] = bonferroni_pval
    df['Bonferroni significant'] = df['Bonferroni P-value'].apply(lambda x: pval_sig(x))

    logs.append("Number of tissue passing Bonferroni: {:,}".format(len(df[df['Bonferroni significant'] == "Yes"])))

    def fdr_pval_cal(p):
        return multipletests(p, alpha=0.05, method='fdr_bh')[1]
    
    fdr_pval = fdr_pval_cal(df['P-value'].tolist())
    df['FDR P-value'] = fdr_pval
    df['FDR significant'] = df['FDR P-value'].apply(lambda x: pval_sig(x))

    df = df[["Tissue", "N genes", "Beta", "Beta STD", "SE", "P-value", "Bonferroni P-value", "Bonferroni significant", "FDR P-value", "FDR significant"]]

    logs.append("Number of tissue passing FDR: {:,}".format(len(df[df['FDR significant'] == "Yes"])))

    df.sort_values(by=["Bonferroni P-value"], inplace=True)

    df['Tissue'] = df['Tissue'].apply(lambda x: x.replace("_", " "))

    df.to_csv("MAGMA_geneproperty_GTExV8_generaltissue_table.csv", sep=",", index=False,
            quoting=csv.QUOTE_ALL)

    logs.append("")
    
    if verbose:
        [print(l) for l in logs]

    with open("MAGMA_geneproperty_GTExV8_generaltissue_table.log", 'w') as f:
        f.writelines([l+"\n" for l in logs])


def MAGMA_geneproperty_BrainSpan_Age(result_dir, verbose=False):
    file = os.path.join(result_dir, "magma_exp_bs_age_avg_log2RPKM.gsa.out")
    assert FileExists(file), "Place `magma_exp_bs_age_avg_log2RPKM.gsa.out` in the result directory."

    logs = []
    logs.append("Running MAGMA property analysis of BrainSpan age to table...")

    logs.append("Column description: Beta and Beta.Std is the non-standardized and semi-standardized regression coefficient from generalized linear regression described in [1]. SE is the standard error of the non-standardized regression coefficient. P value is derived with a one-tailed Z-tests of Beta divided by SE (not adjusted for multiple testing).")
    logs.append("Source: Karlsson Linnér, R. et al. Multivariate analysis of 1.5 million people identifies genetic associations with traits related to self-regulation and addiction. Nat Neurosci 24, 1367–1376 (2021).")

    logs.append("First 6 lines in magma_exp_bs_age_avg_log2RPKM.gsa.out:")

    with open(file, 'r') as f:
        lines = f.readlines()
        logs += ["\t" + l.strip() for l in lines[:6]]
        lines = lines[6:]
        lines = [row.strip() for row in lines]
        lines = [row.split(sep=" ") for row in lines]
        lines = [[l for l in row if l] for row in lines]
    
    df = pd.DataFrame(lines, columns=["Brain age", "Type", "N genes", "Beta", "Beta STD", "SE", "P-value"])
    
    df['P-value'] = df['P-value'].astype(float)

    logs.append("Number of Brain age tested: {:,}".format(len(df)))
    logs.append("Bonferroni correction threshold: {}".format(0.05 / len(df)))

    def pval_sig(p):
        return 'Yes' if p < 0.05 else 'No'

    def bonferroni_pval_cal(p):
        return multipletests(p, alpha=0.05, method='bonferroni')[1]
    
    bonferroni_pval = bonferroni_pval_cal(df['P-value'].tolist())
    df['Bonferroni P-value'] = bonferroni_pval
    df['Bonferroni significant'] = df['Bonferroni P-value'].apply(lambda x: pval_sig(x))

    logs.append("Number of Brain age passing Bonferroni: {:,}".format(len(df[df['Bonferroni significant'] == "Yes"])))

    def fdr_pval_cal(p):
        return multipletests(p, alpha=0.05, method='fdr_bh')[1]
    
    fdr_pval = fdr_pval_cal(df['P-value'].tolist())
    df['FDR P-value'] = fdr_pval
    df['FDR significant'] = df['FDR P-value'].apply(lambda x: pval_sig(x))

    df = df[["Brain age", "N genes", "Beta", "Beta STD", "SE", "P-value", "Bonferroni P-value", "Bonferroni significant", "FDR P-value", "FDR significant"]]

    logs.append("Number of Brain age passing FDR: {:,}".format(len(df[df['FDR significant'] == "Yes"])))

    df.sort_values(by=["Bonferroni P-value"], inplace=True)

    df['Brain age'] = df['Brain age'].apply(lambda x: x.replace("_", " "))

    df.to_csv("MAGMA_geneproperty_BrainSpan_age_table.csv", sep=",", index=False,
            quoting=csv.QUOTE_ALL)

    logs.append("")
    
    if verbose:
        [print(l) for l in logs]

    with open("MAGMA_geneproperty_BrainSpan_age_table.log", 'w') as f:
        f.writelines([l+"\n" for l in logs])


def MAGMA_geneproperty_BrainSpan_Developmental(result_dir, verbose=False):
    file = os.path.join(result_dir, "magma_exp_bs_dev_avg_log2RPKM.gsa.out")
    assert FileExists(file), "Place `magma_exp_bs_dev_avg_log2RPKM.gsa.out` in the result directory."

    logs = []
    logs.append("Running MAGMA property analysis of BrainSpan developmental stages to table...")

    logs.append("Column description: Beta and Beta.Std is the non-standardized and semi-standardized regression coefficient from generalized linear regression described in [1]. SE is the standard error of the non-standardized regression coefficient. P value is derived with a one-tailed Z-tests of Beta divided by SE (not adjusted for multiple testing).")
    logs.append("Source: Karlsson Linnér, R. et al. Multivariate analysis of 1.5 million people identifies genetic associations with traits related to self-regulation and addiction. Nat Neurosci 24, 1367–1376 (2021).")

    logs.append("First 6 lines in magma_exp_bs_dev_avg_log2RPKM.gsa.out:")

    with open(file, 'r') as f:
        lines = f.readlines()
        logs += ["\t" + l.strip() for l in lines[:6]]
        lines = lines[6:]
        lines = [row.strip() for row in lines]
        lines = [row.split(sep=" ") for row in lines]
        lines = [[l for l in row if l] for row in lines]
    
    df = pd.DataFrame(lines, columns=["Brain developmental stage", "Type", "N genes", "Beta", "Beta STD", "SE", "P-value"])
    
    df['P-value'] = df['P-value'].astype(float)

    logs.append("Number of Brain developmental stage tested: {:,}".format(len(df)))
    logs.append("Bonferroni correction threshold: {}".format(0.05 / len(df)))

    def pval_sig(p):
        return 'Yes' if p < 0.05 else 'No'

    def bonferroni_pval_cal(p):
        return multipletests(p, alpha=0.05, method='bonferroni')[1]
    
    bonferroni_pval = bonferroni_pval_cal(df['P-value'].tolist())
    df['Bonferroni P-value'] = bonferroni_pval
    df['Bonferroni significant'] = df['Bonferroni P-value'].apply(lambda x: pval_sig(x))

    logs.append("Number of Brain developmental stage passing Bonferroni: {:,}".format(len(df[df['Bonferroni significant'] == "Yes"])))

    def fdr_pval_cal(p):
        return multipletests(p, alpha=0.05, method='fdr_bh')[1]
    
    fdr_pval = fdr_pval_cal(df['P-value'].tolist())
    df['FDR P-value'] = fdr_pval
    df['FDR significant'] = df['FDR P-value'].apply(lambda x: pval_sig(x))

    df = df[["Brain developmental stage", "N genes", "Beta", "Beta STD", "SE", "P-value", "Bonferroni P-value", "Bonferroni significant", "FDR P-value", "FDR significant"]]

    logs.append("Number of Brain developmental stage passing FDR: {:,}".format(len(df[df['FDR significant'] == "Yes"])))

    df.sort_values(by=["Bonferroni P-value"], inplace=True)

    df['Brain developmental stage'] = df['Brain developmental stage'].apply(lambda x: x.replace("_", " "))


    df.to_csv("MAGMA_geneproperty_BrainSpan_developmental_table.csv", sep=",", index=False,
            quoting=csv.QUOTE_ALL)

    logs.append("")
    
    if verbose:
        [print(l) for l in logs]

    with open("MAGMA_geneproperty_BrainSpan_developmental_table.log", 'w') as f:
        f.writelines([l+"\n" for l in logs])
