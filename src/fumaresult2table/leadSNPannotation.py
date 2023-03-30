import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *


def leadSNPannotation(result_dir, outd, verbose=False):
    logs = []
    logs.append("Running lead SNP annotation to table...")

    # leadSNPs.txt 읽기
    leadsnp_result = os.path.join(result_dir, "leadSNPs.txt")
    assert FileExists(leadsnp_result), "Place `leadSNPs.txt` in the result directory."

    df_leadsnp = pd.read_csv(leadsnp_result, sep="\t", index_col=False)

    # 필요한 column들만 extract: rsID, uniqID
    df = df_leadsnp[['GenomicLocus', 'rsID', 'uniqID']].copy()

    # snps.txt 읽기
    snps_result = os.path.join(result_dir, "snps.txt")
    assert FileExists(snps_result), "Place `snps.txt` in the result directory."

    df_snps = pd.read_csv(snps_result, sep="\t", index_col=False, low_memory=False)

    # 필요한 column들만 extract: rsID, chr, pos, effect_allele, non_effect_allele, MAF, gwasP, beta, se, nearestGene, dist, func, CADD, RDB, minChrState, commonChrState
    df_snps = df_snps[["rsID", "chr", "pos", "effect_allele", "non_effect_allele", "MAF", "gwasP", "beta", "se", "nearestGene", "func", "CADD", "RDB", "minChrState", "commonChrState"]].copy()

    # Merge with main dataframe
    df = pd.merge(df, df_snps, how="left", on="rsID")

    # gwascatalog.txt 읽기
    catalog_result = os.path.join(result_dir, "gwascatalog.txt")
    assert FileExists(catalog_result), "Place `gwascatalog.txt` in the result directory."

    df_catalog = pd.read_csv(catalog_result, sep="\t", index_col=False)

    # 필요한 column들만 extract: snp, PMID, Trait
    df_catalog = df_catalog[['snp', "PMID", "Trait"]]

    # Duplicated SNPs는 하나로 합치기
    data = dict() # {snp: [[], []]}
    for _, row in tqdm(df_catalog.iterrows(), leave=False):
        snp, pmid, trait = row['snp'], row['PMID'], row['Trait']
        if snp in data:
            data[snp][0].append(pmid)
            data[snp][1].append(trait)
        else:
            data[snp] = [[pmid], [trait]]

    data_rows = []
    for snp, v in tqdm(data.items(), leave=False):
        data_rows.append([snp, ", ".join([str(e) for e in v[0]]), ", ".join(v[1])])

    df_catalog_nodup = pd.DataFrame(data_rows, columns=['rsID', "PMID", "Trait"])

    # Merge with main dataframe
    df = pd.merge(df, df_catalog_nodup, how="left", on="rsID")

    # Rename columns
    df.rename(columns={'GenomicLocus':'Genomic locus', 'chr':"Chr", 'pos':"Pos", 'effect_allele':"EA", 'non_effect_allele':"NEA",
        'gwasP':"P-value", 'beta':"Beta", 'se':"SE", 'nearestGene':"Nearest gene", 'func':"Functional consequence of SNP on gene", 
        'RDB':"RegulomeDB score", 'minChrState':"Minimum chromatin state", 'commonChrState':"Commmon chromatin state", 
        'PMID':"PUBMED ID"}, inplace=True)

    # Match functional consequence on the genes for exonic SNPs
    file_annov = os.path.join(result_dir, "annov.txt")
    df_annov = pd.read_csv(file_annov, sep="\t", index_col=False)
    df_annov_sub = df_annov[['uniqID', 'exonic_func']]
    df_annov_sub.columns = ['uniqID', 'Functional consequence of exonic SNP on gene']
    
    df = df.merge(df_annov_sub, how="left", on="uniqID")

    # Re-order columns
    df = df[['Genomic locus', 'rsID', 'uniqID', 'Chr', 'Pos', 'EA', 'NEA', 'MAF',
            'P-value', 'Beta', 'SE', 'Nearest gene', 'Functional consequence of SNP on gene', 
            'Functional consequence of exonic SNP on gene', 'CADD',
            'RegulomeDB score', 'Minimum chromatin state',
            'Commmon chromatin state', 'PUBMED ID', 'Trait']]

    # Replace Nan to -
    df.fillna("-")

    # Save the result
    df.to_csv(os.path.join(outd, "leadSNP_annotation_table.csv"), sep=",", index=False, quoting=csv.QUOTE_ALL)

    logs.append("Table saved at: {}".format(os.path.join(outd, "leadSNP_annotation_table.csv")))
    logs.append("")

    if verbose:
        [print(l) for l in logs]
