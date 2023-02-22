from util import *

import pandas as pd
import os
import csv

def leadSNP2GeneMapping(result_dir, verbose=False):
	snp2gene_file = os.path.join(result_dir, "genes.txt")
	assert FileExists(snp2gene_file), "Place `genes.txt` in the result directory."

	logs = []
	logs.append("Running lead SNP to gene mapping to table...")

	df = pd.read_csv(snp2gene_file, sep="\t", index_col=False)

	df = df[['ensg', 'symbol', 'type', 'chr', 'start', 'end', 'ciMap', 'posMapSNPs', 'posMapMaxCADD', 'eqtlMapSNPs', 'eqtlMapminP', 'eqtlMapminQ', 'eqtlMapts', 'ciMapts']]

	df.rename(columns={'ensg':'Ensembl ID', 
					'symbol':'Gene', 
					'type': 'Gene type', 
					'chr': 'Chr', 
					'start': 'Start', 
					'end': 'End', 
					'ciMap': 'ciMap', 
					'posMapSNPs': 'N posMap', 
					'posMapMaxCADD': 'Max CADD posMap', 
					'eqtlMapSNPs': 'N eqtlMap', 
					'eqtlMapminP': 'P eqtlMap', 
					'eqtlMapts': 'Tissue eqtlMap', 
					'ciMapts': 'Tissue ciMap'}, inplace=True)

	df['posMap'] = df['N posMap'].apply(lambda x: 'Yes' if x > 0 else 'No')
	df['eqtlMap'] = df['N eqtlMap'].apply(lambda x: 'Yes' if x > 0 else 'No')

	df = df[['Ensembl ID', 'Gene', 'Gene type', 'Chr', 'Start', 'End', 'posMap', 'eqtlMap', 'ciMap',
		'N posMap', 'Max CADD posMap', 'N eqtlMap', 'P eqtlMap', 'eqtlMapminQ', 'Tissue eqtlMap',
		'Tissue ciMap']]


	logs.append("Number of mapped gene by positional mapping: {:,}".format(len(df[df['posMap'] == 'Yes'])))
	logs.append("Number of mapped gene by eQTL mapping: {:,}".format(len(df[df['eqtlMap'] == 'Yes'])))
	logs.append("Number of mapped gene by chromatin interaction mapping: {:,}".format(len(df[df['ciMap'] == 'Yes'])))

	df.fillna("-", inplace=True)
	df['Tissue eqtlMap'] = df['Tissue eqtlMap'].apply(lambda x: x.replace(":", ", "))
	df['Tissue eqtlMap'] = df['Tissue eqtlMap'].apply(lambda x: x.replace("_", " "))
	df['Tissue ciMap'] = df['Tissue ciMap'].apply(lambda x: x.replace(":", ", "))
	df['Tissue ciMap'] = df['Tissue ciMap'].apply(lambda x: x.replace("_", " "))


	df.to_csv("leadSNP_to_gene_mapping_table.csv", sep=",", index=False, quoting=csv.QUOTE_ALL)
	
	logs.append("Number of genes identified in all 3 strategies: {:,}".format(len(df[(df['posMap'] == 'Yes') & (df['eqtlMap'] == 'Yes') & (df['ciMap'] == 'Yes')])))
	logs.append("")

	if verbose:
		[print(l) for l in logs]
	
	with open("leadSNP_to_gene_mapping_table.log", 'w') as f:
		f.writelines([l+"\n" for l in logs])
