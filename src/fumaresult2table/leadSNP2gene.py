import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *

def leadSNP2GeneMapping(result_dir, outd, max_dis = 10, verbose=False):
	"""
	max_dis: in kb
	"""
	max_dis = max_dis * 1000

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
					'eqtlMapminP': 'min P eqtlMap', 
					'eqtlMapts': 'Tissue eqtlMap', 
					'ciMapts': 'Tissue ciMap'}, inplace=True)

	df['posMap'] = df['N posMap'].apply(lambda x: 'Yes' if x > 0 else 'No')
	df['eqtlMap'] = df['N eqtlMap'].apply(lambda x: 'Yes' if x > 0 else 'No')

	## Match which SNP the gene was mapped to
	# positional mapping: `N posMap`
	snps_file = os.path.join(result_dir, "snps.txt")
	df_snps = pd.read_csv(snps_file, sep="\t", index_col=False)
	df_snps = df_snps[['uniqID', 'rsID', 'chr', 'pos', 'posMapFilt']]

	def match_snps4posmap(df_snps, chr, start, end, max_dis):
		df_snps_subset = df_snps[df_snps['chr'] == chr]
		df_snps_subset = df_snps_subset[df_snps_subset['pos'].between(start - max_dis, end + max_dis)]
		if df_snps_subset.empty:
			return '-'
		snps = df_snps_subset[df_snps_subset['posMapFilt'] == 1]['rsID'].tolist()
		return ', '.join(snps)

	df['SNPs posMap'] = df.apply(lambda row: match_snps4posmap(df_snps, row['Chr'], row['Start'], row['End'], max_dis), axis=1)
	
	# eqtl mapping: 
	eqtl_file = os.path.join(result_dir, "eqtl.txt")
	df_eqtl = pd.read_csv(eqtl_file, sep="\t", index_col=False)
	df_eqtl = df_eqtl[['uniqID', 'symbol', 'eqtlMapFilt']]

	def match_snps4eqtlmap(df_eqtl, df_snps, gene):
		df_eqtl_subset = df_eqtl[(df_eqtl['symbol'] == gene) & (df_eqtl['eqtlMapFilt'] == 1)]
		if df_eqtl_subset.empty:
			return '-'
		snps = df_eqtl_subset['uniqID'].tolist()
		df_snps = df_snps[df_snps['uniqID'].isin(snps)][['uniqID', 'rsID']]
		df_eqtl_subset = df_eqtl_subset.merge(df_snps, how="left", on='uniqID')
		df_eqtl_subset.fillna({'rsID':'-9'}, inplace=True)
		df_eqtl_subset['rsID'] = df_eqtl_subset.apply(lambda row: row['uniqID'] if row['rsID'] == '-9' else row['rsID'], axis=1)
		return ', '.join(df_eqtl_subset['rsID'].tolist())

	df['SNPs eqtlMap'] = df.apply(lambda row: match_snps4eqtlmap(df_eqtl, df_snps, row['Gene']), axis=1)

	# chromatin interaction mapping: 
	ci_file = os.path.join(result_dir, 'ci.txt')
	df_ci = pd.read_csv(ci_file, sep="\t", index_col=False)
	df_ci = df_ci[['SNPs', 'genes', 'ciMapFilt']]
	df_ci.fillna({'genes':'NA'}, inplace=True)

	def match_snps4cimap(df_ci, gene_id):
		df_ci_subset = df_ci[(df_ci['genes'].str.contains(gene_id)) & (df_ci['ciMapFilt'] == 1)]
		if df_ci_subset.empty:
			return '-'
		snps = [v.split(sep=":") for v in df_ci_subset['SNPs'].tolist()]
		snps = list(set([item for sublist in snps for item in sublist]))
		return ', '.join(snps)

	df['SNPs ciMap'] = df.apply(lambda row: match_snps4cimap(df_ci, row['Ensembl ID']), axis=1)

	## Re-order columns
	df = df[['Ensembl ID', 'Gene', 'Gene type', 'Chr', 'Start', 'End', 'posMap', 'eqtlMap', 'ciMap',
		'N posMap', 'Max CADD posMap', 'N eqtlMap', 'min P eqtlMap', 'eqtlMapminQ', 'Tissue eqtlMap',
		'Tissue ciMap', 'SNPs posMap', 'SNPs eqtlMap', 'SNPs ciMap']]

	logs.append("Number of mapped gene by positional mapping: {:,}".format(len(df[df['posMap'] == 'Yes'])))
	logs.append("Number of mapped gene by eQTL mapping: {:,}".format(len(df[df['eqtlMap'] == 'Yes'])))
	logs.append("Number of mapped gene by chromatin interaction mapping: {:,}".format(len(df[df['ciMap'] == 'Yes'])))

	df.fillna("-", inplace=True)
	df['Tissue eqtlMap'] = df['Tissue eqtlMap'].apply(lambda x: x.replace(":", ", "))
	df['Tissue eqtlMap'] = df['Tissue eqtlMap'].apply(lambda x: x.replace("_", " "))
	df['Tissue ciMap'] = df['Tissue ciMap'].apply(lambda x: x.replace(":", ", "))
	df['Tissue ciMap'] = df['Tissue ciMap'].apply(lambda x: x.replace("_", " "))


	df.to_csv(os.path.join(outd, "leadSNP_to_gene_mapping_table.csv"), sep=",", index=False, quoting=csv.QUOTE_ALL)
	
	logs.append("Number of genes identified in all 3 strategies: {:,}".format(len(df[(df['posMap'] == 'Yes') & (df['eqtlMap'] == 'Yes') & (df['ciMap'] == 'Yes')])))
	logs.append("")

	if verbose:
		[print(l) for l in logs]
	
	with open(os.path.join(outd, "leadSNP_to_gene_mapping_table.log"), 'w') as f:
		f.writelines([l+"\n" for l in logs])
