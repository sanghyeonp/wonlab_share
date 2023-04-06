from packages import *

# NCBI Entrez gene ID to gene name mapping file
GENEINFO = "/data1/sanghyeon/wonlab_contribute/combined/software/MAGMA/resources/gene_info.gz"

def read_geneinfo():
    print("Reading gene_info.gz...")
    df = pd.read_csv(GENEINFO, sep="\t", index_col=False, compression='gzip', low_memory=False)

    df = df[['GeneID', 'Symbol']]
    df.columns = ['GeneID', 'GeneName']
    return df

# NCBI dbSNP rsID mapping file
rsID_MAP_dbSNP_GRCh37 = "/data1/sanghyeon/Projects/MetabolicSyndrome/commonFiles/reference/NCBI_dbsnp/human_9606_b151_GRCh37p13.tsv.gz"

rsID_MAP_dbSNP = {37: rsID_MAP_dbSNP_GRCh37}