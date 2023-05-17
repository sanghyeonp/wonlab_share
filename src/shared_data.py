from packages import *

### NCBI Entrez gene ID to gene name mapping file
GENEINFO = "/data1/sanghyeon/wonlab_contribute/combined/data_common/gene_info.gz"

def read_geneinfo():
    df = pd.read_csv(GENEINFO, sep="\t", index_col=False, compression='gzip', low_memory=False)

    df = df[['GeneID', 'Symbol']]
    df.columns = ['GeneID', 'GeneName']
    return df

### NCBI dbSNP rsID mapping file
rsID_MAP_dbSNP_GRCh37 = "/data1/sanghyeon/wonlab_contribute/combined/data_common/human_9606_b151_GRCh37p13.tsv.gz"
rsID_MAP_dbSNP_GRCh38 = "/data1/sanghyeon/wonlab_contribute/combined/data_common/human_9606_b151_GRCh38p7.tsv.gz"
rsID_MAP_dbSNP = {37: rsID_MAP_dbSNP_GRCh37,
                38: rsID_MAP_dbSNP_GRCh38
                }


### ANNOVAR
ANNOVAR_software = "/data1/sanghyeon/wonlab_contribute/combined/software/annovar/annovar/table_annovar.pl"
ANNOVAR_humandb = "/data1/sanghyeon/wonlab_contribute/combined/software/annovar/annovar/humandb/"
# ANNOVAR_humandb = "/data/software/Annovar/annovar_20221005/humandb/"

### Ensembl gene information
ENSEMBL_GENE_INFO = "/data1/sanghyeon/wonlab_contribute/combined/data_common/Homo_sapiens.GRCh37.87.gtf.gz"

### NCBI variant annotation
NCBI_VARIANT_ANNOT = {37: "/data1/sanghyeon/wonlab_contribute/combined/data_common/00_All_b37.vcf.parquet.gz",
                    38: "/data1/sanghyeon/wonlab_contribute/combined/data_common/00_All_b38.vcf.parquet.gz"
                    }