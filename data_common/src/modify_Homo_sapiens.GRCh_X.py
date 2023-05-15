import pandas as pd


def read_ensembl(ENSEMBL_GENE_INFO):
    df_gene_info_ = pd.read_csv(ENSEMBL_GENE_INFO, sep="\t", index_col=False, compression="gzip", skiprows=5,
                            names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'],
                            low_memory=False)
    
    def attribute_map(x):
        x = x.split(sep="; ")
        x_map = {x1.split(sep=" ")[0]:x1.split(sep=" ")[1] for x1 in x}
        return x_map

    df_gene_info_['Probe_Ensembl'] = df_gene_info_['attribute'].apply(lambda x: attribute_map(x)['gene_id'].replace('"', ''))
    df_gene_info_['gene_symbol'] = df_gene_info_['attribute'].apply(lambda x: attribute_map(x)['gene_name'].replace('"', ''))
    df_gene_info_ = df_gene_info_[df_gene_info_['feature'] == 'gene']

    df_gene_info = df_gene_info_[['Probe_Ensembl', 'gene_symbol', 'seqname', 'start', 'strand']]
    df_gene_info.columns = ['ensembl_gene', 'gene_symbol', 'chr_gene', 'TSS', 'strand']

    return df_gene_info


if __name__ == "__main__":
    ENSEMBL_GENE_INFO = "../Homo_sapiens.GRCh37.87.gtf.gz"
    df_gene_info = read_ensembl(ENSEMBL_GENE_INFO)
    df_gene_info.to_csv("../modified.Homo_sapiens.GRCh37.87.gtf.gz", sep="\t", index=False, compression="gzip")

    ENSEMBL_GENE_INFO = "../Homo_sapiens.GRCh38.109.gtf.gz"
    df_gene_info = read_ensembl(ENSEMBL_GENE_INFO)
    df_gene_info.to_csv("../modified.Homo_sapiens.GRCh38.109.gtf.gz", sep="\t", index=False, compression="gzip")