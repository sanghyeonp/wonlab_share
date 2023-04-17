import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from packages import *
from util import map_delim
from shared_data import ENSEMBL_GENE_INFO

def parse_args():
    parser = argparse.ArgumentParser(description=":: Map Ensembl gene ID to gene symbol, TSS, and strand ::")
    parser.add_argument('--file', required=True,
                        help='Path to the input file.')
    parser.add_argument('--gene-id-col', dest="gene_id_col", required=True, 
                        help='Column name of Ensembl gene ID in the input file.')
    parser.add_argument('--delim-in', dest="delim_in", required=False, type=str, default="tab",
                        help='Delimiter of the input file. Default = "tab". Choices = ["tab", "comma", "whitespace", "formatted"].')
    parser.add_argument('--compression-in', dest="compression_in", required=False, default="infer",
                        help="Specify compression type from the following ['zip', 'gzip', 'bz2', 'zstd', 'tar']. Default='infer'.")
    parser.add_argument('--skip-rows', dest="skip_rows", required=False, type=int, default=0,
                        help='Specify the number of first lines in the input file to skip. Default = 0.')

    parser.add_argument('--gene-symbol-only', dest="gene_symbol_only", action='store_true',
                    help='Specify to map only the gene symbol. Default = False.')

    parser.add_argument('--outf', required=False, type=str, default='NA',
                        help='Name of output file. Default = geneannot.<Input file name>')
    parser.add_argument('--outd', required=False, type=str, default='NA',
                        help='Directory where the output will be saved. Default = Current working directory.')
    parser.add_argument('--delim-out', dest="delim_out", required=False, type=str, default="NA",
                        help='Delimiter of the output file. Default = "NA", then identical delimiter as the input file. \
                            Choices = ["tab", "comma", "whitespace", "formatted"].')

    args = parser.parse_args()
    return args


def read_ensembl(gene_symbol_only):
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
    if gene_symbol_only:
        df_gene_info = df_gene_info[['ensembl_gene', 'gene_symbol']]

    return df_gene_info


def main(file, delim_in, compression_in, gene_id_col, gene_symbol_only, outf, outd, delim_out):
    df_ = pd.read_csv(file, sep=delim_in, index_col=False, compression=compression_in, low_memory=False)

    df_gene_info = read_ensembl(gene_symbol_only)

    ### Mapping
    df = df_.merge(df_gene_info, how="left", left_on = gene_id_col, right_on = "ensembl_gene")
    df.drop(columns=['ensembl_gene'], inplace=True)

    ### Fill NA 
    df.fillna(".", inplace=True)

    ### Save output    
    df.to_csv(os.path.join(outd, outf), sep=delim_out, index=False)


if __name__ == "__main__":
    args = parse_args()

    if args.outf == "NA":
        args.outf = "geneannot.{}".format(os.path.split(args.file)[-1])

    if args.outd == "NA":
        args.outd = os.getcwd()
    
    if args.delim_out == "NA":
        args.delim_out = args.delim_in

    main(file=args.file, 
        delim_in=map_delim(args.delim_in), 
        compression_in=args.compression_in, 
        gene_id_col=args.gene_id_col, 
        gene_symbol_only=args.gene_symbol_only,
        outf=args.outf, 
        outd=args.outd, 
        delim_out=map_delim(args.delim_out)
        )
