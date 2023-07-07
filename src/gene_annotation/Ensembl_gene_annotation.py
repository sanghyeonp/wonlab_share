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
    
    parser.add_argument('--single-gene', dest="single_gene", action='store_true',
                    help='Specify get annotation for a single gene. Default = False.')
    parser.add_argument('--ensembl-id', dest="ensembl_id", required=False,
                        help='ENSEMBL ID for a single gene. Default = False.')

    parser.add_argument('--file', required=False, default="NA",
                        help='Path to the input file.')
    parser.add_argument('--gene-id-col', dest="gene_id_col", required=False, default="NA",
                        help='Column name of Ensembl gene ID in the input file.')
    parser.add_argument('--genome-build', dest='genome_build', required=False, default=37,
                        help='Specify genome build. Choices = [37, 38]. Default=37.')
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


def main(single_gene, ensembl_id, file, delim_in, compression_in, gene_id_col, genome_build, gene_symbol_only, outf, outd, delim_out):
    if file != "NA" and not single_gene:
        df_ = pd.read_csv(file, sep=delim_in, index_col=False, compression=compression_in, low_memory=False)

        ### Drop version in gene ID if present
        df_[gene_id_col] = df_[gene_id_col].astype(str)
        df_['gene_id_new'] = df_[gene_id_col].apply(lambda x: x.split(sep=".")[0] if "ENSG" in x and x != "nan" else x)

    ### Read gene annotation table
    df_gene_info = pd.read_csv(ENSEMBL_GENE_INFO[genome_build], sep="\t", index_col=False, compression="gzip")
    if gene_symbol_only:
        df_gene_info = df_gene_info[['ensembl_gene', 'gene_symbol']]

    ### Mapping
    if file == "NA" and single_gene:
        print(df_gene_info[df_gene_info['ensembl_gene'] == ensembl_id])
        return None
    
    df = df_.merge(df_gene_info, how="left", left_on = "gene_id_new", right_on = "ensembl_gene")
    df.drop(columns=['ensembl_gene', 'gene_id_new'], inplace=True)

    ### TSS
    if not gene_symbol_only:
        df['TSS'] = df['TSS'].astype('Int64')

    ### Fill NA
    if not gene_symbol_only:
        df[['gene_symbol', 'chr_gene', 'strand']] = df[['gene_symbol', 'chr_gene', 'strand']].fillna(value=".")
        df[['TSS']] = df[['TSS']].fillna(value=-9)
    else:
        df[['gene_symbol']] = df[['gene_symbol']].fillna(value=".")
        
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

    main(single_gene=args.single_gene,
        ensembl_id=args.ensembl_id,
        file=args.file, 
        delim_in=map_delim(args.delim_in), 
        compression_in=args.compression_in, 
        gene_id_col=args.gene_id_col, 
        genome_build=args.genome_build,
        gene_symbol_only=args.gene_symbol_only,
        outf=args.outf, 
        outd=args.outd, 
        delim_out=map_delim(args.delim_out)
        )
