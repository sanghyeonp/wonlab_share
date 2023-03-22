import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from packages import *
from util import map_delim, read_formatted_file
from shared_data import read_geneinfo


def parse_args():
    parser = argparse.ArgumentParser(description=":: Map NCBI Entrez gene ID to gene name ::")
    parser.add_argument('--file', required=True,
                        help='Path to the input file.')
    parser.add_argument('--gene_id_col', required=True, 
                        help='Column name of gene ID in the input file.')
    parser.add_argument('--delim', required=False, type=str, default="tab",
                        help='Delimiter of the input file. Default = "tab". Choices = ["tab", "comma", "whitespace", "formatted"].')
    parser.add_argument('--skip_rows', required=False, type=int, default=0,
                        help='Specify the number of first lines in the input file to skip. Default = 0.')

    parser.add_argument('--outf', required=False, type=str, default='NA',
                        help='Name of output file. Default = id2name.<Input file name>')
    parser.add_argument('--outd', required=False, type=str, default='NA',
                        help='Directory where the output will be saved. Default = Current working directory.')

    args = parser.parse_args()
    return args


def main(file, skip_rows, delim, gene_id_col, outf, outd):
    if delim == "formatted":
        df = read_formatted_file(file)
    else:
        df = pd.read_csv(file, sep=delim, skiprows=skip_rows)

    df[gene_id_col] = df[gene_id_col].astype('int64')

    df_map = read_geneinfo()

    df2 = df.merge(df_map, how="left", left_on=gene_id_col, right_on='GeneID')

    df2.drop(columns=['GeneID'], inplace=True)

    df.fillna({'GeneName':"NA"}, inplace=True)

    if delim == "formatted":
        delim = "\t"
    df2.to_csv(os.path.join(outd, outf), sep=delim, index=False)


if __name__ == "__main__":
    args = parse_args()

    if args.outf == "NA":
        args.outf = "id2name.{}".format(os.path.split(args.file)[-1])

    if args.outd == "NA":
        args.outd = os.getcwd()

    main(file=args.file, 
        skip_rows=args.skip_rows,
        delim=map_delim(args.delim),
        gene_id_col=args.gene_id_col,  
        outf=args.outf, 
        outd=args.outd
        )
    