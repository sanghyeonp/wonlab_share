import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *
from shared_data import rsID_MAP_dbSNP


def parse_args():
    parser = argparse.ArgumentParser(description=":: Map rsID using chromosome and base position ::")
    parser.add_argument('--file', required=True,
                        help='Path to the input file.')
    parser.add_argument('--build', required=True, type=int,
                        help="Genome build number (GRCh). Choices = [18, 37, 38].")
    
    parser.add_argument('--file_compression', required=False, default="infer",
                        help="Specify compression type from the following ['zip', 'gzip', 'bz2', 'zstd', 'tar']. Default='infer'.")
    parser.add_argument('--delim_in', required=False, default="tab",
                    help="Delimiter used in the input file. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.")
    
    parser.add_argument('--infer_chr_pos', nargs='+', required=False,
                        help="Specify column name, data format, and separator to infer chr and pos from specified column.\
                            For example, column named 'variant' have variant name '2:179816990:C:T' where chromosome and position can be inferred as 2 and 179816990, respectively.\
                            Then, specify as follows: --infer_chr_pos variant CHR:POS:REF:ALT :")
    parser.add_argument('--chr_col', required=False, default="CHR",
                    help="Name of the chromosome column in the input file. Default = 'CHR'.")
    parser.add_argument('--pos_col', required=False, default="POS",
                    help="Name of the base position column in the input file. Default = 'POS'.")
    parser.add_argument('--new_snp_col', required=False, default="SNP_ref",
                    help="Name of the new SNP column mapped from the reference mapping file. Default = 'SNP_ref'.")

    # 
    parser.add_argument('--keep_unmapped', action='store_true',
                        help='Specify to keep unmapped SNPs in the output file.')
    parser.add_argument('--keep_refmap_position', action='store_true',
                        help='Specify to keep the genome position of the reference mapping file (columns: CHR_ref, POS_ref) in the output file rather than dropping them. \
                                Keeping theses columns allow comparisons of input CHR and POS to the reference mapping file CHR and POS. Default = False.')
    parser.add_argument('--keep_refmap_allele', action='store_true',
                        help='Specify to keep the reference and alternative allele from the reference mapping file (columns: REF_ref, ALT_ref) in the output file rather than dropping them. \
                                Keeping theses columns allow comparisons of input REF and ALT to the reference mapping file REF and ALT. Default = False.')
    parser.add_argument('--keep_col', nargs='+', required=False,
                        help='Keep only the specified columns in the output file. Multiple column names can be specified as following: --keep_col COL1 COL2 COL3')

    # Output options.
    parser.add_argument('--outf', required=False, default="NA",
                        help="Specify the name of the output file. Default = 'rsidmapped.<file>'.")
    parser.add_argument('--outd', required=False, default="NA",
                        help="Specify the path to output directory. Default = Current working directory.")
    parser.add_argument('--delim_out', required=False, default="NA",
                        help="Delimiter for the output file. Choices = ['NA', 'tab', 'comma', 'whitespace']. If 'NA', identical delimiter as input delimiter will be used. Default = 'NA'.")
    parser.add_argument('--out_compression', required=False, default="NA",
                        help="Specify compression type from the following ['NA', 'zip', 'gzip', 'bz2', 'zstd', 'tar']. Default='NA'.")

    args = parser.parse_args()
    return args


def main(file, file_compression, delim_in, 
        mapping_ref,
        infer_chr_pos, chr_col, pos_col, new_snp_col,
        outf, outd, delim_out, out_compression, 
        keep_unmapped, keep_col,
        keep_refmap_position, keep_refmap_allele):
    
    # Read input file
    print("## Reading the input file...")
    df_ = pd.read_csv(file, sep=delim_in, index_col=False, compression=file_compression, low_memory=False)
    print("Number of SNPs in the input file: {:,}".format(len(df_)))

    # Infer chromosome and position
    if infer_chr_pos:
        infer_col = infer_chr_pos[0]
        data_format = infer_chr_pos[1]
        separator = infer_chr_pos[2]

        print("## Inferring CHR and POS from the column '{}' with data structure '{}'...".format(infer_col, data_format))

        data_format = data_format.split(sep=separator)
        idx = {v:data_format.index(v) for v in ['CHR', 'POS']}

        df_[[chr_col, pos_col]] = df_.apply(lambda row: [row[infer_col].split(sep=separator)[idx['CHR']], 
                                                    row[infer_col].split(sep=separator)[idx['POS']]
                                                    ], axis=1, result_type='expand')


    df_[chr_col] = df_[chr_col].astype(str)
    df_[pos_col] = df_[pos_col].astype(str)

    # Check CHR
    chr_check = 'chr' in df_.loc[0, 'CHR']
    if chr_check:
        df_['CHR'] = df_['CHR'].apply(lambda x: x.replace('chr', ''))

    # Read mapping reference file
    print("## Reading the reference mapping file...")
    df_ref = pd.read_csv(mapping_ref, sep="\t", index_col=False, compression='gzip', low_memory=False)
    df_ref.rename(columns={'SNP':'SNP_ref',
                        'CHR':'CHR_ref',
                        'POS':'POS_ref',
                        'REF':'REF_ref',
                        'ALT':'ALT_ref'}, inplace=True)
    print("Number of SNPs in the reference mapping file: {:,}".format(len(df_ref)))

    df_ref['CHR_ref'] = df_ref['CHR_ref'].astype(str)
    df_ref['POS_ref'] = df_ref['POS_ref'].astype(str)

    # Map rsID
    print("## Mapping the rsID...")
    df = df_.merge(df_ref, how="left", left_on=[chr_col, pos_col], right_on=['CHR_ref', 'POS_ref'])

    print("Number of SNPs mapped to rsID / not mapped to rsID: {:,} / {:,}".format(len(df[~df['SNP_ref'].isna()]), len(df[df['SNP_ref'].isna()])))

    if not keep_unmapped:
        print("## Dropping SNPs with unmapped rsID...")
        df = df[~df['SNP_ref'].isna()]
    # If unmapped SNPs are kept, let it have the following rsID
    df['SNP_ref'] = df.apply(lambda row: "{}:{}:{}:{}".format(row[chr_col], row[pos_col]) if row['SNP_ref'].isna() else row['SNP_ref'], axis=1)

    if not keep_refmap_position:
        print("## Dropping CHR and POS from the reference mapping file...")
        df.drop(columns=['CHR_ref', 'POS_ref'], inplace=True)
    if not keep_refmap_allele:
        print("## Dropping REF and ALT from the reference mapping file...")
        df.drop(columns=['REF_ref', 'ALT_ref'], inplace=True)

    df.rename(columns={'SNP_ref': new_snp_col}, inplace=True)

    if keep_col:
        assert all(col in df.columns.tolist() for col in keep_col), "Specified column(s) is/are not present in the dataframe. Columns present in the dataframe: {}".join(', '.join(df.columns.tolist()))
        print("## Leaving only the specified columns: {}".format(keep_col))
        df = df[keep_col]

    # Save output
    if out_compression == "NA":
        out_compression = None
    df.to_csv(os.path.join(outd, outf), sep=delim_out, index=False, compression=out_compression)



if __name__ == "__main__":
    args = parse_args()

    mapping_ref = rsID_MAP_dbSNP[args.build]

    if args.outf == "NA":
        args.outf = "rsidmapped.{}".format(os.path.split(args.file)[-1])

    if args.outd == "NA":
        args.outd = os.getcwd()

    main(file=args.file, 
        file_compression=args.file_compression, 
        delim_in=map_delim(args.delim_in), 
        mapping_ref=mapping_ref,
        infer_chr_pos=args.infer_chr_pos,
        chr_col=args.chr_col, 
        pos_col=args.pos_col, 
        new_snp_col=args.new_snp_col,
        outf=args.outf, 
        outd=args.outd, 
        out_compression=args.out_compression,
        delim_out=map_delim(args.delim_out), 
        keep_unmapped=args.keep_unmapped, 
        keep_col=args.keep_col,
        keep_refmap_position=args.keep_refmap_position, 
        keep_refmap_allele=args.keep_refmap_allele
        )
