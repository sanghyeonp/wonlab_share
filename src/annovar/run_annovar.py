import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *
from make_annovar_input import make_annovar_input
from annovar import annovar


def parse_args():
    parser = argparse.ArgumentParser(description=":: Make ANNOVAR input file from user-specified input file ::")
    parser.add_argument('--file', required=True,
                        help='Path to the input file.')

    parser.add_argument('--file_compression', required=False, default="infer",
                        help="Specify compression type from the following ['zip', 'gzip', 'bz2', 'zstd', 'tar']. Default='infer'.")
    parser.add_argument('--delim_in', required=False, default="tab",
                    help="Delimiter used in the input file. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.")
    
    parser.add_argument('--infer_chr_pos_ref_alt', nargs='+', required=False,
                        help="Specify column name, data format, and separator to infer chr and pos from specified column.\
                            For example, column named 'variant' have variant name '2:179816990:C:T' where chromosome and position can be inferred as 2 and 179816990, respectively.\
                            Then, specify as follows: --infer_chr_pos_ref_alt variant CHR:POS:REF:ALT :")
    parser.add_argument('--chr_col', required=False, default="CHR",
                    help="Name of the chromosome column in the input file. Default = 'CHR'.")
    parser.add_argument('--pos_col', required=False, default="POS",
                    help="Name of the base position column in the input file. Default = 'POS'.")
    parser.add_argument('--ref_col', required=False, default="REF",
                    help="Name of the reference allele column in the input file. Default = 'REF'.")
    parser.add_argument('--alt_col', required=False, default="ALT",
                    help="Name of the alternative allele column in the input file. Default = 'ALT'.")
    
    args = parser.parse_args()
    return args


def main(file, delim_in, file_compression,
        infer_chr_pos_ref_alt, chr_col, pos_col, ref_col, alt_col):
    
    # Make ANNOVAR input file
    annov_input1, annov_input2 = make_annovar_input(file, delim_in, file_compression, infer_chr_pos_ref_alt, chr_col, pos_col, ref_col, alt_col)

    # Run ANNOVAR
    annovar(annov_input1, annov_input2)

    # Map ANNOVAR result
    

if __name__ == "__main__":
    args = parse_args()

    main(file=args.file, 
        delim_in=map_delim(args.delim_in), 
        file_compression=args.file_compression,
        infer_chr_pos_ref_alt=args.infer_chr_pos_ref_alt, 
        chr_col=args.chr_col, 
        pos_col=args.pos_col, 
        ref_col=args.ref_col, 
        alt_col=args.alt_col
        )

