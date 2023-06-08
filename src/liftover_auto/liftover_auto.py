import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from packages import *
from util import logger, save_log, map_delim
from fnc import *


def parse_args():
    parser = argparse.ArgumentParser(description=":: Make a neat summarized results of FUMA and MAGMA ::")

    # Input
    parser.add_argument('--file', required=True,
                        help='Path to the input file being lifted over.')
    parser.add_argument('--in-compression', dest="in_compression", required=False, default="infer",
                        help="Specify compression type from the following ['zip', 'gzip', 'bz2', 'zstd', 'tar']. Default='infer'.")
    parser.add_argument('--delim', required=False, default="tab",
                    help="Delimiter used in the input file. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.")
    parser.add_argument('--infer-chr-pos', dest="infer_chr_pos", nargs='+', required=False,
                        help="Specify column name, data format, and separator to infer chr and pos from specified column.\
                            For example, column named 'variant' have variant name '2:179816990:C:T' where chromosome and position can be inferred as 2 and 179816990, respectively.\
                            Then, specify as follows: --infer-chr-pos variant CHR:POS :")
    parser.add_argument('--snp-col', dest="snp_col", required=False, default="SNP",
                    help="Name of the SNP column in the input file. Default = 'SNP'.")
    parser.add_argument('--chr-col', dest="chr_col", required=False, default="CHR",
                    help="Name of the chromosome column in the input file. Default = 'CHR'.")
    parser.add_argument('--pos-col', dest="pos_col", required=False, default="POS",
                    help="Name of the base position column in the input file. Default = 'POS'.")

    # New column names
    parser.add_argument('--chr-col-new', dest="chr_col_new", required=False, default="NA",
                    help="New column name for chromosome in the lifted file. Default = 'NA'.")
    parser.add_argument('--pos-col-new', dest="pos_col_new", required=False, default="NA",
                    help="New column name for base position in the lifted file. Default = 'NA'.")

    # Liftover options.
    parser.add_argument('--build-from', dest="build_from", required=True, type=int,
                    help="Genome build number initial. Choices = [18, 37, 38].")
    parser.add_argument('--build-to', dest="build_to", required=True, type=int,
                    help="Genome build number after performing liftover. Choices = [19, 37, 38].")
    
    # Merging lifted result options.
    parser.add_argument('--keep-initial-pos', dest="keep_initial_pos", action='store_true',
                help='Specify to save both previous and lifted Chr and Pos columns. Default = False.')
    parser.add_argument('--keep-unlifted', dest="keep_unlifted", action='store_true',
                help='Specify to retain unlifted SNPs with their base position value as -9. Default = False.')
    parser.add_argument('--keep-intermediate', dest="keep_intermediate", action='store_true',
                help='Specify to keep the intermediate files generated; <>.bed, <>.liftover.lifted, and <>.liftover.unlifted. Default = False.')
    parser.add_argument('--unlifted-snplist', dest="unlifted_snplist", action='store_true',
                help='Specify to save the SNP list that has been unlifed. <>.unlifted.snplist. Default = False.')

    # Output options.
    parser.add_argument('--outf', required=False, default="NA",
                        help="Specify the name of the output file. Default = 'lifted.<file>'.")
    parser.add_argument('--outd', required=False, default="NA",
                        help="Specify the path to output directory. Default = Current working directory.")
    parser.add_argument('--out-compression', dest="out_compression", required=False, default=None,
                        help="Specify compression type from the following ['zip', 'gzip', 'bz2', 'zstd', 'tar']. Default='infer'.")

    # Other optional.
    parser.add_argument('--verbose', action='store_true',
                help='Specify see the summary in the terminal. Default = False')

    args = parser.parse_args()
    return args


def main(file, in_compression, delim, snp_col, infer_chr_pos, chr_col, pos_col, 
        chr_col_new, pos_col_new,
        build_from, build_to, 
        keep_initial_pos, keep_unlifted, keep_intermediate, unlifted_snplist, 
        outf, outd, out_compression,
        verbose):
    
    logs_ = []
    log_ = "Performing liftover for:\n\t{}\n\tGenome build from GRCh{} to GRCh{}".format(file, build_from, build_to); logger(logs_, log_, verbose)

    # 1. Make bed file
    _, filename = os.path.split(file)
    
    ## << 수정하기 (1) 
    """
    script를 다시 돌리는 건 bed 파일이 잘못 만들어진게 있는데, 다시 안 만들면 계속 liftover 결과 잘못된걸 테니까. 
    """
    # bed_file_exists = os.path.exists(os.path.join(outd, filename+".liftover.bed"))
    ## >> 수정하기 (1)
    input_bed = make_bed(file, in_compression, delim, snp_col, infer_chr_pos, chr_col, pos_col, outd, bed_file_exists=False)

    # 2. Perform liftover
    if not os.path.exists(os.path.join(outd, filename+".lifted")):
        run_liftover(input_bed, build_from, build_to, outd)

    # 3. Merge the result
    logs = merge_lifted(file, delim, snp_col, chr_col, pos_col, chr_col_new, pos_col_new,
                unlifted_snplist, keep_initial_pos, keep_unlifted, keep_intermediate,
                outf, outd, out_compression, verbose); logs_ += logs

    # 4. Save the log
    save_log(logs_, os.path.join(outd, outf+".liftover.log"))


if __name__ == "__main__":
    args = parse_args()

    if args.outf == "NA":
        args.outf = os.path.split(args.file)[1]
    if args.outd == "NA":
        args.outd = os.getcwd()

    main(file=args.file,
        in_compression=args.in_compression,
        delim=map_delim(args.delim),
        snp_col=args.snp_col,
        infer_chr_pos=args.infer_chr_pos,
        chr_col=args.chr_col,
        pos_col=args.pos_col,
        chr_col_new=args.chr_col_new,
        pos_col_new=args.pos_col_new,
        build_from=args.build_from,
        build_to=args.build_to,
        keep_initial_pos=args.keep_initial_pos, 
        keep_unlifted=args.keep_unlifted,
        keep_intermediate=args.keep_intermediate, 
        unlifted_snplist=args.unlifted_snplist,
        outf=args.outf,
        outd=args.outd,
        out_compression=args.out_compression,
        verbose=args.verbose
        )
