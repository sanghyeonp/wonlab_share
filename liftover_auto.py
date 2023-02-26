from src.packages import *
from src.fnc import *
from src.util import logger, save_log, map_delim


def parse_args():
    parser = argparse.ArgumentParser(description=":: Make a neat summarized results of FUMA and MAGMA ::")

    # Input
    parser.add_argument('--file', required=True,
                        help='Path to the input file being lifted over.')
    parser.add_argument('--delim', required=False, default="tab",
                    help="Delimiter used in the input file. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.")
    parser.add_argument('--snp_col', required=False, default="SNP",
                    help="Name of the SNP column in the input file. Default = 'SNP'.")
    parser.add_argument('--chr_col', required=False, default="CHR",
                    help="Name of the chromosome column in the input file. Default = 'CHR'.")
    parser.add_argument('--pos_col', required=False, default="POS",
                    help="Name of the base position column in the input file. Default = 'POS'.")

    # Liftover options.
    parser.add_argument('--build_from', required=True, type=int,
                    help="Genome build number initial. Choices = [18, 37, 38].")
    parser.add_argument('--build_to', required=True, type=int,
                    help="Genome build number initial. Choices = [19, 37, 38].")
    
    # Merging lifted result options.
    parser.add_argument('--keep_all_col', action='store_true',
                help='Specify to save both previous and lifted Chr and Pos columns. Default = False')
    parser.add_argument('--keep_intermediate', action='store_true',
                help='Specify to keep the intermediate files generated; <>.bed, <>.liftover.lifted, and <>.liftover.unlifted. Default = False')
    parser.add_argument('--unlifted_snplist', action='store_true',
                help='Specify to save the SNP list that have been unlifed. Default = False')

    # Output options.
    parser.add_argument('--outf', required=False, default="NA",
                        help="Specify the name of the output file. Default = 'lifted.<file>'.")
    parser.add_argument('--outd', required=False, default="NA",
                        help="Specify the path to output directory. Default = Current working directory.")
    
    # Other optional.
    parser.add_argument('--verbose', action='store_true',
                help='Specify see the summary in the terminal. Default = False')

    args = parser.parse_args()
    return args


def main(file, delim, snp_col, chr_col, pos_col, 
        build_from, build_to, 
        keep_all_col, keep_intermediate, unlifted_snplist, 
        outf, outd,
        verbose):
    
    logs_ = []
    log_ = "Performing liftover for:\n\t{}\n\tGenome build from GRCh{} to GRCh{}".format(file, build_from, build_to); logger(logs_, log_, verbose)

    # 1. Make bed file
    input_bed = make_bed(file, delim, snp_col, chr_col, pos_col, outd)

    # 2. Perform liftover
    run_liftover(input_bed, build_from, build_to, outd)

    # 3. Merge the result
    logs = merge_lifted(file, delim, snp_col, chr_col, pos_col, 
                unlifted_snplist, keep_all_col, keep_intermediate,
                outf, outd, verbose); logs_ += logs

    # 4. Save the log
    save_log(logs_, os.path.join(outd, outf+".liftover.log"))


if __name__ == "__main__":
    args = parse_args()

    if args.outf == "NA":
        args.outf = os.path.split(args.file)[1]
    if args.outd == "NA":
        args.outd = os.getcwd()

    main(file=args.file,
        delim=map_delim(args.delim),
        snp_col=args.snp_col,
        chr_col=args.chr_col,
        pos_col=args.pos_col,
        build_from=args.build_from,
        build_to=args.build_to,
        keep_all_col=args.keep_all_col, 
        keep_intermediate=args.keep_intermediate, 
        unlifted_snplist=args.unlifted_snplist,
        outf=args.outf,
        outd=args.outd,
        verbose=args.verbose
        )
