import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *
from compute_ld_single import bfile_selection, compute_ld_single

from itertools import combinations
from timeit import default_timer as timer
from datetime import timedelta
tqdm.pandas(leave=False, bar_format='{l_bar}{bar:30}{r_bar}{bar:-30b}')


def parse_args():
    parser = argparse.ArgumentParser(description=":: Compute LD between two SNPs ::")
    
    # Required arguments
    parser.add_argument('--file-in', dest="file_in", required=True,
                        help='File containing two columns for the pairs of SNP to compute LD.')
    parser.add_argument('--delim-in', dest="delim_in", required=False, default="tab",
                    help="Delimiter used in the input file. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.")
    
    parser.add_argument('--all-possible-pairs', dest="all_possible_pairs", action='store_true',
                    help='Specify to compute LD for all possible pairs from SNP1 and SNP2 list. Default = False.')

    parser.add_argument('--threads', required=False, type=int, default=1,
                    help='Number of threads to use.')
    parser.add_argument('--variant-identifier', dest="variant_identifier", required=True,
                        help='Type of variant identifier. Choices = ["rsid", "chr:pos"]')
    parser.add_argument('--reference', required=True,
                        help='Type of reference. Choices = ["1kg", "ukb"]')
    parser.add_argument('--ancestry', required=True,
                        help='Ancestry of the input file. \
                            Choices from 1kg = ["European", "East-asian"] \
                            Choices from ukb = ["European"]')

    parser.add_argument('--outf', required=False, default="NA",
                        help="Specify the name of the output file. Default = 'r2.<file>'.")
    parser.add_argument('--outd', required=False, default="NA",
                        help="Specify the path to output directory. Default = Current working directory.")
    parser.add_argument('--delim-out', dest="delim_out", required=False, default="NA",
                        help="Delimiter for the output file. Choices = ['NA', 'tab', 'comma', 'whitespace']. If 'NA', identical delimiter as input delimiter will be used. Default = 'NA'.")
    args = parser.parse_args()
    return args


def compute_ld_parallel(bfile, snp1_list, snp2_list, fnc, threads):
    m = Manager()
    computed_r2 = m.list()

    pool = Pool(processes=threads)
    inputs = [(snp1, snp2_list[idx], bfile, computed_r2) for idx, snp1 in enumerate(snp1_list)]

    pool.starmap(fnc, tqdm(inputs, total=len(inputs), leave=False, desc="R2 computation", 
                        bar_format='{l_bar}{bar:30}{r_bar}{bar:-30b}'))
    pool.close()
    pool.join()

    run_bash("rm plink.log plink.nosex")

    return list(computed_r2)


def main(file_in, delim_in, all_possible_pairs,
        reference, variant_identifier, ancestry,
        threads,
        outf, outd, delim_out):
    
    ### Read input file
    df_ = pd.read_csv(file_in, sep=delim_in, index_col=False)
    col_snp1, col_snp2 = list(df_.columns)
    snp1_list, snp2_list = df_[col_snp1].tolist(), df_[col_snp2].tolist()

    ### All pair-wise specified
    if all_possible_pairs:
        snp_list = snp1_list + snp2_list
        unique_combinations = list(combinations(snp_list, 2))
        snp1_list, snp2_list = [], []
        for comb in unique_combinations:
            snp1_list.append(comb[0])
            snp2_list.append(comb[1])

        print("Number of unique combinations: {:,}".format(len(unique_combinations)))

    ### Specify bfile
    bfile = bfile_selection(reference, variant_identifier, ancestry)

    ### Compute LD by parallel computing
    computed_r2 = compute_ld_parallel(bfile = bfile,
                                    snp1_list = snp1_list,
                                    snp2_list = snp2_list,
                                    fnc = compute_ld_single,
                                    threads = threads
                                    )

    ### Make resultant dataframe
    df = pd.DataFrame(computed_r2, columns=[col_snp1, col_snp2, "R2"])

    if all_possible_pairs:
        # Make 
        df2 = df.copy()
        for _, row in df2.iterrows():
            if row[col_snp1] != row[col_snp2]:
                df = pd.concat([df, pd.DataFrame(list(zip([row[col_snp2]], [row[col_snp1]], [row['R2']])),
                                                    columns =[col_snp1, col_snp2, 'R2'])])
            df = pd.concat([df, pd.DataFrame(list(zip([row[col_snp1]], [row[col_snp1]], [1])), 
                                                columns =[col_snp1, col_snp2, 'R2'])])

        df = df.pivot_table(index=col_snp1, columns=col_snp2, values='R2', fill_value='.')

    ### Save the result
    if outf == "NA" and not all_possible_pairs:
        outf = "r2." + os.path.split(file_in)[-1]
    elif outf == "NA" and all_possible_pairs:
        outf = "r2_allpair." + os.path.split(file_in)[-1]
    if outd == "NA":
        outd = "."

    if all_possible_pairs:
        df.to_csv(os.path.join(outd, outf), sep=delim_out, index=True)
    else:
        df.to_csv(os.path.join(outd, outf), sep=delim_out, index=False)


if __name__ == "__main__":
    args = parse_args()

    if args.delim_out == "NA":
        args.delim_out = args.delim_in

    main(file_in=args.file_in, 
        delim_in=map_delim(args.delim_in), 
        all_possible_pairs=args.all_possible_pairs,
        reference=args.reference, 
        variant_identifier=args.variant_identifier, 
        ancestry=args.ancestry,
        threads=args.threads,
        outf=args.outf, 
        outd=args.outd, 
        delim_out=map_delim(args.delim_out)
        )
