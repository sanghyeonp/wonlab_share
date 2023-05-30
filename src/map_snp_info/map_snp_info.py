import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *
from shared_data import NCBI_VARIANT_ANNOT_DIR

import_or_install('dask')
import dask
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
pbar = ProgressBar()
pbar.register()
from timeit import default_timer as timer
from datetime import timedelta


def parse_args():
    parser = argparse.ArgumentParser(description=":: Map CHR, POS, REF, and ALT to rsID ::")
    
    # Required arguments
    parser.add_argument('--file', required=True,
                        help='Path to the input file.')
    parser.add_argument('--rsid-col', dest="rsid_col", required=True,
                        help="Name of the rsID column in the input file.")
    parser.add_argument('--build', required=True, type=int,
                        help="Genome build to map CHR, POS, REF, and ALT. Choices = [37, 38]")
    
    # Optional arguments
    parser.add_argument('--delim-in', dest="delim_in", required=False, default="tab",
                        help="Delimiter used in the input file. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.")
    parser.add_argument('--compression-in', dest="compression_in", required=False, default="NA",
                        help="Specify compression type from the following ['zip', 'gzip', 'bz2', 'zstd', 'tar']. Default='NA'.")
    parser.add_argument('--chromosome-filter', dest="chromosome_filter", required=False, default="autosome",
                        help="Specify which chromosomes to include. Default = 'autosome'. If 'all' is specified, then chromosome X, Y, MT are included. \
                        Choices = ['autosome', 'all']")

    parser.add_argument('--outf', required=False, default="NA",
                        help="Specify the name of the output file. Default = 'abc.<file>'.")
    parser.add_argument('--outd', required=False, default="NA",
                        help="Specify the path to output directory. Default = Current working directory.")
    parser.add_argument('--delim-out', dest="delim_out", required=False, default="NA",
                        help="Delimiter for the output file. Choices = ['NA', 'tab', 'comma', 'whitespace']. If 'NA', identical delimiter as input delimiter will be used. Default = 'NA'.")
    parser.add_argument('--compression-out', dest="compression_out", required=False, default="NA",
                        help="Specify compression type from the following ['NA', 'zip', 'gzip', 'bz2', 'zstd', 'tar']. Default='NA'.")


    args = parser.parse_args()
    return args



def main(file, delim_in, compression_in, 
        rsid_col, genome_build, chromosome_filter,
        outf, outd, delim_out, compression_out
        ):
    
    ### Read input file
    start0 = timer()
    print(":: Reading input file ::")
    if compression_in == "NA":
        compression_in = None
    df_ = dd.read_csv(file, blocksize=None, sep=delim_in, compression=compression_in)

    snp_list = df_[rsid_col].compute().tolist()

    print("Elapsed: {}".format(timedelta(seconds=timer() - start0)))

    ### Read annotation data
    start1 = timer()
    print(":: Reading annotation file ::")
    print(NCBI_VARIANT_ANNOT_DIR[genome_build])
    df_annot = dd.read_parquet(path=NCBI_VARIANT_ANNOT_DIR[genome_build])
    
    if chromosome_filter == "autosome":
        df_annot = df_annot[~df_annot['CHR'].isin(['X', 'Y', 'MT'])]
    print("Elapsed: {}".format(timedelta(seconds=timer() - start1)))

    # Rename the annotation dataframe columns
    col_rename = {"CHR":"CHR_b{}".format(genome_build), 
                    "POS":"POS_b{}".format(genome_build),
                    "REF":"REF_b{}".format(genome_build),
                    "ALT":"ALT_b{}".format(genome_build)
                    }
    df_annot = df_annot.rename(columns=col_rename)

    ### Merge
    print(":: Mapping ::")
    
    df_annot_filter = df_annot[df_annot['SNP'].isin(snp_list)]

    df = df_.merge(df_annot_filter, how="left", left_on=rsid_col, right_on="SNP")

    df2 = df.compute()

    ### Handle missing value
    # CHR, POS: -9 | ALT, REF: .
    start2 = timer()
    print(":: Handling missing values ::")
    df2['CHR_b{}'.format(genome_build)] = df2['CHR_b{}'.format(genome_build)].fillna('-9')
    df2['POS_b{}'.format(genome_build)] = df2['POS_b{}'.format(genome_build)].fillna(-9)
    df2['CHR_b{}'.format(genome_build)] = df2['CHR_b{}'.format(genome_build)].astype(str)
    df2['POS_b{}'.format(genome_build)] = df2['POS_b{}'.format(genome_build)].astype(int)
    df2.fillna(".", inplace=True)
    print("Elapsed: {}".format(timedelta(seconds=timer() - start2)))

    ### Sort values by chromosome and base position
    start3 = timer()
    print(":: Sort by chromosome and base position ::")
    df_missing = df2[df2['CHR_b{}'.format(genome_build)] == '-9']
    df3 = df2[df2['CHR_b{}'.format(genome_build)] != '-9']
    
    df_automsome = df3[~df3['CHR_b{}'.format(genome_build)].isin(['X', 'Y', 'MT'])]
    df_automsome['CHR_b{}'.format(genome_build)] = df_automsome['CHR_b{}'.format(genome_build)].astype(int)
    df_automsome.sort_values(by=['CHR_b{}'.format(genome_build), 'POS_b{}'.format(genome_build)], ascending=[True, True], inplace=True)

    
    if chromosome_filter == "all":
        df_other = df3[df3['CHR_b{}'.format(genome_build)].isin(['X', 'Y', 'MT'])]
        custom_dict = {'X':1, 'Y':2, 'MT':3}  
        df_other['code'] = df_other['CHR_b{}'.format(genome_build)].map(custom_dict)
        df_other.sort_values(by=['code', 'POS_b{}'.format(genome_build)], ascending=[True, True], inplace=True)
        df_other.drop(columns=['code'], inplace=True)
    
        df4 = pd.concat([df_automsome, df_other])
    else:
        df4 = df_automsome.copy()
        del df_automsome
    df4 = pd.concat([df4, df_missing])

    print("Elapsed: {}".format(timedelta(seconds=timer() - start3)))

    ### Save the result
    start4 = timer()
    print(":: Saving the annotated result ::")
    if outf == "NA":
        
        outf = "abc." + os.path.split(file)[-1]
    if outd == "NA":
        outd = "."
    if compression_out == "NA":
        compression_out = None
    print("Path: {}".format(os.path.join(outd, outf)))

    df4.to_csv(os.path.join(outd, outf), sep=delim_out, 
                index=False, compression=compression_out)
    print("Elapsed: {}".format(timedelta(seconds=timer() - start4)))


if __name__ == "__main__":
    args = parse_args()

    dask.config.set(scheduler="threads", num_workers=1, threads_per_worker=2)

    if args.delim_out == "NA":
        args.delim_out = args.delim_in

    main(file=args.file,
        delim_in=map_delim(args.delim_in),
        compression_in=args.compression_in,
        rsid_col=args.rsid_col,
        genome_build=args.build,
        chromosome_filter=args.chromosome_filter,
        outf=args.outf,
        outd=args.outd,
        delim_out=map_delim(args.delim_out),
        compression_out=args.compression_out
        )
