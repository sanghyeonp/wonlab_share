"""
이거 참고해서 해보기: https://www.coiled.io/blog/dask-dataframe-merge-join

MAF는 1kG 데이터에서 가져오는 걸로 준비해보기.
"""

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
    print("Elapsed: {}".format(timedelta(seconds=timer() - start0)))

    ### Read mapping data
    start1 = timer()
    print(":: Reading annotation file ::")
    print(NCBI_VARIANT_ANNOT_DIR[genome_build])
    df_annot = dd.read_parquet(path=NCBI_VARIANT_ANNOT_DIR[genome_build])
    
    if chromosome_filter == "autosome":
        df_annot = df_annot[~df_annot['CHR'].isin(['X', 'Y', 'MT'])]
    print("Elapsed: {}".format(timedelta(seconds=timer() - start1)))

    ### Merge
    print(":: Mapping ::")
    df = df_annot.merge(
            df_, 
            how="right", 
            left_on="SNP",
            right_on=rsid_col
        )

    df2 = df.compute()

    ### Handle missing value
    # CHR, POS: -9 | ALT, REF: .
    start2 = timer()
    print(":: Handling missing values ::")
    df2['CHR'] = df2['CHR'].fillna('-9')
    df2['POS'] = df2['POS'].fillna(-9)
    df2['CHR'] = df2['CHR'].astype(str)
    df2['POS'] = df2['POS'].astype(int)
    df2.fillna(".", inplace=True)
    print("Elapsed: {}".format(timedelta(seconds=timer() - start2)))

    ### Sort values by chromosome and base position
    start3 = timer()
    print(":: Sort by chromosome and base position ::")
    df3 = df2[~df2['CHR'].isin(['X', 'Y', 'MT'])]
    df3['CHR'] = df3['CHR'].astype(int)
    df3.sort_values(by=['CHR', 'POS'], ascending=[True, True], inplace=True)

    if chromosome_filter == "all":
        df_other = df2[df2['CHR'].isin(['X', 'Y', 'MT'])]
        custom_dict = {'X':1, 'Y':2, 'MT':3}  
        df_other['code'] = df_other['CHR'].map(custom_dict)
        df_other.sort_values(by=['code', 'POS'], ascending=[True, True], inplace=True)
        df_other.drop(columns=['code'], inplace=True)
    
        df3 = pd.concat([df3, df_other])
    print("Elapsed: {}".format(timedelta(seconds=timer() - start3)))

    ### Save the result
    start4 = timer()
    print(":: Saving the annotated result ::")
    if outf == "NA":
        outf = "abc." + file
    if outd == "NA":
        outd = "."
    if compression_out == "NA":
        compression_out = None
    print("Path: {}".format(os.path.join(outd, outf)))

    df3.to_csv(os.path.join(outd, outf), sep=delim_out, 
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




# ############################################################################

# # import pandas as pd
# import os
# import pickle
# from tqdm import tqdm
# # tqdm.pandas()
# import code
# # code.interact(local=dict(globals(), **locals()))
# import pyarrow as pa
# import pyarrow.parquet as pq
# import dask.dataframe as dd
# from timeit import default_timer as timer
# from datetime import timedelta


# ########################################
# # cis-eQTL summary statistics all
# ########################################
# # dir = "/data1/sanghyeon/Projects/TAS1R3_SNU/obesity/data/eQTL/INTERVAL/cis_eQTLs_INTERVAL"

# # files = [os.path.join(dir, "tensorqtl_cis_MAF0.005_cisNominal_chr{}.csv".format(chr)) for chr in range(1, 23)]

# # df_ = None
# # for file in files:
# #     print(file)
# #     temp = pd.read_csv(file, sep=",", index_col=False)
# #     df_ = pd.concat([df_, temp])

# # df_.to_csv("INTERVAL_eQTL_summary_statistics_temp.csv", sep=",", index=False)

# ########################################################################################################################

# ########################################
# # Read cis-eQTL summary statistics all
# ########################################
# df_ = dd.read_csv("INTERVAL_eQTL_summary_statistics_temp.csv.gz", blocksize=None, compression='gzip')
# ########################################################################################################################

# ########################################
# # Read annotation file
# ########################################
# n_thread = 5
# dd.config.set(scheduler='processes', num_workers = n_thread)  
# try:
#     start = timer()
#     parquet_dir = "/data1/sanghyeon/wonlab_contribute/combined/data_common/00_All_b37_vcf_CHR"
#     parquet_files = os.listdir(parquet_dir)
#     for idx, parquet_file in enumerate(tqdm(parquet_files, desc="Chromosome", leave=False)):
#         # df_parquet = pq.read_table(source=os.path.join(parquet_dir, parquet_file)).to_pandas()
#         df_parquet = dd.read_parquet(os.path.join(parquet_dir, parquet_file), engine="pyarrow")

#         df_parquet = df_parquet.astype({'SNP': 'object',
#                             'CHR':'object',
#                             'POS':'int64',
#                             'REF':'object',
#                             'ALT':'object'})
#         if idx == 0:
#             df_ = dd.merge(df_, df_parquet, how="left", left_on='variant_id', right_on='SNP')
#             # df_.merge(df_parquet, , how='left')
#             df_ = df_.drop(columns=['SNP'])
#         else:
#             df_ = df_.merge(df_parquet, left_on='variant_id', right_on='SNP', how='left', suffixes=['', '_'])
#             df_ = df_.drop(columns=['SNP'])
            
#             df_['CHR'] = df_['CHR'].fillna(df_.pop('CHR_'))
#             df_['POS'] = df_['POS'].fillna(df_.pop('POS_'))
#             df_['REF'] = df_['REF'].fillna(df_.pop('REF_'))
#             df_['ALT'] = df_['ALT'].fillna(df_.pop('ALT_'))
# except:
#     print("SOMETHING'S WRONG!!!!")
#     code.interact(local=dict(globals(), **locals()))

# print("Merging elapse: ")
# print(timedelta(seconds=timer()-start))

# values = {"CHR": "-9", "POS": -9, "REF": ".", "ALT": "."}  
# df_ = df_.fillna(value=values)

# start2 = timer()
# df_.compute().to_csv("./INTERVAL_eQTL_summary_statistics_b37.csv.gz", index=False, compression='gzip')
# print("Saving elapsed:")
# print(timedelta(seconds=timer() - start2))
# code.interact(local=dict(globals(), **locals()))
