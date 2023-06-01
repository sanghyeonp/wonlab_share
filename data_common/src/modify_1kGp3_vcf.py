"""
Multi-allelic 인 경우.
AF for multi-allelic variants are reported for each allele independently, separated by ",".

"""

import pandas as pd
import os
import gzip
import dask.dataframe as dd
from tqdm import tqdm
tqdm.pandas(leave=False)
from timeit import default_timer as timer
from datetime import timedelta
import code
# code.interact(local=dict(globals(), **locals()))

dir_1kg = "/data/public/1kG/phase3/vcf"

### Get VCF path
# ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
file_path_list = [os.path.join(dir_1kg, "ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz".format(chr)) 
                for chr in range(1, 23)]
file_path_list += [os.path.join(dir_1kg, file) for file in ["ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz",
                                                            "ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz",
                                                            "ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz"
                                                            ]
                    ]

## < 지우기
file_path_list = [os.path.join(dir_1kg, file) for file in ["ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz",
                                                            "ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz",
                                                            "ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz"
                                                            ]
                    ]

## > 지우기

### Read VCF
def read_gz_vcf(vcf_path):
    start0 = timer()
    print("\n:: Reading VCF ::")
    with gzip.open(vcf_path, mode="rt") as f:
        file_content = f.readlines()
    print("Elapsed: {}".format(timedelta(seconds=timer() - start0)))
    
    start1 = timer()
    print(":: Modify VCF contents ::")
    file_content2 = [row.replace("\n", "").split(sep="\t")[:9] for row in tqdm(file_content, leave=False, desc="Input processing") if "##" not in row]

    del file_content

    print("Elapsed: {}".format(timedelta(seconds=timer() - start1)))
    df = pd.DataFrame(file_content2[1:], columns=file_content2[0])
    
    del file_content2
    
    # Drop sample genotype columns
    retain_col_list = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    df = df[retain_col_list]

    # Filter any rows with None
    df = df.loc[~df['POS'].isin([None])]
    df = df.loc[~df['ID'].isin([None])]
    df = df.loc[~df['REF'].isin([None])]
    df = df.loc[~df['ALT'].isin([None])]
    return df


def reformat_vcf(vcf_df):
    def INFO_col_dict(INFO, value):
        content = INFO.split(sep=";")
        content_map = dict()
        for ele in content:
            if "=" in ele:
                content_map[ele.split(sep="=")[0]] = ele.split(sep="=")[1]
            else:
                content_map[ele] = 'True'
        if value in content_map:
            return content_map[value]
        if 'AF' in value:
            return '-9'
        if value in ['IMPRECISE', 'EX_TARGET', 'MULTI_ALLELIC']:
            return 'False'
        return 'NA'
    
    # Extract information from INFO column
    vcf_df['AF_global'] = vcf_df['INFO'].progress_apply(lambda x: INFO_col_dict(x, 'AF'))
    vcf_df['AF_EAS'] = vcf_df['INFO'].progress_apply(lambda x: INFO_col_dict(x, 'EAS_AF'))
    vcf_df['AF_EUR'] = vcf_df['INFO'].progress_apply(lambda x: INFO_col_dict(x, 'EUR_AF'))
    vcf_df['AF_AMR'] = vcf_df['INFO'].progress_apply(lambda x: INFO_col_dict(x, 'AMR_AF'))
    vcf_df['AF_AFR'] = vcf_df['INFO'].progress_apply(lambda x: INFO_col_dict(x, 'AFR_AF'))
    vcf_df['AF_SAS'] = vcf_df['INFO'].progress_apply(lambda x: INFO_col_dict(x, 'SAS_AF'))
    vcf_df['VT'] = vcf_df['INFO'].progress_apply(lambda x: INFO_col_dict(x, 'VT'))
    vcf_df['IMPRECISE'] = vcf_df['INFO'].progress_apply(lambda x: INFO_col_dict(x, 'IMPRECISE'))
    vcf_df['EX_TARGET'] = vcf_df['INFO'].progress_apply(lambda x: INFO_col_dict(x, 'EX_TARGET'))
    vcf_df['MULTI_ALLELIC'] = vcf_df['INFO'].progress_apply(lambda x: INFO_col_dict(x, 'MULTI_ALLELIC'))

    return vcf_df



for vcf_path in tqdm(file_path_list, desc="VCF", leave=False):
    chr = [ele for ele in vcf_path.split(sep=".") if 'chr' in ele][0].replace('chr', '')
    
    print("\n##################### Chromosome {} #####################".format(chr))
    
    temp = read_gz_vcf(vcf_path)
    df = reformat_vcf(temp)
    del temp
    
    # retain_col = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 
    #             'AF_global', 'AF_EAS', 'AF_EUR', 'AF_AMR', 'AF_AFR', 'AF_SAS', 
    #             'IMPRECISE', 'EX_TARGET', 'MULTI_ALLELIC'
    #             ]

    # Make all values to be in string format
    df=df.astype(str)
    
    df.to_csv("/data1/sanghyeon/wonlab_contribute/combined/data_common/1kGp3_vcf/1kgp3_chr{}.tsv.gz".format(chr), 
            sep="\t", index=False, compression="gzip")
    del df
    
    ### 굳이 parquet으로 저장안해도 될 것 같음. (chr1을 tsv로 gzip 하면 약 120M)
    # Convert to dask dataframe and partition it
    # ddf = dd.from_pandas(df, npartitions=100)
    # del temp; del df
    
    # Make all values in dask dataframe into object
    # ddf = ddf.astype('object') # 이미 다 string으로 바꿔놔서 자동으로 됨.
    
    # ddf.to_parquet("/data1/sanghyeon/wonlab_contribute/combined/data_common/1kGp3_vcf/1kGp3_chr{}_parquet_partitioned".format(chr), 
    #                 engine="pyarrow", 
    #                 compression="snappy"
    #                 )
    # del ddf
    