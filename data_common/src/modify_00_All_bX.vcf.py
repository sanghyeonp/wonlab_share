"""
참고: https://www.coiled.io/blog/writing-parquet-files-with-dask-using-to-parquet
"""

import gzip
import pandas as pd
import code
import dask.dataframe as dd
# code.interact(local=dict(globals(), **locals()))

def vcf_to_parquet(vcf, build):
    try:
        print("[0]")
        with gzip.open(vcf, mode="rt") as f:
            file_content = f.read()

        file_content = file_content.split(sep="\n")

        file_content = [row.split(sep="\t") for row in file_content if "##" not in row]

        print("[1]")
        df = pd.DataFrame(file_content[1:], columns=file_content[0])
        
        del file_content

        print("[2]")
        df.dropna(axis=0, how="any", inplace=True)
        
        print("[3]")
        df = df[['ID', '#CHROM', 'POS', 'REF', 'ALT']]
        df.columns = ['SNP', 'CHR', 'POS', 'REF', 'ALT']
        
        print(df.dtypes)

        print("[4]")
        df['SNP'] = df['SNP'].astype("string")
        df["CHR"] = df["CHR"].astype("string")
        df["POS"] = df["POS"].astype("int64")
        df["REF"] = df["REF"].astype("string")
        df["ALT"] = df["ALT"].astype("string")
        
        print(df.dtypes)

        print("[5]")
        ddf = dd.from_pandas(df, npartitions=300)
        
        print("[6]")
        
        ddf = ddf.astype({'SNP': 'object',
                            'CHR':'object',
                            'POS':'int64',
                            'REF':'object',
                            'ALT':'object'})
        ddf.to_parquet("../00_All_b{}_parquet_partitioned".format(build), 
                    engine="pyarrow", 
                    compression="snappy"
                    )
    except:
        code.interact(local=dict(globals(), **locals()))


if __name__ == "__main__":
    vcf = "../00_All_b37.vcf.gz"
    build = 37
    vcf_to_parquet(vcf, build)

    vcf = "../00_All_b38.vcf.gz"
    build = 38
    vcf_to_parquet(vcf, build)
