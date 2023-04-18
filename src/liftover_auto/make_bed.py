import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from packages import *

def make_bed(file, file_compression, delim, snp_col, chr_col, pos_col, outd, bed_file_exists):
    _, filename = os.path.split(file)
    if bed_file_exists:
        return os.path.join(outd, filename+".liftover.bed")
    
    df_ = pd.read_csv(file, sep=delim, index_col=False, na_values=["NA", "NaN", " "], compression=file_compression)
    df_[chr_col] = df_[chr_col].astype(str)
    df = df_[[chr_col, pos_col, snp_col]]
    if 'chr' not in df.loc[0, chr_col]:
        df[chr_col] = df[chr_col].apply(lambda x: "chr{}".format(x))

    df['pos-1'] = df[pos_col].apply(lambda x: int(x) - 1)
    df = df[[chr_col, 'pos-1', pos_col, snp_col]]

    df.dropna(axis='index', how="any", subset=[snp_col], inplace=True)
    df.drop_duplicates(subset=snp_col, keep=False, inplace=True)

    df.to_csv(os.path.join(outd, filename+".liftover.bed"), sep="\t", header=False, index=False)

    return os.path.join(outd, filename+".liftover.bed")
