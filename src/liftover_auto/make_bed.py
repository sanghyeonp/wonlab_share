import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from packages import *

def make_bed(file, file_compression, delim, snp_col, infer_chr_pos, chr_col, pos_col, outd, bed_file_exists):
    _, filename = os.path.split(file)
    if bed_file_exists:
        return os.path.join(outd, filename+".liftover.bed")
    
    df_ = pd.read_csv(file, sep=delim, index_col=False, na_values=["NA", "NaN", " "], compression=file_compression)

    ## Handle CHR:POS infer
    if infer_chr_pos:
        # 지금은 CHR:POS:REF/ALT:REF/ALT format이라고 생각하고 코딩한 것.
        # 하지만 CHR:POS 이런 식 일 수도 있음.
        # 이거 참고. https://www.terality.com/post/why-pandas-apply-method-is-slow
        infer_col = infer_chr_pos[0]
        data_format = infer_chr_pos[1]
        separator = infer_chr_pos[2]

        data_format = data_format.split(sep=separator)
        idx = {v:data_format.index(v) for v in ['CHR', 'POS']}

        df_[chr_col] = df_[infer_col].apply(lambda x: x.split(sep=separator)[idx['CHR']])
        df_[pos_col] = df_[infer_col].apply(lambda x: x.split(sep=separator)[idx['POS']])

    df_[chr_col] = df_[chr_col].astype(str)
    df = df_[[chr_col, pos_col, snp_col]]
    if 'chr' not in df.loc[0, chr_col]:
        df[chr_col] = df[chr_col].apply(lambda x: "chr{}".format(x))

    def fnc(pos):
        try:
            return int(pos) - 1
        except:
            return -9
    
    df['pos-1'] = df[pos_col].apply(lambda x: fnc(x))
    df = df[[chr_col, 'pos-1', pos_col, snp_col]]

    df.dropna(axis='index', how="any", subset=[snp_col], inplace=True)
    df.drop_duplicates(subset=snp_col, keep=False, inplace=True)

    df.to_csv(os.path.join(outd, filename+".liftover.bed"), sep="\t", header=False, index=False)

    return os.path.join(outd, filename+".liftover.bed")
