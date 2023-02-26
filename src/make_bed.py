from src.packages import *

def make_bed(file, delim, snp_col, chr_col, pos_col, outd):
    df_ = pd.read_csv(file, delim, index_col=False, na_values=["NA", "NaN", " "])

    df = df_[[chr_col, pos_col, snp_col]]
    df[chr_col] = df[chr_col].apply(lambda x: "chr{}".format(x))

    df['pos-1'] = df[pos_col].apply(lambda x: int(x) - 1)
    df = df[[chr_col, 'pos-1', pos_col, snp_col]]

    _, filename = os.path.split(file)

    df.to_csv(os.path.join(outd, filename+".liftover.bed"), sep="\t", header=False, index=False)

    return os.path.join(outd, filename+".liftover.bed")
