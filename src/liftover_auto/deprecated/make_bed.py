import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from packages import *

def convert_to_int(x):
    try:
        x = int(float(x))
        return x
    except:
        return x


def make_bed(file, file_compression, delim, 
            outd,
            snp_col='NA', chr_col='CHR', pos_col='POS',
            infer_col=None):
    
    ### File name extraction.
    _, filename = os.path.split(file)
    
    ### Read the input file
    df_ = pd.read_csv(file, 
                    sep=delim, 
                    index_col=False, 
                    na_values=["NA", "NaN", " "], 
                    compression=file_compression)

    ### Handle infer column.
    if infer_col:
        # If `variant` column has `2:123423:SNP`, then specify as follows: --infer-col variant CHR:POS:X : CHR,POS
        col_to_infer = infer_col[0]
        data_format = infer_col[1]
        separator = infer_col[2]
        col_to_take = infer_col[3].split(sep=",")

        data_format = data_format.split(sep=separator)
        infercol_idx_map = {v:data_format.index(v) for v in col_to_take}

        for inferring_col, idx in infercol_idx_map.items(): 
            df_[inferring_col] = df_[col_to_infer].apply(lambda x: x.split(sep=separator)[idx])
        
        snp_col = col_to_infer
        chr_col = 'CHR'; pos_col; 'POS'

    ### Double check chromosome column.
    df_copy = df_.copy(deep=True)
    df_[chr_col] = df_[chr_col].astype(str)
    df_[chr_col] = df_[chr_col].apply(lambda x: x.upper())

    ### Handle SNP column.
    if snp_col == 'NA':
        df_['SNP'] = df_.apply(lambda row:"{}:{}".format(str(row[chr_col]), str(row[pos_col])), axis=1)
        snp_col = 'SNP'

    ### Fill NAN in SNP column with chr:pos
    df_[snp_col].fillna("NA", inplace=True)
    df_[snp_col] = df_.apply(lambda row: "{}:{}".format(row[chr_col], row[pos_col]) if row[snp_col] == 'NA' else row[snp_col], axis=1)

    ### Drop NA
    # CHR nad POS가 'NAN'인 경우. => Line51에서 string으로 변경하였었음.
    # df_na = df_[(df_[chr_col] == "NAN") | (df_[pos_col].isna())]
    df_ = df_[~(df_[chr_col] == "NAN") | ~(df_[pos_col].isna())]

    ### Handle duplicates.
    # Check if duplicates are present.
    duplicates = df_[snp_col].duplicated()
    has_dup = sum(duplicates) > 0
    if has_dup:
        # df_dropped_duplicates = df_[df_[snp_col].duplicated(keep='first')]
        df_.drop_duplicates(subset=[snp_col], keep="first", inplace=True)


    ### Add chr to chromosome column.
    df_['CHR_new'] = df_[chr_col].apply(lambda x: "chr{}".format(convert_to_int(x)))

    ### Convert to integer
    df_[pos_col] = df_[pos_col].apply(lambda x: convert_to_int(x))

    ### Add POS-1 column.
    df_['POS-1'] = df_[pos_col].apply(lambda x: convert_to_int(int(x) - 1))

    ### Retain only CHR, POS-1, POS, SNP columns
    df = df_[['CHR_new', 'POS-1', pos_col, snp_col]]

    ### Save bed format liftOver input file.
    bed_input = os.path.join(outd, "{}.bed".format(filename))
    with open(bed_input, 'w') as f:
        row_list = df.values.tolist()
        rows2write = ["\t".join([str(v) for v in row]) + "\n" for row in row_list]
        f.writelines(rows2write)
    
    ### Return bed format liftOver input file path.
    return bed_input, df_copy, snp_col, chr_col, pos_col
