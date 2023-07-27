import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *


def read_input(file, delim_in, file_compression,
        infer_chr_pos_ref_alt, chr_col, pos_col, ref_col, alt_col,
        log_list=[]):
    # Read input file
    log_list = logger(log_list, log="## Reading the input file...")
    df_ = pd.read_csv(file, sep=delim_in, index_col=False, compression=file_compression, low_memory=False)
    log_list = logger(log_list, log="Number of SNPs in the input file: {:,}".format(len(df_)))

    # Infer chromosome and position
    if infer_chr_pos_ref_alt:
        # 지금은 CHR:POS:REF/ALT:REF/ALT format이라고 생각하고 코딩한 것.
        # 하지만 CHR:POS 이런 식 일 수도 있음.
        infer_col = infer_chr_pos_ref_alt[0]
        data_format = infer_chr_pos_ref_alt[1]
        separator = infer_chr_pos_ref_alt[2]

        log_list = logger(log_list, log="## Inferring CHR and POS from the column '{}' with data structure '{}'...".format(infer_col, data_format))

        data_format = data_format.split(sep=separator)
        idx = {v:data_format.index(v) for v in ['CHR', 'POS', 'REF', 'ALT']}

        df_[chr_col] = df_[infer_col].apply(lambda x: x.split(sep=separator)[idx['CHR']])
        df_[pos_col] = df_[infer_col].apply(lambda x: x.split(sep=separator)[idx['POS']])
        df_[ref_col] = df_[infer_col].apply(lambda x: x.split(sep=separator)[idx['REF']])
        df_[alt_col] = df_[infer_col].apply(lambda x: x.split(sep=separator)[idx['ALT']])

    df_[chr_col] = df_[chr_col].astype(str)
    # Drop any inappropriate base position
    def check_pos(pos):
        try:
            int(pos)
            return False
        except:
            return True
    df_['Drop'] = df_[pos_col].apply(lambda x: check_pos(x))
    df_ = df_[df_['Drop'] == False]
    df_.drop(columns=['Drop'], inplace=True)
    df_[pos_col] = df_[pos_col].astype(int)
    df_[pos_col] = df_[pos_col].astype(str)
    df_[ref_col] = df_[ref_col].astype(str)
    df_[alt_col] = df_[alt_col].astype(str)

    # Check CHR
    chr_check = 'chr' in df_.loc[0, chr_col]
    if chr_check:
        df_[chr_col] = df_[chr_col].apply(lambda x: x.replace('chr', ''))
    
    return df_, log_list


def make_annovar_input(input_df, file, chr_col, pos_col, ref_col, alt_col, outd):
    # Subset columns
    df = input_df[[chr_col, pos_col, pos_col, ref_col, alt_col]]

    df.columns = ['CHR', 'BP1', 'BP2', 'A2', 'A1']

    with open(os.path.join(outd, os.path.split(file)[-1]+".annovin"), 'w') as f:
        rows = df.values.tolist()
        rows2write =["\t".join([str(v) for v in row]) + "\n" for row in rows]
        f.writelines(rows2write)
    
    with open(os.path.join(outd, os.path.split(file)[-1]+".flip.annovin"), 'w') as f:
        rows = df[['CHR', 'BP1', 'BP2', 'A1', 'A2']].values.tolist()
        rows2write =["\t".join([str(v) for v in row]) + "\n" for row in rows]
        f.writelines(rows2write)

    return os.path.split(file)[-1]+".annovin", os.path.split(file)[-1]+".flip.annovin"
