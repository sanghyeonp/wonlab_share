import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *


def read_input(file, delim_in, file_compression,
        infer_chr_pos_ref_alt, chr_col, pos_col, ref_col, alt_col):
    # Read input file
    print("## Reading the input file...")
    df_ = pd.read_csv(file, sep=delim_in, index_col=False, compression=file_compression, low_memory=False)
    print("Number of SNPs in the input file: {:,}".format(len(df_)))

    # Infer chromosome and position
    if infer_chr_pos_ref_alt:
        # 지금은 CHR:POS:REF/ALT:REF/ALT format이라고 생각하고 코딩한 것.
        # 하지만 CHR:POS 이런 식 일 수도 있음.
        infer_col = infer_chr_pos_ref_alt[0]
        data_format = infer_chr_pos_ref_alt[1]
        separator = infer_chr_pos_ref_alt[2]

        print("## Inferring CHR and POS from the column '{}' with data structure '{}'...".format(infer_col, data_format))

        data_format = data_format.split(sep=separator)
        idx = {v:data_format.index(v) for v in ['CHR', 'POS', 'REF', 'ALT']}

        df_[[chr_col, pos_col, ref_col, alt_col]] = df_.apply(lambda row: [row[infer_col].split(sep=separator)[idx['CHR']], 
                                                    row[infer_col].split(sep=separator)[idx['POS']],
                                                    row[infer_col].split(sep=separator)[idx['REF']],
                                                    row[infer_col].split(sep=separator)[idx['ALT']]
                                                    ], axis=1, result_type='expand')

    df_[chr_col] = df_[chr_col].astype(str)
    df_[pos_col] = df_[pos_col].astype(str)
    df_[ref_col] = df_[ref_col].astype(str)
    df_[alt_col] = df_[alt_col].astype(str)

    # Check CHR
    chr_check = 'chr' in df_.loc[0, chr_col]
    if chr_check:
        df_[chr_col] = df_[chr_col].apply(lambda x: x.replace('chr', ''))
    
    return df_


def make_annovar_input(input_df, file, chr_col, pos_col, ref_col, alt_col):
    # Subset columns
    df = input_df[[chr_col, pos_col, pos_col, ref_col, alt_col]]

    df.columns = ['CHR', 'BP', 'BP', 'A2', 'A1']

    df.to_csv(os.path.split(file)[-1]+".annovin", sep="\t", index=False, header=False)
    df[['CHR', 'BP', 'BP', 'A1', 'A2']].to_csv(os.path.split(file)[-1]+".flip.annovin", sep="\t", index=False, header=False)

    return os.path.split(file)[-1]+".annovin", os.path.split(file)[-1]+".flip.annovin"