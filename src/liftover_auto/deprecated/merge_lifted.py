import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from packages import *
from util import run_bash, logger


def convert_to_int(x):
    try:
        x = int(float(x))
        return x
    except:
        return x


def merge_lifted(lifted_file,
                df_input_gwas, snp_col
                ):
    """
    lifted_file (str) : liftOver lifted file path

    df_input_gwas (pd.DataFrame) : pandas dataframe of the input GWAS summary statistics

    input_bed (str) : liftOver input bed file path
    """
    ### Read in the lifted output
    df_lifted = pd.read_csv(lifted_file, sep="\t", index_col=False, 
                                names=['CHR_lifted', 'BP-1', 'POS_lifted', 'SNP'])
    df_lifted.drop(columns=['BP-1'], inplace=True)
    df_lifted['CHR_lifted'] = df_lifted['CHR_lifted'].apply(lambda x: x.replace("chr", ""))

    ### Merge lifted.
    df_infile_merged = pd.merge(df_input_gwas, df_lifted, how="left", left_on=snp_col, right_on='SNP')
    if snp_col != 'SNP':
        df_infile_merged.drop(columns=['SNP'], inplace=True)
    
    ### Fill NA.
    df_infile_merged.fillna({'CHR_lifted':-9, 'POS_lifted':-9}, inplace=True)

    ### Convert to integer.
    df_infile_merged['CHR_lifted'] = df_infile_merged['CHR_lifted'].apply(lambda x: convert_to_int(x))
    df_infile_merged['POS_lifted'] = df_infile_merged['POS_lifted'].apply(lambda x: convert_to_int(x))

    return df_infile_merged
