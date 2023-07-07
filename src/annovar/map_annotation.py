import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *


def combine_annov_out(annov_input1, annov_input2, outd, log_list=[]):
    log_list = logger(log_list, log="## Merging multianno.txt and flipped multianno.txt...")
    df_annov1 = pd.read_csv(os.path.join(outd, annov_input1), sep="\t", index_col=False, dtype=str)
    df_annov1 = df_annov1.loc[df_annov1['avsnp150'] != '.', ]
    log_list = logger(log_list, log="Number of SNPs annotated with multianno.txt: {:,}".format(len(df_annov1)))
    
    if not df_annov1.empty:
        df_annov1['chr_pos_ref_alt_new'] = df_annov1.apply(lambda row: "{}:{}:{}:{}".format(row['Chr'], row['End'], row['Ref'], row['Alt']), axis=1)

    df_annov2 = pd.read_csv(os.path.join(outd, annov_input2), sep="\t", index_col=False, dtype=str)
    df_annov2 = df_annov2.loc[df_annov2['avsnp150'] != '.', ]

    df_annov = None
    if not df_annov2.empty:
        # This is flipped version. Initially REF and ALT should be ALT and REF, respectively.
        df_annov2.columns = ['Chr', 'Start', 'End', 'Alt', 'Ref', 'Func.refGene', 'Gene.refGene',
                            'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
                            'avsnp150']
        df_annov2['chr_pos_ref_alt_new'] = df_annov2.apply(lambda row: "{}:{}:{}:{}".format(row['Chr'], row['End'], row['Ref'], row['Alt']), axis=1)
        
        # Get annotated SNPs not in df_annov1.
        # If annotated SNP is present in df_annov2 after filtering the ones in df_annov1, it means this SNP has flipped allele. 
        if not df_annov1.empty:
            df_annov2 = df_annov2.loc[~df_annov2['chr_pos_ref_alt_new'].isin(df_annov1['chr_pos_ref_alt_new'].tolist()), ]

        if not df_annov2.empty:
            log_list = logger(log_list, log="Number of SNPs annotated with flipped multianno.txt: {:,}".format(len(df_annov2)))

            df_annov1['flipped'] = '.'
            df_annov2['flipped'] = '1'

            df_annov = pd.concat([df_annov1, df_annov2])
    
    if df_annov is None and df_annov1.empty:
        return None, None, False

    if df_annov is None:
        df_annov1['flipped'] = '.'
        df_annov = df_annov1

    log_list = logger(log_list, log="Total number of SNPs annotated: {:,}".format(len(df_annov)))

    df_annov = df_annov[['chr_pos_ref_alt_new', 'avsnp150', 'flipped']]
    df_annov.columns = ['chr_pos_ref_alt_new', 'rsID_annov', 'flipped_annov']

    return df_annov, log_list, True


def map_annovar_out(input_df, annov_input1, annov_input2, chr_col, pos_col, ref_col, alt_col, outd, log_list=[]):
    annov_output1 = annov_input1.replace(".annovin", ".hg19_multianno.txt")
    annov_output2 = annov_input2.replace(".annovin", ".hg19_multianno.txt")
    df_annov, log_list, have_result = combine_annov_out(annov_output1, annov_output2, outd)

    if not have_result:
        return None, None, False

    log_list = logger(log_list, log="## Mapping annotated SNPs to the input file...")
    input_df['chr_pos_ref_alt'] = input_df.apply(lambda row: "{}:{}:{}:{}".format(row[chr_col], row[pos_col], row[ref_col], row[alt_col]), axis=1)

    output_df = input_df.merge(df_annov, how="left", left_on="chr_pos_ref_alt", right_on="chr_pos_ref_alt_new")

    return output_df, log_list, have_result

