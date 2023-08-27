import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *
from shared_data import VCF_1kGp3_modified_dir
AF_SOURCE = {'1kgp3': VCF_1kGp3_modified_dir
            }

from timeit import default_timer as timer
from datetime import timedelta
tqdm.pandas(leave=False, bar_format='{l_bar}{bar:30}{r_bar}{bar:-30b}')


def parse_args():
    parser = argparse.ArgumentParser(description=":: Map AF by chromosome and base position ::")
    
    # Required arguments
    parser.add_argument('--file', required=True,
                        help='Path to the input file.')
    parser.add_argument('--chr-col', dest="chr_col", required=True,
                        help="Name of the chr column in the input file. Must be in genome build GRCh37.")
    parser.add_argument('--pos-col', dest="pos_col", required=True,
                        help="Name of the pos column in the input file. Must be in genome build GRCh37.")
    # parser.add_argument('--rsid-col', dest="rsid_col", required=True,
    #                     help="Name of the rsID column in the input file.")
    parser.add_argument('--a1-col', dest="a1_col", required=True,
                        help="Name of the a1 column, which is effect allele, in the input file.")
    parser.add_argument('--ancestry', required=True, nargs='+',
                        help="Specify ancestry (or ancestries if MAF from multiple ancestries are required.).\
                            Choices = ['GLOBAL', 'EAS', 'EUR', 'AMR', 'AFR', 'SAS'].\
                            EAS: East-asian; EUR: European; AMR: American; AFR: African; SAS: South-asian."
                        )
    parser.add_argument('--source', required=False, default="1kgp3",
                        help="Source to annotate MAF. 아직 1000 Genome phase 3 밖에 없습니다. Choices = ['1kgp3'].")
    
    
    # Optional arguments
    parser.add_argument('--delim-in', dest="delim_in", required=False, default="tab",
                        help="Delimiter used in the input file. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.")
    parser.add_argument('--compression-in', dest="compression_in", required=False, default="NA",
                        help="Specify compression type from the following ['zip', 'gzip', 'bz2', 'zstd', 'tar']. Default='NA'.")
    parser.add_argument('--chromosome-filter', dest="chromosome_filter", required=False, default="autosome",
                        help="Specify which chromosomes to include. Default = 'autosome'. If 'all' is specified, then chromosome X, Y, MT are included. \
                        Choices = ['autosome', 'all', either from 1 to 23]")

    parser.add_argument('--outf', required=False, default="NA",
                        help="Specify the name of the output file. Default = 'af.<file>'.")
    parser.add_argument('--outd', required=False, default="NA",
                        help="Specify the path to output directory. Default = Current working directory.")
    parser.add_argument('--delim-out', dest="delim_out", required=False, default="NA",
                        help="Delimiter for the output file. Choices = ['NA', 'tab', 'comma', 'whitespace']. If 'NA', identical delimiter as input delimiter will be used. Default = 'NA'.")
    parser.add_argument('--compression-out', dest="compression_out", required=False, default="NA",
                        help="Specify compression type from the following ['NA', 'zip', 'gzip', 'bz2', 'zstd', 'tar']. Default='NA'.")


    args = parser.parse_args()
    return args


# < VCF파일에 rsID mappind하면 rsid_col argument 추가한거 사용하기
# def main(file, delim_in, compression_in, 
#         chr_col, pos_col, rsid_col, ancestry, source, chromosome_filter,
#         outf, outd, delim_out, compression_out
#         ):
# > VCF파일에 rsID mappind하면 rsid_col argument 추가한거 사용하기
def main(file, delim_in, compression_in, 
        chr_col, pos_col, a1_col, ancestry, source, chromosome_filter,
        outf, outd, delim_out, compression_out
        ):
    global AF_SOURCE
    
    ### Read input file
    start0 = timer()
    print(":: Reading input file ::")
    if compression_in == "NA":
        compression_in = None
    df_ = pd.read_csv(file, sep=delim_in, index_col=False, compression=compression_in, dtype=str)

    print("Elapsed: {}".format(timedelta(seconds=timer() - start0)))


    ### Prepare annotation file path
    print(":: Prepare annotation file path ::")
    if source == '1kgp3':
        print("MAF from 1000 Genome phase 3 VCF")
        vcf_files = [os.path.join(AF_SOURCE[source], f) for f in os.listdir(AF_SOURCE[source])]
    
    if chromosome_filter == 'autosome':
        vcf_files = [f for f in vcf_files if not isin_list(f.split(sep="chr")[-1], ['X', 'Y', 'MT'])]
    
    if chromosome_filter.isdigit():
        vcf_files = [f for f in vcf_files if chromosome_filter == f.split(sep="chr")[-1].replace(".tsv.gz", "")]

    # code.interact(local=dict(globals(), **locals()))

    ### Columns to retain in annotation dataframe
    ancestry_col_mapper = {'GLOBAL':'AF_global', 
                            'EAS':'AF_EAS', 
                            'EUR':'AF_EUR', 
                            'AMR':'AF_AMR', 
                            'AFR':'AF_AFR', 
                            'SAS':'AF_SAS'
                            }
    col_retain = [ancestry_col_mapper[anc] for anc in ancestry]

    ### Annotation
    start1 = timer()
    print(":: Annotating AF ::")

    df = None
    for idx, vcf in enumerate(tqdm(vcf_files, desc="Chromosome", leave=False, bar_format='{l_bar}{bar:30}{r_bar}{bar:-30b}')):
        df_annot = pd.read_csv(vcf, sep="\t", index_col=False, compression="gzip", dtype=str)
        df_annot.rename(columns={'#CHROM':'CHR_1kG',
                                'POS':'POS_1kG',
                                'REF':'REF_1kG',
                                'ALT':'ALT_1kG'}, inplace=True)

        df_annot = df_annot[['CHR_1kG', 'POS_1kG', 'REF_1kG', 'ALT_1kG'] + col_retain]
        if idx == 0:
            df = pd.merge(df_, df_annot, how="left", left_on=[chr_col, pos_col], right_on=['CHR_1kG', 'POS_1kG'])
            df.drop(columns=['CHR_1kG', 'POS_1kG'], inplace=True)
        else:            
            df = pd.merge(df, df_annot, how="left", left_on=[chr_col, pos_col], right_on=['CHR_1kG', 'POS_1kG'], suffixes=['', '_'])
            df.drop(columns=['CHR_1kG', 'POS_1kG'], inplace=True)
            # code.interact(local=dict(globals(), **locals()))
            for col in col_retain + ['REF_1kG', 'ALT_1kG']:
                df[col] = df[col].fillna(df[col + '_'])
                df.drop(columns=[col + "_"], inplace=True)
            
    del df_
    print("Elapsed: {}".format(timedelta(seconds=timer() - start1)))

    ### Handle missing value
    start2 = timer()
    print(":: Handling missing values ::")
    for col in col_retain:
        df[col] = df[col].fillna('-9')
    
    df['REF_1kG'] = df['REF_1kG'].fillna('NA')
    df['ALT_1kG'] = df['ALT_1kG'].fillna('NA')
    print("Elapsed: {}".format(timedelta(seconds=timer() - start2)))

    ### Align allele frequency to a1 column
    start3 = timer()
    print(":: Aligning allele frequency to A1 column ::")
    def align(a1, alt, af):
        if ("," in af) or (alt == "NA") or (a1 == alt): 
            # Multi-allelic인 경우, 따로 aligning 하지 않기.
            # Reference에 없는 SNP이면 -9로 return.
            # Reference ALT랑 A1이 matching되면 AF 그대로 return.
            return af
        return str(1 - float(af))

    for col in tqdm(col_retain, desc="Align", leave=False, bar_format='{l_bar}{bar:30}{r_bar}{bar:-30b}'):
        df[col + "_aligned"] = df.apply(lambda row: align(row[a1_col], row['ALT_1kG'], row[col]), axis=1)
        df.drop(columns=[col], inplace=True)
        df.rename(columns={col + "_aligned":col}, inplace=True)

    print("Elapsed: {}".format(timedelta(seconds=timer() - start3)))

    ### Save the result
    start4 = timer()
    print(":: Saving the annotated result ::")
    if outf == "NA":
        outf = "af." + os.path.split(file)[-1]
    if outd == "NA":
        outd = "."
    if compression_out == "NA":
        compression_out = None
    print("Path: {}".format(os.path.join(outd, outf)))

    df.to_csv(os.path.join(outd, outf), sep=delim_out, 
                index=False, compression=compression_out)
    print("Elapsed: {}".format(timedelta(seconds=timer() - start4)))


if __name__ == "__main__":
    args = parse_args()

    if args.delim_out == "NA":
        args.delim_out = args.delim_in

    main(file=args.file,
        delim_in=map_delim(args.delim_in),
        compression_in=args.compression_in,
        chr_col=args.chr_col,
        pos_col=args.pos_col,
        a1_col=args.a1_col,
        chromosome_filter=args.chromosome_filter,
        ancestry=args.ancestry,
        source=args.source,
        outf=args.outf,
        outd=args.outd,
        delim_out=map_delim(args.delim_out),
        compression_out=args.compression_out
        )
