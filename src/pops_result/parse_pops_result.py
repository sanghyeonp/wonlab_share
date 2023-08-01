import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *

FILE_GENE_ANNOT = "/data1/sanghyeon/tutorials/PoPs_tut/pops/example/data/utils/gene_annot_jun10.txt"

def parse_args():
    parser = argparse.ArgumentParser(description=":: Parse PoPS result ::")
    ### PoPS result directory parser
    parser.add_argument('--dir-pops-out', dest="dir_pops_out", required=True,
                        help='Path to the output directory of PoPS.')
    
    ### Lead SNP file parser
    parser.add_argument('--file-leadsnp', dest="file_leadsnp", required=True,
                        help='Path to the file that contains information of the lead SNPs. Best if lead SNP file from FUMA leveraged directly.')
    parser.add_argument('--delim-in', dest="delim_in", required=False, default="tab",
                    help="Delimiter used in the input lead SNP file. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.")
    parser.add_argument('--snp-col', dest="snp_col", required=False, default="rsID",
                    help="Name of the chromosome column in the input lead SNP file. Default = 'rsID'.")
    parser.add_argument('--chr-col', dest="chr_col", required=False, default="chr",
                    help="Name of the chromosome column in the input lead SNP file. Default = 'chr'.")
    parser.add_argument('--pos-col', dest="pos_col", required=False, default="pos",
                    help="Name of the base position column in the input lead SNP file. Default = 'pos'.")
    parser.add_argument('--window', required=False, default=500000,
                    help="Window used to define a locus centered by the lead SNP. Default = 500000.")
    
    parser.add_argument('--out-prefix', dest="out_prefix", required=False, default="NA",
                        help="Specify the prefix of the output file name. '<output prefix>.pops_result.*'. Default = 'trait.pops_result.*'.")
    parser.add_argument('--outd', required=False, default="NA",
                        help="Specify the path to output directory. Default = Current working directory.")

    args = parser.parse_args()
    return args


def main(dir_pops_out,
        file_leadsnp, delim_in, snp_col, chr_col, pos_col, window,
        out_prefix, outd):
    log_list = []

    msg = ":: Parsing PoPS result ::\n\tPoPS result directory: {}".format(dir_pops_out); log_list.append(msg); print(msg)

    ### Information on feature selection
    file_feature = [f for f in os.listdir(dir_pops_out) if ".marginals" in f][0]
    df_feature = pd.read_csv(file_feature, sep="\t", index_col=False)

    df_feature.columns = ['feature', 'beta', 'se', 'pval', 'r2', 'selected']

    print("Total number of features included: {:,}\n\tNumber of features selected: {:,} ({:.1%})\n\tNumber of features discarded: {:,} ({:.1%})".format(
                len(df_feature),
                len(df_feature[df_feature['selected'] == True]),
                len(df_feature[df_feature['selected'] == True]) / len(df_feature),
                len(df_feature[df_feature['selected'] == False]),
                len(df_feature[df_feature['selected'] == False]) / len(df_feature)
                ))

    ### Read lead SNP file
    df_leadsnp = pd.read_csv(file_leadsnp, sep=delim_in, index_col=False)
    df_leadsnp = df_leadsnp[[snp_col, chr_col, pos_col]]
    df_leadsnp[chr_col] = df_leadsnp[chr_col].astype(int)
    df_leadsnp[pos_col] = df_leadsnp[pos_col].astype(int)
    df_leadsnp['cis region'] = df_leadsnp.progress_apply(lambda row: pd.Interval(row[pos_col] - window, row[pos_col] + window), axis=1)
    n_leadsnp = len(df_leadsnp)
    msg = "Number of lead SNPs: {:,}".format(n_leadsnp); log_list.append(msg); print(msg)
    
    ### Gene annotation 정보 가져오기.
    df_gene_annot = pd.read_csv(FILE_GENE_ANNOT, sep="\t", index_col=False)
    df_gene_annot['CHR'] = df_gene_annot['CHR'].astype(int)
    df_gene_annot['START'] = df_gene_annot['START'].astype(int)
    df_gene_annot['END'] = df_gene_annot['END'].astype(int)
    df_gene_annot['gene region'] = df_gene_annot.apply(lambda row: pd.Interval(row['START'], row['END']), axis=1)

    ### Find overlapping genes in cis-region of the lead SNP
    row_list = []
    for _, row in tqdm(df_leadsnp.iterrows(), leave=False):
        chr = row[chr_col]
        df_gene_annot_filter = df_gene_annot[df_gene_annot['CHR'] == chr]
        df_gene_annot_filter['overlap'] = df_gene_annot_filter['gene region'].progress_apply(lambda x: x.overlaps(row['cis region']))
        df_gene_overlap = df_gene_annot_filter[df_gene_annot_filter['overlap'] == True]
        
        genes_symbol = ""
        genes_ensgid = ""
        n_genes = 0
        if not df_gene_overlap.empty:
            genes_symbol = ";".join(df_gene_overlap['NAME'].tolist())
            genes_ensgid = ";".join(df_gene_overlap['ENSGID'].tolist())
            n_genes = len(df_gene_overlap)

        row['overlap_gene_N'] = n_genes
        row['overlap_gene_name'] = genes_symbol
        row['overlap_gene_ENSGID'] = genes_ensgid

        row_list.append(row)

    df_leadsnp_genemapped = pd.DataFrame(row_list)

    msg = "Number of lead SNP with no overlapping genes in the locus: {:,}".format(len(df_leadsnp_genemapped[df_leadsnp_genemapped['overlap_gene_N'] == 0])); log_list.append(msg); print(msg)

    ### Read PoPS result
    file_pops = [file for file in os.listdir(dir_pops_out) if ".preds" in file][0]
    df_pops = pd.read_csv(os.path.join(dir_pops_out, file_pops), sep="\t", index_col=False)

    ### Map PoP score for each genes
    row_list = []
    for idx, row in tqdm(df_leadsnp_genemapped.iterrows(), leave=False):
        if row['overlap_gene_N'] > 0:
            overlapping_gene_list = row['overlap_gene_ENSGID'].split(sep=";")
            df_pops_filter = df_pops[df_pops['ENSGID'].isin(overlapping_gene_list)]
            gene_pops_list = [float(df_pops_filter[df_pops_filter['ENSGID'] == gene]['PoPS_Score'].values[0]) 
                                for gene in overlapping_gene_list]
            row['PoPS_score'] = ";".join([str(score) for score in gene_pops_list])

            max_pops_score = max(gene_pops_list)
            max_idx = gene_pops_list.index(max_pops_score)
            max_pops_score_gene = overlapping_gene_list[max_idx]

            row['Max_PoPS_score'] = max_pops_score
            if max_pops_score > 0:
                prioritized = True
            else:
                prioritized = False
            row['Prioritized'] = prioritized
            if prioritized:
                row['Prioritized_gene_name'] = df_gene_annot[df_gene_annot['ENSGID'] == max_pops_score_gene]['NAME'].values[0]
                row['Prioritized_gene_ENSGID'] = max_pops_score_gene
            else:
                row['Prioritized_gene_name'] = "."
                row['Prioritized_gene_ENSGID'] = "."

        else:
            row['PoPS_score'] = "."
            row['Prioritized'] = False
            row['Max_PoPS_score'] = "."
            row['Prioritized_gene_name'] = "."
            row['Prioritized_gene_ENSGID'] = "."

        row_list.append(row)

    df_result = pd.DataFrame(row_list)
    df_result.drop(columns=['cis region'], inplace=True)

    out_path_all_result = os.path.join(outd, out_prefix + ".pops_result_all.csv")
    with open(out_path_all_result, 'w') as f:
        rows = df_result.values.tolist()
        col_list = df_result.columns.tolist()
        rows2write = [",".join(col_list) + "\n"] + [",".join([str(v) for v in row]) + "\n" for row in rows]
        f.writelines(rows2write)
        msg = "All results from PoPS written at {}".format(out_path_all_result); log_list.append(msg); print(msg)
    
    ### Deal with only prioritized variant
    df_result_prioritized = df_result[df_result['Prioritized'] == True].copy()
    msg = "Number of lead SNP with prioritized genes: {:,} ({:.1%})".format(len(df_result_prioritized),
                                                                            len(df_result_prioritized) / n_leadsnp); log_list.append(msg); print(msg)
    df_result_prioritized.sort_values(by=['Max_PoPS_score'], ascending=[False], inplace=True)
    out_path_prioritized = os.path.join(outd, out_prefix + ".pops_result_prioritized.csv")
    with open(out_path_prioritized, 'w') as f:
        rows = df_result_prioritized.values.tolist()
        col_list = df_result_prioritized.columns.tolist()
        rows2write = [",".join(col_list) + "\n"] + [",".join([str(v) for v in row]) + "\n" for row in rows]
        f.writelines(rows2write)
        msg = "Prioritized results from PoPS written at {}".format(out_path_prioritized); log_list.append(msg); print(msg)
    
    ### Save log
    out_path_log = os.path.join(outd, out_prefix + ".pops_result.log")
    with open(out_path_log, 'w') as f:
        f.writelines([row + "\n" for row in log_list])
        msg = "Parsing log written at {}".format(out_path_log)


if __name__ == "__main__":
    args = parse_args()

    log_list = []

    if args.out_prefix == "NA":
        args.out_prefix = "trait"

    if args.outd == "NA":
        args.outd = os.getcwd()
    
    main(dir_pops_out=args.dir_pops_out,
        file_leadsnp=args.file_leadsnp, 
        delim_in=map_delim(args.delim_in), 
        snp_col=args.snp_col, 
        chr_col=args.chr_col, 
        pos_col=args.pos_col,
        window=args.window,
        out_prefix=args.out_prefix, 
        outd=args.outd
        )
