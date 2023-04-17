import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *


def parse_args():
    parser = argparse.ArgumentParser(description=":: Multiple comparisons correction ::")
    parser.add_argument('--file', required=True,
                        help='Path to the input file.')
    parser.add_argument('--p-col', dest="p_col", nargs="+", 
                        help='P-value column name in the input file. If there are multiple P-values to compute, those column name can be specified subsequently. \
                            For example, --p_col P1 P2 P3')
    parser.add_argument('--delim-in', dest="delim_in", required=False, type=str, default="tab",
                        help='Delimiter of the input file. Default = "tab". Choices = ["tab", "comma", "whitespace"].')
    parser.add_argument('--skip-rows', dest="skip_rows", required=False, type=int, default=0,
                        help='Specify the number of first lines in the input file to skip. Default = 0.')

    parser.add_argument('--bonferroni-only', dest='bonferroni_only', action='store_true',
                    help='Specify to conduct Bonferroni correction only. Default = False.')
    parser.add_argument('--fdr-only', dest='fdr_only', action='store_true',
                        help='Specify to conduct FDR correction only. Default = False.')
    

    parser.add_argument('--outf', required=False, type=str, default='NA',
                        help='Name of output file. Default = mcc.<Input file name>')
    parser.add_argument('--delim-out', dest="delim_out", required=False, type=str, default="NA",
                        help='Delimiter of the output file. Default = "NA" then identical delimiter as the input file. Choices = ["tab", "comma", "whitespace"].')
    parser.add_argument('--outd', required=False, type=str, default='NA',
                        help='Directory where the output will be saved. Default = Current working directory.')
    parser.add_argument('--quoting', action='store_true',
                        help='Specify to quote the data in the output file. Default = False.')

    parser.add_argument('--log', action='store_true',
                        help='Specify to save the log. Default = False.')
    parser.add_argument('--verbose', action='store_true',
                        help='Specify to see output from the terminal. Default = False.')

    args = parser.parse_args()
    return args


def pval_sig(p):
    return 'Yes' if p < 0.05 else 'No'


def fdr_pval_cal(p):
    return multipletests(p, alpha=0.05, method='fdr_bh')[1]


def bonferroni_pval_cal(p):
    return multipletests(p, alpha=0.05, method='bonferroni')[1]


def log_action(logs, log_, verbose):
    if verbose:
        print(log_)
    logs.append(log_)
    return logs


def compute_fdr(df, p_col, logs, verbose):
    fdr_pval = fdr_pval_cal(df[p_col].tolist())
    df['FDR {}'.format(p_col)] = fdr_pval
    df['FDR significant {}'.format(p_col)] = df['FDR {}'.format(p_col)].apply(lambda x: pval_sig(x))
    log_ = "Number of FDR significant in {}: {:,}".format(p_col,
                                                        len(df[df['FDR significant {}'.format(p_col)] == 'Yes'])); logs = log_action(logs, log_, verbose)
    return df, logs


def compute_bonferroni(df, p_col, logs, verbose):
    log_ = "Number of comparisons in {}: {:,}".format(p_col, len(df)); logs = log_action(logs, log_, verbose)
    log_ = "Bonferroni correction threshold (0.05/{:,}) in {}: {}".format(len(df), 
                                                                        p_col, 
                                                                        0.05/len(df)); logs = log_action(logs, log_, verbose)

    bonferroni_pval = bonferroni_pval_cal(df[p_col].tolist())
    df['Bonferroni {}'.format(p_col)] = bonferroni_pval
    df['Bonferroni significant {}'.format(p_col)] = df['Bonferroni {}'.format(p_col)].apply(lambda x: pval_sig(x))
    log_ = "Number of Bonferroni significant in {}: {:,}".format(p_col,
                                                                len(df[df['Bonferroni significant {}'.format(p_col)] == 'Yes'])); logs = log_action(logs, log_, verbose)
    return df, logs


def main(file, skip_rows, delim_in, bonferroni, fdr, p_col_list, outf, delim_out, outd, quoting, log, verbose):
    logs = []
    log_ = "Conducting multiple comparisons correction for:\n\t{}".format(file); logs = log_action(logs, log_, verbose)
    log_ = "Reading the input file...".format(file); logs = log_action(logs, log_, verbose)

    df = pd.read_csv(file, sep=delim_in, index_col=False, skiprows=skip_rows)

    if fdr:
        for p_col in p_col_list:
            df, logs = compute_fdr(df, p_col, logs, verbose)

    if bonferroni:
        for p_col in p_col_list:
            df, logs = compute_bonferroni(df, p_col, logs, verbose)

    if fdr is False and bonferroni is False:
        for p_col in p_col_list:
            df, logs = compute_fdr(df, p_col, logs, verbose)
            df, logs = compute_bonferroni(df, p_col, logs, verbose)


    df.sort_values(by=['{}'.format(p_col_list[0])], inplace=True)
        

    if quoting:
        df.to_csv(os.path.join(outd, outf), sep=delim_out, index=False, quoting=csv.QUOTE_ALL)
    else:
        df.to_csv(os.path.join(outd, outf), sep=delim_out, index=False)

    if log:
        with open(os.path.join(outd, outf + ".log"), 'w') as f:
            f.writelines([l+'\n' for l in logs])


if __name__ == "__main__":
    args = parse_args()

    if args.outf == "NA":
        args.outf = "mcc.{}".format(os.path.split(args.outf)[-1])

    if args.outd == "NA":
        args.outd = os.getcwd()
    
    if args.delim_out == "NA":
        args.delim_out = args.delim_in
    
    main(file=args.file, 
        skip_rows=args.skip_rows,
        delim_in=map_delim(args.delim_in),
        p_col_list=args.p_col, 
        bonferroni=args.bonferroni_only, 
        fdr=args.fdr_only,
        outf=args.outf, 
        delim_out=map_delim(args.delim_out),
        outd=args.outd, 
        quoting=args.quoting, 
        log=args.log, 
        verbose=args.verbose
        )
    