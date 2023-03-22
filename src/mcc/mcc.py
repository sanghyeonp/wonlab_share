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
    parser.add_argument('--p_col', required=True, 
                        help='P-value column name in the input file.')
    parser.add_argument('--delim', required=False, type=str, default="tab",
                        help='Delimiter of the input file. Default = "tab". Choices = ["tab", "comma", "whitespace"].')
    parser.add_argument('--skip_rows', required=False, type=int, default=0,
                        help='Specify the number of first lines in the input file to skip. Default = 0.')

    parser.add_argument('--outf', required=False, type=str, default='NA',
                        help='Name of output file. Default = mcc.<Input file name>')
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


def main(file, skip_rows, delim, p_col, outf, outd, quoting, log, verbose):
    logs = []
    log_ = "Conducting multiple comparisons correction for:\n\t{}".format(file); logs = log_action(logs, log_, verbose)
    log_ = "Reading the input file...".format(file); logs = log_action(logs, log_, verbose)

    df = pd.read_csv(file, sep=delim, index_col=False, skiprows=skip_rows)

    fdr_pval = fdr_pval_cal(df[p_col].tolist())
    df['FDR P-value'] = fdr_pval
    df['FDR significant'] = df['FDR P-value'].apply(lambda x: pval_sig(x))

    log_ = "Number of comparisons: {:,}".format(len(df)); logs = log_action(logs, log_, verbose)
    log_ = "Bonferroni correction threshold (0.05/{}): {:,}".format(len(df), 0.05/len(df)); logs = log_action(logs, log_, verbose)

    bonferroni_pval = bonferroni_pval_cal(df[p_col].tolist())
    df['Bonferroni P-value'] = bonferroni_pval
    df['Bonferroni significant'] = df['Bonferroni P-value'].apply(lambda x: pval_sig(x))

    df.sort_values(by=['Bonferroni P-value'], inplace=True)

    # Summary
    log_ = "Number of FDR significant: {:,}".format(len(df[df['FDR significant'] == 'Yes'])); logs = log_action(logs, log_, verbose)
    log_ = "Number of Bonferroni significant: {:,}".format(len(df[df['Bonferroni significant'] == 'Yes'])); logs = log_action(logs, log_, verbose)

    if quoting:
        df.to_csv(os.path.join(outd, outf), sep=delim, index=False, quoting=csv.QUOTE_ALL)
    else:
        df.to_csv(os.path.join(outd, outf), sep=delim, index=False)

    if log:
        with open(os.path.join(outd, outf + ".log"), 'w') as f:
            f.writelines([l+'\n' for l in logs])


if __name__ == "__main__":
    args = parse_args()

    if args.outf == "NA":
        args.outf = "mcc.{}".format(os.path.split(args.file))

    if args.outd == "NA":
        args.outd = os.getcwd()
    
    main(file=args.file, 
        skip_rows=args.skip_rows,
        delim=map_delim(args.delim),
        p_col=args.p_col,  
        outf=args.outf, 
        outd=args.outd, 
        quoting=args.quoting, 
        log=args.log, 
        verbose=args.verbose
        )
    