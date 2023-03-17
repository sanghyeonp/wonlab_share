from src.packages import *
from src.util import run_bash, logger

def is_chr(x):
    try:
        chr = int(str(x).replace('chr', ''))
    except:
        return False
    return chr


def merge_lifted(file, delim, snp_col, chr_col, pos_col, 
                unlifted_snplist, keep_initial_pos, keep_unlifted, keep_intermediate,
                outf, outd, verbose=False):
    logs = []
    _ ,filename = os.path.split(file)
    lifted_file = os.path.join(outd, "{}.liftover.lifted".format(filename))
    unlifted_file = os.path.join(outd, "{}.liftover.unlifted".format(filename))

    ## Match bp from liftover
    liftover_out = pd.read_csv(lifted_file, sep='\t', 
                            index_col=False, names=['CHR', 'BP-1', 'BP', 'SNP'])
    liftover_out.drop(columns=['BP-1'], inplace=True)
    liftover_out.columns = ['Chr_lifted', 'Pos_lifted', 'SNP_lifted']
    liftover_out['Chr_lifted'] = liftover_out['Chr_lifted'].apply(lambda x: is_chr(x))
    liftover_out = liftover_out[liftover_out['Chr_lifted'] != False]

    ## keep_unlifted
    if keep_unlifted:
        retain = []
        with open(unlifted_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for idx, row in enumerate(reader):
                if idx % 2 == 1:
                    retain.append([row[0], -9, row[3]])

        unlifted_retained = pd.DataFrame(retain, columns=['Chr_lifted', 'Pos_lifted', 'SNP_lifted'])
        unlifted_retained['Chr_lifted'] = unlifted_retained['Chr_lifted'].apply(lambda x: is_chr(x))
        unlifted_retained = unlifted_retained[unlifted_retained['Chr_lifted'] != False]

        liftover_out = pd.concat([liftover_out, unlifted_retained])
        
    ## Merge
    df = pd.read_csv(file, sep=delim, index_col=False)

    df = pd.merge(df, liftover_out, how="left", left_on=snp_col, right_on="SNP_lifted")

    log_ = "Number of SNPs initially: {:,}".format(len(df)); logger(logs, log_, verbose)
    log_ = "Number of SNPs lifted: {:,}".format(len(df[~df['SNP_lifted'].isna()])); logger(logs, log_, verbose)
    if keep_unlifted:
        log_ = "Number of SNPs unlifted: {:,}".format(len(df[(df['Pos_lifted'] == -9) | (df['SNP_lifted'].isna())])); logger(logs, log_, verbose)
    else:
        log_ = "Number of SNPs unlifted: {:,}".format(len(df[df['SNP_lifted'].isna()])); logger(logs, log_, verbose)

    if unlifted_snplist:
        df[df['SNP_lifted'].isna()][[snp_col]].to_csv("{}.unlifted.snplist".format(filename), sep="\t", index=False, header=False)

    df = df[~df['SNP_lifted'].isna()]

    df.drop(columns=['SNP_lifted'], inplace=True)

    if not keep_initial_pos:
        df.drop(columns=[chr_col, pos_col], inplace=True)
        df.rename(columns={'Chr_lifted':chr_col, 'Pos_lifted':pos_col}, inplace=True)

    if not keep_intermediate:
        run_bash(bash_cmd="rm {}".format(lifted_file))
        run_bash(bash_cmd="rm {}".format(unlifted_file))
        run_bash(bash_cmd="rm {}".format(os.path.join(outd, filename+".liftover.bed")))

    df.to_csv(os.path.join(outd, outf), sep=delim, index=False)

    return logs