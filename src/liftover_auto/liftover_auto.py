import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from packages import *
from util import logger, save_log, map_delim, run_bash
from make_bed import make_bed, convert_to_int
from run_liftover import run_liftover, liftOver_log, reformat_unlifted
from merge_lifted import merge_lifted



def parse_args():
    parser = argparse.ArgumentParser(description=":: liftOver automation ::")

    # Input
    parser.add_argument('--file', required=True,
                        help='Path to the input file being lifted over.')
    parser.add_argument('--in-compression', dest="in_compression", required=False, default="infer",
                        help="Specify compression type from the following ['zip', 'gzip', 'bz2', 'zstd', 'tar']. Default='infer'.")
    parser.add_argument('--delim', required=False, default="tab",
                    help="Delimiter used in the input file. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.")

    # Column information
    parser.add_argument('--snp-col', dest="snp_col", required=False, default="NA",
                    help="Name of the SNP column in the input file. Default = 'SNP'.")
    parser.add_argument('--chr-col', dest="chr_col", required=False, default="CHR",
                    help="Name of the chromosome column in the input file. Default = 'CHR'.")
    parser.add_argument('--pos-col', dest="pos_col", required=False, default="POS",
                    help="Name of the base position column in the input file. Default = 'POS'.")

    # Infer column information
    parser.add_argument('--infer-col', dest="infer_col", nargs='+', required=False,
                        help="Specify `column name`, `data format`, `separator`, and `columns to infer` to infer necessary columns from the specified column.\
                            For example, column named 'variant' have variant name '2:179816990:C:T' where chromosome and position can be inferred as 2 and 179816990, respectively.\
                            Then, specify as follows: --infer-col variant CHR:POS:REF:ALT : CHR,POS,REF,ALT\
                            If `variant` column has `2:123423:SNP`, then specify as follows: --infer-col variant CHR:POS:X : CHR,POS")

    # New column names after litfOver
    parser.add_argument('--chr-col-new', dest="chr_col_new", required=False, default="CHR_lifted",
                    help="New column name for chromosome in the lifted file. Default = 'CHR_lifted'.")
    parser.add_argument('--pos-col-new', dest="pos_col_new", required=False, default="POS_lifted",
                    help="New column name for base position in the lifted file. Default = 'POS_lifted'.")

    # liftOver options.
    parser.add_argument('--build-from', dest="build_from", required=True, type=int,
                    help="Genome build number initial. Choices = [36, 37, 38].")
    parser.add_argument('--build-to', dest="build_to", required=True, type=int,
                    help="Genome build number after performing liftover. Choices = [36, 37, 38].")
    
    # Output options.
    parser.add_argument('--outf', required=False, default="NA",
                        help="Specify the name of the output file. Default = 'lifted.<file>'.")
    parser.add_argument('--outd', required=False, default="NA",
                        help="Specify the path to output directory. Default = Current working directory.")

    # Other optional.
    parser.add_argument('--verbose', action='store_true',
                help='Specify see the summary in the terminal. Default = False')
    parser.add_argument('--do-not-save-unlifted', action='store_true',
                help='Specify to avoid saving unlifted SNPs and corresponding errors. Default = False')
    parser.add_argument('--save-mapping-file', action='store_true',
                help='Specify to save mapping file that can be used to map genomic positions of GWAS from the identical provider. Default = False')
    parser.add_argument('--rm-intermediate-file', action='store_true',
                help='Specify to delete all the intermediate files (i.e., bed file, lifted, and unlifted files). Default = False')
    parser.add_argument('--do-not-save-lifted-gwas', action='store_true',
                help='Specify to avoid saving a new GWAS with the lifted result. Default = False')
    parser.add_argument('--drop-pos-build-before', action='store_true',
                help='Specify to drop the initial build genomic position in the new GWAS summary statistics. Default = False')

    args = parser.parse_args()
    return args


def main(file, in_compression, delim,
        snp_col, chr_col, pos_col,
        infer_col,
        chr_col_new, pos_col_new,
        build_from, build_to,
        outf, outd,
        verbose,
        do_not_save_unlifted,
        save_mapping_file,
        rm_intermediate_file,
        do_not_save_lifted_gwas,
        drop_pos_build_before):
    
    ### Initialize the logger
    log_list = []
    log_ = ":: liftOver INPUT ::\n\
    Input: {}\n\
    Genome build from: GRCh{}\n\
    Genome build to: GRCh{}\n\
    --do-not-save-unlifted: {}\n\
    --save-mapping-file: {}\n\
    --rm-intermediate-file: {}\n\
    --do-not-save-lifted-gwas: {}\n\
    --drop-pos-build-before: {}".format(file, build_from, build_to, 
                                            do_not_save_unlifted, save_mapping_file,
                                            rm_intermediate_file, do_not_save_lifted_gwas,
                                            drop_pos_build_before); log_list = logger(log_list, log_, verbose = verbose)

    ### Make bed format liftOver input file
    input_bed_file, df_input_file, \
        snp_col, chr_col, pos_col = make_bed(file=file, 
                                            file_compression=in_compression, 
                                            delim=delim, 
                                            outd=outd,
                                            snp_col=snp_col, 
                                            chr_col=chr_col, 
                                            pos_col=pos_col,
                                            infer_col=infer_col
                                            )

    ### Run liftOver.
    lifted_file, unlifted_file = run_liftover(input_bed_file=input_bed_file, 
                                            build_from=build_from, 
                                            build_to=build_to, 
                                            outd=outd
                                            )
    
    log_ = liftOver_log(input_bed_file=input_bed_file, 
                                lifted_file=lifted_file, 
                                unlifted_file=unlifted_file
                                ); log_list = logger(log_list, log_, verbose = verbose)

    ### Merge liftOver output.
    df_merged = merge_lifted(lifted_file=lifted_file,
                            df_input_gwas=df_input_file, 
                            snp_col=snp_col
                            )
    df_merged[chr_col] = df_merged[chr_col].apply(lambda x: convert_to_int(x))
    df_merged[pos_col] = df_merged[pos_col].apply(lambda x: convert_to_int(x))

    ### Additional manipulation and save the file.
    log_saved_file = []

    # Save unlifted SNPs
    if not do_not_save_unlifted:
        unlifted_rows = reformat_unlifted(unlifted_file=unlifted_file)
        unlifted_reformat_file = os.path.join(outd, "unlifted_reformat.{}".format(outf)); log_saved_file.append("Unlifted reformat file: {}".format(unlifted_reformat_file))
        with open(unlifted_reformat_file, 'w') as f:
            f.writelines(['\t'.join(['CHR', 'POS', 'SNP', 'ERROR']) + '\n'] + ['\t'.join(row) + '\n' for row in unlifted_rows])

    # Save the mapping file
    if save_mapping_file:
        df_temp = df_merged[[snp_col, chr_col, pos_col, 'CHR_lifted', 'POS_lifted']]
        df_temp.rename(columns={
                            'CHR_lifted':'CHR_b{}'.format(build_to),
                            'POS_lifted':'POS_b{}'.format(build_to)
                            }, inplace=True
                        )
        df_temp.drop_duplicates(subset=[snp_col], keep="first", inplace=True)

        mapping_file = os.path.join(outd, "mapping.{}".format(outf)); log_saved_file.append("Mapping file: {}".format(mapping_file))
        with open(mapping_file, 'w') as f:
            row_list = df_temp.values.tolist()
            col_list = df_temp.columns.tolist()
            rows2write = ['\t'.join(col_list) + "\n"] + ['\t'.join([str(v) for v in row]) + "\n" for row in row_list]
            f.writelines(rows2write)
    
    # Delete all the intermediate files: input bed, lifted, unlifted
    if rm_intermediate_file:
        run_bash(bash_cmd="rm {}".format(lifted_file))
        run_bash(bash_cmd="rm {}".format(unlifted_file))
        run_bash(bash_cmd="rm {}".format(input_bed_file))

    # Save the lifted GWAS summary statistics
    if not do_not_save_lifted_gwas:
        df_merged.rename(columns={"CHR_lifted":chr_col_new,
                                "POS_lifted":pos_col_new}, inplace=True)

        if drop_pos_build_before:
            df_merged.drop(columns=[chr_col, pos_col], inplace=True)

        merged_file = os.path.join(outd, "lifted_b{}.{}".format(build_to, outf)); log_saved_file.append("Lifted GWAS file: {}".format(merged_file))
        with open(merged_file, 'w') as f:
            row_list = df_merged.values.tolist()
            col_list = df_merged.columns.tolist()
            rows2write = [delim.join(col_list) + "\n"] + [delim.join([str(v) for v in row]) + "\n" for row in row_list]
            f.writelines(rows2write)

    # Log where the files have been saved.
    log_ = ":: liftOver OUTPUT ::\n" + "".join(["\t" + l + "\n" for l in log_saved_file]); log_list = logger(log_list, log_, verbose = verbose)

    # Save the log.
    log_file = os.path.join(outd, "liftOver.{}.log".format(outf))
    save_log(log_list=log_list,
            out=log_file
            )


if __name__ == "__main__":
    args = parse_args()

    if args.outf == "NA":
        args.outf = os.path.split(args.file)[1]
    if args.outd == "NA":
        args.outd = os.getcwd()
    
    main(file=args.file, 
        in_compression=args.in_compression, 
        delim=map_delim(args.delim),
        snp_col=args.snp_col, 
        chr_col=args.chr_col, 
        pos_col=args.pos_col,
        infer_col=args.infer_col,
        chr_col_new=args.chr_col_new,
        pos_col_new=args.pos_col_new,
        build_from=args.build_from, 
        build_to=args.build_to,
        outf=args.outf, 
        outd=args.outd,
        verbose=args.verbose,
        do_not_save_unlifted=args.do_not_save_unlifted,
        save_mapping_file=args.save_mapping_file,
        rm_intermediate_file=args.rm_intermediate_file,
        do_not_save_lifted_gwas=args.do_not_save_lifted_gwas,
        drop_pos_build_before=args.drop_pos_build_before
        )
