import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *
from make_annovar_input import read_input, make_annovar_input
from annovar import annovar
from map_annotation import map_annovar_out


def parse_args():
    parser = argparse.ArgumentParser(description=":: Make ANNOVAR input file from user-specified input file ::")
    parser.add_argument('--file', required=True,
                        help='Path to the input file.')

    parser.add_argument('--in_compression', required=False, default="infer",
                        help="Specify compression type from the following ['zip', 'gzip', 'bz2', 'zstd', 'tar']. Default='infer'.")
    parser.add_argument('--delim_in', required=False, default="tab",
                    help="Delimiter used in the input file. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.")
    
    parser.add_argument('--infer_chr_pos_ref_alt', nargs='+', required=False,
                        help="Specify column name, data format, and separator to infer chr and pos from specified column.\
                            For example, column named 'variant' have variant name '2:179816990:C:T' where chromosome and position can be inferred as 2 and 179816990, respectively.\
                            Then, specify as follows: --infer_chr_pos_ref_alt variant CHR:POS:REF:ALT :")
    parser.add_argument('--chr_col', required=False, default="CHR",
                    help="Name of the chromosome column in the input file. Default = 'CHR'.")
    parser.add_argument('--pos_col', required=False, default="POS",
                    help="Name of the base position column in the input file. Default = 'POS'.")
    parser.add_argument('--ref_col', required=False, default="REF",
                    help="Name of the reference allele column in the input file. Default = 'REF'.")
    parser.add_argument('--alt_col', required=False, default="ALT",
                    help="Name of the alternative allele column in the input file. Default = 'ALT'.")
    
    parser.add_argument('--outf', required=False, default="NA",
                        help="Specify the name of the output file. Default = 'rsidmapped.<file>'.")
    parser.add_argument('--outd', required=False, default="NA",
                        help="Specify the path to output directory. Default = Current working directory.")
    parser.add_argument('--delim_out', required=False, default="NA",
                        help="Delimiter for the output file. Choices = ['NA', 'tab', 'comma', 'whitespace']. If 'NA', identical delimiter as input delimiter will be used. Default = 'NA'.")
    parser.add_argument('--out_compression', required=False, default="NA",
                        help="Specify compression type from the following ['NA', 'zip', 'gzip', 'bz2', 'zstd', 'tar']. Default='NA'.")

    parser.add_argument('--do_not_save_annotated_input', action='store_true',
                    help='Specify to avoid saving annotated input file. Default = False.')
    parser.add_argument('--save_unannotated_snp', action='store_true',
                    help='Specify to save un-annotated SNP list as chr:pos:ref:alt format with following file name: "unannotated_variant.list". Default = False.')
    parser.add_argument('--save_flipped_snp', action='store_true',
                    help='Specify to save flipped SNP list with columns chr:pos:ref:alt and chr:pos:ref:alt_new and with following file name: "flipped_variant.list". Default = False.')
    parser.add_argument('--save_mapping_file', action='store_true',
                    help='Specify to save the SNP mapping file that can be used to map GWAS summary statistics from identical cohort. \
                        Columns are ["chr_pos_ref_alt", "chr_pos_ref_alt_new", "rsid"] and tab-delimited. \
                        Saved with the following file name: "mapping.file". Default = False.')
    parser.add_argument('--delete_intermediate_files', action='store_true',
                    help='Specify to delete all the intermediate files generated during annotation. Default = False.')

    args = parser.parse_args()
    return args


def main(file, delim_in, in_compression,
        infer_chr_pos_ref_alt, chr_col, pos_col, ref_col, alt_col,
        outf, outd, delim_out, out_compression,
        save_unannotated_snp, save_flipped_snp, save_mapping_file, do_not_save_annotated_input,
        delete_intermediate_files):
    
    # Read input file
    input_df_ = read_input(file, delim_in, in_compression, infer_chr_pos_ref_alt, chr_col, pos_col, ref_col, alt_col)

    # Make ANNOVAR input file
    annov_input1, annov_input2 = make_annovar_input(input_df_, file, chr_col, pos_col, ref_col, alt_col)

    # Run ANNOVAR
    annovar(annov_input1, annov_input2)

    # Map ANNOVAR result
    df_out_ = map_annovar_out(input_df_, annov_input1, annov_input2, chr_col, pos_col, ref_col, alt_col)

    ## Save necessary data
    n_tot = len(df_out_)
    print("Number of SNPs in the input file: {:,}".format(n_tot))
    
    # Annotated.
    df_out = df_out_[~df_out_['rsID_annov'].isna()]
    print("Number of SNPs annotated: {:,} ({:.1%})".format(len(df_out), len(df_out)/n_tot))
    
    ## Un-annotated.
    # Not annotated.
    df_out_no_annot = df_out_[df_out_['rsID_annov'].isna()]
    print("Number of SNPs not annotated: {:,} ({:.1%})".format(len(df_out_no_annot), len(df_out_no_annot)/n_tot))

    # Save SNPs that were not annotated
    if save_unannotated_snp:
        df_out_no_annot[['chr_pos_ref_alt']].to_csv("unannotated_variant.list", sep="\t", index=False, header=False)

    ## Flipped.
    if not df_out[df_out['flipped_annov'] == "1"].empty:
        print("Number of SNPs flipped: {:,}".format(len(df_out[df_out['flipped_annov'] == "1"])))
        if save_flipped_snp:
            df_out[df_out['flipped_annov'] == "1"][['chr_pos_ref_alt', 'chr_pos_ref_alt_new']].to_csv("flipped_variant.list", sep="\t", index=False)
    else:
        print("Number of SNPs flipped: 0")

    # Save mapping file
    if save_mapping_file:
        df_mapping = df_out[['chr_pos_ref_alt', 'chr_pos_ref_alt_new', 'rsID_annov']]
        df_mapping.columns = ['chr_pos_ref_alt', 'chr_pos_ref_alt_new', 'rsid']
        df_mapping.to_csv("mapping.file", sep="\t", index=False)

    # Save output file with annotated input file
    if not do_not_save_annotated_input:
        df_out.drop(columns=['chr_pos_ref_alt', chr_col, pos_col, ref_col, alt_col, 'flipped_annov'], inplace=True)
        df_out[[chr_col, pos_col, 'REF', 'ALT']] = df_out.apply(lambda row: row['chr_pos_ref_alt_new'].split(sep=":"), axis=1, result_type='expand')
        df_out.drop(columns=['chr_pos_ref_alt_new'], inplace=True)
        df_out.to_csv(os.path.join(outd, outf), sep=delim_out, index=False, compression=out_compression)

    # Remove intermediate files
    if delete_intermediate_files:
        input1_prefix = annov_input1.replace(".annovin", "")
        input2_prefix = annov_input2.replace(".annovin", "")
        cmd = "rm {}.annovin {}.annovin {}.hg19_multianno.txt {}.hg19_multianno.txt {}.invalid_input {}.invalid_input {}.refGene.invalid_input {}.refGene.invalid_input".format(
            input1_prefix, input2_prefix, input1_prefix, input2_prefix, input1_prefix, input2_prefix, input1_prefix, input2_prefix
        )
        run_bash(cmd)


if __name__ == "__main__":
    args = parse_args()

    if args.outf == "NA":
        args.outf = "annovarmapped.{}".format(os.path.split(args.file)[-1])

    if args.outd == "NA":
        args.outd = os.getcwd()
    
    if args.delim_out == "NA":
        args.delim_out = args.delim_in

    if args.out_compression == "NA":
        args.out_compression = None

    main(file=args.file, 
        delim_in=map_delim(args.delim_in), 
        in_compression=args.in_compression,
        infer_chr_pos_ref_alt=args.infer_chr_pos_ref_alt, 
        chr_col=args.chr_col, 
        pos_col=args.pos_col, 
        ref_col=args.ref_col, 
        alt_col=args.alt_col,
        outf=args.outf, 
        outd=args.outd, 
        delim_out=map_delim(args.delim_out), 
        out_compression=args.out_compression,
        save_unannotated_snp=args.save_unannotated_snp, 
        save_flipped_snp=args.save_flipped_snp, 
        save_mapping_file=args.save_mapping_file, 
        do_not_save_annotated_input=args.do_not_save_annotated_input,
        delete_intermediate_files=args.delete_intermediate_files
        )