import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *

import_or_install("gzip")
import gzip

from timeit import default_timer as timer
from datetime import timedelta
tqdm.pandas(leave=False, bar_format='{l_bar}{bar:30}{r_bar}{bar:-30b}')


def parse_args():
    parser = argparse.ArgumentParser(description=":: Convert VCF file to table ::")
    
    # Required arguments
    parser.add_argument('--file', required=True,
                        help='Path to the input file.')
    
    # Optional arguments
    parser.add_argument('--delim-in', dest="delim_in", required=False, default="tab",
                        help="Delimiter used in the input file. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.")
    parser.add_argument('--compression-in', dest="compression_in", required=False, default="NA",
                        help="Currently, only supports gzip compression. Specify from the following: ['gzip']. Default='NA'.")

    parser.add_argument('--outf', required=False, default="NA",
                        help="Specify the name of the output file. Default = 'vcf2table.<file>'.")
    parser.add_argument('--outd', required=False, default="NA",
                        help="Specify the path to output directory. Default = Current working directory.")
    parser.add_argument('--delim-out', dest="delim_out", required=False, default="NA",
                        help="Delimiter for the output file. Choices = ['NA', 'tab', 'comma', 'whitespace']. If 'NA', identical delimiter as input delimiter will be used. Default = 'NA'.")
    parser.add_argument('--compression-out', dest="compression_out", required=False, default="NA",
                        help="Currently, only supports gzip compression. Specify from the following: ['gzip']. Default='NA'.")


    args = parser.parse_args()
    return args


def main(file, delim_in, compression_in, 
        outf, outd, delim_out, compression_out
        ):
    
    print(":: Converting VCF to table ::\n\tFile: {}".format(file))

    ### Read input vcf file
    if compression_in == "gzip":
        with gzip.open(file, 'rt') as f:
            file_content = f.readlines()
    elif compression_in == "NA":
        with open(file, 'r') as f:
            file_content = f.readlines()

    ### Modify the file content
    ## Drop next line and separate by delimiter
    file_content = [line.replace("\n", "").split(sep=delim_in) for line in file_content]
    ## Drop lines starting with '##'
    file_content = [line for line in file_content if line[0][:2] != "##"]
    ## Extract and drop column name
    # col_names = file_content.pop(0)

    ### Save the table
    file_content = [delim_out.join(line) + "\n" for line in file_content]
    outp = os.path.join(outd, outf)
    with open(outp, "w") as f:
        f.writelines(file_content)
    
    if compression_out == "gzip":
        run_bash("mv {} {}".format(outp, outp.replace(".gz", "")))
        run_bash("gzip {}".format(outp.replace(".gz", "")))


if __name__ == "__main__":
    args = parse_args()

    if args.delim_out == "NA":
        args.delim_out = args.delim_in

    if args.outf == "NA":
        args.outf = "vcf2table." + os.path.split(args.file)[-1]

    if args.outd == "NA":
        args.outd = "."

    main(file=args.file,
        delim_in=map_delim(args.delim_in),
        compression_in=args.compression_in,
        outf=args.outf,
        outd=args.outd,
        delim_out=map_delim(args.delim_out),
        compression_out=args.compression_out
        )
