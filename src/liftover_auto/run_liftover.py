import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from packages import *
from util import run_bash

LIFTOVER_DIR = "/data1/sanghyeon/wonlab_contribute/combined/software/liftover" 

liftover_software = os.path.join(LIFTOVER_DIR, "liftOver")

chain_mapping = {(36, 37): os.path.join(LIFTOVER_DIR,"chainfile", "hg18ToHg19.over.chain.gz"),
                (36, 38): os.path.join(LIFTOVER_DIR,"chainfile", "hg18ToHg38.over.chain.gz"),
                (37, 38): os.path.join(LIFTOVER_DIR, "chainfile", "hg19ToHg38.over.chain.gz"),
                (38, 37): os.path.join(LIFTOVER_DIR, "chainfile", "hg38ToHg19.over.chain.gz")
                }


def run_liftover(input_bed_file, 
                build_from: int, build_to: int, 
                outd):
    global liftover_software, chain_mapping

    ### Extract file name.
    _, filename = os.path.split(input_bed_file)
    if filename[-4:] == '.bed':
        filename = filename[:-4]
    
    ### Generate liftOver execution command line.
    lifted_file = os.path.join(outd, filename+".lifted")
    unlifted_file = os.path.join(outd, filename+".unlifted")
    bash_cmd = "{} {} {} {} {}".format(liftover_software, 
                                    input_bed_file, 
                                    chain_mapping[(build_from, build_to)],
                                    lifted_file,
                                    unlifted_file
                                    )

    ### Run liftOver.
    stdout = run_bash(bash_cmd)

    return lifted_file, unlifted_file


def liftOver_log(input_bed_file, lifted_file, unlifted_file):
    with open(input_bed_file, 'r') as f:
        n_snp_input = len(f.readlines())

    with open(lifted_file, 'r') as f:
        n_snp_lifted = len(f.readlines())

    with open(unlifted_file, 'r') as f:
        n_snp_unlifted = len(f.readlines()) // 2
    
    return ":: liftOver RESULT ::\n\tNumber of SNPs initially: {:,}\n\tNumber of SNPs lifted: {:,}\n\tNumber of SNPs unlifted: {:,}".format(n_snp_input, n_snp_lifted, n_snp_unlifted)


def reformat_unlifted(unlifted_file):
    with open(unlifted_file, 'r') as f:
        rows = [row.strip().split(sep='\t') for row in f.readlines()]

    n_row = len(rows)
    idx = 1
    new_rows = []
    while idx < n_row:
        new_row = rows[idx] + rows[idx - 1]
        new_rows.append([new_row[0].replace('chr', '')] + new_row[2:])
        idx += 2
    
    return new_rows
