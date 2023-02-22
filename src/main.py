from util import *
from leadSNPannotation import leadSNPannotation
from leadSNP2gene import leadSNP2GeneMapping
from magma_genebased import MAGMA_genebased2Table

import argparse


def parse_args():
    parser = argparse.ArgumentParser(description=":: Make a neat summarized results of FUMA and MAGMA ::")

    parser.add_argument('--result_dir', required=True,
                        help='Specify path to directory where the necessary result files are present.')
    # SNP annotation
    parser.add_argument('--snp_annotation', action='store_true',
                    help='Specify to make SNP annotation result table.')
    # SNP2GENE
    parser.add_argument('--snp2gene_mapping', action='store_true',
                    help='Specify to make SNP2Gene mapping result table.')
    # MAGMA gene-based analysis
    parser.add_argument('--magma_genebased', action='store_true',
                    help='Specify to make MAGMA gene-based analysis result table.')
    # MAGMA gene-set analysis
    
    parser.add_argument('--verbose', action='store_true',
                        help='Specify to print logs.')
    
    
    args = parser.parse_args()
    return args



if __name__ == "__main__":
    args = parse_args()

    if args.snp_annotation:
        leadSNPannotation(result_dir=args.result_dir, 
                        verbose=args.verbose)
    
    if args.snp2gene_mapping:
        leadSNP2GeneMapping(result_dir=args.result_dir, 
                            verbose=args.verbose
                            )
    
    if args.magma_genebased:
        MAGMA_genebased2Table(result_dir=args.result_dir, 
                            snp2gene_table=None, 
                            verbose=args.verbose
                            )
        
    