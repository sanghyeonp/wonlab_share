import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from leadSNPannotation import leadSNPannotation
from leadSNP2gene import leadSNP2GeneMapping
from magma_genebased import MAGMA_genebased2Table
from magma_geneset import MAGMA_geneset2Table
from magma_geneproperty import MAGMA_geneproperty_GTEx_SpecificTissue, MAGMA_geneproperty_GTEx_GeneralTissue, MAGMA_geneproperty_BrainSpan_Age, MAGMA_geneproperty_BrainSpan_Developmental

import argparse


def parse_args():
    parser = argparse.ArgumentParser(description=":: Make a neat summarized results of FUMA and MAGMA ::")

    # Data directory
    parser.add_argument('--result_dir', required=True,
                        help='Specify path to directory where the necessary result files are present.')
    
    # Run all
    parser.add_argument('--run_all', action='store_true',
                    help='Specify to make tables for all the available results. Default = False.')
    
    # SNP annotation
    parser.add_argument('--snp_annotation', action='store_true',
                    help='Specify to make SNP annotation result table. Default = False.')
    # SNP2GENE
    parser.add_argument('--snp2gene_mapping', action='store_true',
                    help='Specify to make SNP2Gene mapping result table. Default = False.')
    # MAGMA gene-based analysis
    parser.add_argument('--magma_genebased', action='store_true',
                    help='Specify to make MAGMA gene-based analysis result table. Default = False.')
    # MAGMA gene-set analysis
    parser.add_argument('--magma_geneset', action='store_true',
                    help='Specify to make MAGMA gene-set analysis result table. Default = False.')
    # MAGMA gene-property analysis
    parser.add_argument('--magma_geneproperty_gtex_specific', action='store_true',
                    help='Specify to make MAGMA gene-property analysis in GTEx v8 specific tissues result table. Default = False.')
    parser.add_argument('--magma_geneproperty_gtex_general', action='store_true',
                    help='Specify to make MAGMA gene-property analysis in GTEx v8 general tissues result table. Default = False.')
    parser.add_argument('--magma_geneproperty_brainspan_age', action='store_true',
                    help='Specify to make MAGMA gene-property analysis in BrainSpan age result table. Default = False.')
    parser.add_argument('--magma_geneproperty_brainspan_dev', action='store_true',
                    help='Specify to make MAGMA gene-property analysis in BrainSpan developmental stages result table. Default = False.')
    
    # Output.
    parser.add_argument('--outd', required=False, default="NA",
                        help='Specify the output directory path. Default = Current working directory.')

    # Others
    parser.add_argument('--verbose', action='store_true',
                        help='Specify to print logs.')
    
    args = parser.parse_args()
    return args



if __name__ == "__main__":
    args = parse_args()

    if args.outd == "NA":
        args.outd = os.getcwd()

    if args.run_all:
        args.snp_annotation = True
        args.snp2gene_mapping = True
        args.magma_genebased = True
        args.magma_geneset = True
        args.magma_geneproperty_gtex_specific = True
        args.magma_geneproperty_gtex_general = True
        args.magma_geneproperty_brainspan_age = True
        args.magma_geneproperty_brainspan_dev = True

    if args.snp_annotation:
        leadSNPannotation(result_dir=args.result_dir, 
                        verbose=args.verbose,
                        outd=args.outd
                        )
    
    if args.snp2gene_mapping:
        leadSNP2GeneMapping(result_dir=args.result_dir, 
                            verbose=args.verbose,
                            outd=args.outd
                            )
    
    if args.magma_genebased:
        MAGMA_genebased2Table(result_dir=args.result_dir, 
                            snp2gene_table=None, 
                            verbose=args.verbose,
                            outd=args.outd
                            )
        
    if args.magma_geneset:
        MAGMA_geneset2Table(result_dir=args.result_dir, 
                            verbose=args.verbose,
                            outd=args.outd
                            )
    
    if args.magma_geneproperty_gtex_specific:
        MAGMA_geneproperty_GTEx_SpecificTissue(result_dir=args.result_dir, 
                                                verbose=args.verbose,
                                                outd=args.outd
                                                )
    if args.magma_geneproperty_gtex_general:
        MAGMA_geneproperty_GTEx_GeneralTissue(result_dir=args.result_dir, 
                                                verbose=args.verbose,
                                                outd=args.outd
                                                )
    if args.magma_geneproperty_brainspan_age:
        MAGMA_geneproperty_BrainSpan_Age(result_dir=args.result_dir, 
                                                verbose=args.verbose,
                                                outd=args.outd
                                                )
    if args.magma_geneproperty_brainspan_dev:
        MAGMA_geneproperty_BrainSpan_Developmental(result_dir=args.result_dir, 
                                                verbose=args.verbose,
                                                outd=args.outd
                                                )

