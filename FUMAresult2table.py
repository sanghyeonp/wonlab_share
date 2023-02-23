from src.util import *
from src.leadSNPannotation import leadSNPannotation
from src.leadSNP2gene import leadSNP2GeneMapping
from src.magma_genebased import MAGMA_genebased2Table
from src.magma_geneset import MAGMA_geneset2Table
from src.magma_geneproperty import MAGMA_geneproperty_GTEx_SpecificTissue, MAGMA_geneproperty_GTEx_GeneralTissue, MAGMA_geneproperty_BrainSpan_Age, MAGMA_geneproperty_BrainSpan_Developmental

import argparse


def parse_args():
    parser = argparse.ArgumentParser(description=":: Make a neat summarized results of FUMA and MAGMA ::")

    # Data directory
    parser.add_argument('--result_dir', required=True,
                        help='Specify path to directory where the necessary result files are present.')
    
    # Run all
    parser.add_argument('--run_all', action='store_true',
                    help='Specify to make tables for all the available results.')
    
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
    parser.add_argument('--magma_geneset', action='store_true',
                    help='Specify to make MAGMA gene-set analysis result table.')
    # MAGMA gene-property analysis
    parser.add_argument('--magma_geneproperty_gtex_specific', action='store_true',
                    help='Specify to make MAGMA gene-property analysis in GTEx v8 specific tissues result table.')
    parser.add_argument('--magma_geneproperty_gtex_general', action='store_true',
                    help='Specify to make MAGMA gene-property analysis in GTEx v8 general tissues result table.')
    parser.add_argument('--magma_geneproperty_brainspan_age', action='store_true',
                    help='Specify to make MAGMA gene-property analysis in BrainSpan age result table.')
    parser.add_argument('--magma_geneproperty_brainspan_dev', action='store_true',
                    help='Specify to make MAGMA gene-property analysis in BrainSpan developmental stages result table.')
    
    # Others
    parser.add_argument('--verbose', action='store_true',
                        help='Specify to print logs.')
    
    args = parser.parse_args()
    return args



if __name__ == "__main__":
    args = parse_args()

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
        
    if args.magma_geneset:
        MAGMA_geneset2Table(result_dir=args.result_dir, 
                            verbose=args.verbose
                            )
    
    if args.magma_geneproperty_gtex_specific:
        MAGMA_geneproperty_GTEx_SpecificTissue(result_dir=args.result_dir, 
                                                verbose=args.verbose
                                                )
    if args.magma_geneproperty_gtex_general:
        MAGMA_geneproperty_GTEx_GeneralTissue(result_dir=args.result_dir, 
                                                verbose=args.verbose
                                                )
    if args.magma_geneproperty_brainspan_age:
        MAGMA_geneproperty_BrainSpan_Age(result_dir=args.result_dir, 
                                                verbose=args.verbose
                                                )
    if args.magma_geneproperty_brainspan_dev:
        MAGMA_geneproperty_BrainSpan_Developmental(result_dir=args.result_dir, 
                                                verbose=args.verbose
                                                )

