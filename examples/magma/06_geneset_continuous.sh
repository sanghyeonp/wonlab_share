#!/bin/bash

magma=../../software/MAGMA/magma

genes_raw=AN_genebased.genes.raw
covar_file=./geneset/tissue_gex.cov

${magma} --gene-results ${genes_raw} --gene-covar ${covar_file} --out AN_geneset_cont
