#!/bin/bash

magma=../../software/MAGMA/magma

reference=../../software/MAGMA/reference/g1000_eur
pval_file=AN.pval
annot_file=AN.genes.annot

${magma} --bfile ${reference} --pval ${pval_file} ncol=N --gene-annot ${annot_file} --out AN_genebased
