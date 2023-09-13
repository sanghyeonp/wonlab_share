#!/bin/bash

# script=/data1/sanghyeon/wonlab_contribute/combined/src/gwas_plot/manhattan_plot.R
script=../../src/gwas_plot/manhattan_plot.R

gwas=./data/MDD_als2023.random10k.txt

Rscript ${script} --gwas ${gwas} \
                --snp-col SNP \
                --chr-col CHR \
                --pos-col POS \
                --p-col PVAL \
                --chr-select 1 2 3 \
                --outf Manhattan.MDD.2
