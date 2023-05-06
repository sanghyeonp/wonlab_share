#!/bin/bash

### Ensembl gene ID to gene symbol, gene chromosome, gene TSS, and gene strand direction (mapping)
# script=/data1/sanghyeon/wonlab_contribute/combined/src/gene_annotation/Ensembl_gene_annotation.py
script=../../src/gene_annotation/Ensembl_gene_annotation.py

python3 ${script} --file ./data/example_data_ensembl.tsv \
                --gene-id-col Probe_Ensembl \
                --delim-in tab \
                --outf example_data_Ensembl_gene_annotation.tsv \
                --delim-out tab
