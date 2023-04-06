#!/bin/bash

script=../../src/geneID2name/geneID2name.py

python3 ${script} --file ./data/example.data.tsv \
                --gene_id_col GENE \
                --delim formatted \
                --outf example.data.genename.tsv
