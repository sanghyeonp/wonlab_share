#!/bin/bash
script=../../src/geneID2name/geneID2name.py

python3 ${script} --file AN_genebased.genes.out \
                    --gene_id_col GENE \
                    --delim formatted