#!/bin/bash

### NCBI Entrez gene ID to gene symbol (mapping)
script=../../src/gene_annotation/NCBI_Entrez_ID_to_symbol.py

python3 ${script} --file ./data/example_data_NCBI_Entrez.tsv \
                --gene_id_col GENE \
                --delim formatted \
                --outf example_data_NCBI_Entrez_genesymbol.tsv
