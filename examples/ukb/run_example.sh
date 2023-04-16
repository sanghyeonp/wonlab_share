#!/bin/bash

Rscript ../../src/ukb/ukb_field_extract.R --field 129 130 \
                                        --retain-all-instances \
                                        --outf example.out.tsv \
                                        --delim-out tab

# Leave only few lines as an example due to storage limit
head -30 example.out.tsv > example.out.tsv.head
rm example.out.tsv
