#!/bin/bash

script=../../src/annovar/run_annovar.py

python3 ${script} --file ./data/example.1k.tsv \
                --delim_in tab \
                --infer_chr_pos_ref_alt SNPID CHR:POS:REF:ALT :
