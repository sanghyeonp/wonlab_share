#!/bin/bash
script=../../src/mcc/mcc.py

python3 ${script} --file id2name.AN_genebased.genes.out \
                    --p_col P \
                    --delim tab \
                    --verbose