#!/bin/bash
# script=/data1/sanghyeon/wonlab_contribute/combined/src/parse_pops_out/parse_pops_result.py
script=../../src/parse_pops_out/parse_pops_result.py

dir_pops_result=./data

file_leadsnp=./data/leadSNPs_scz.txt
delim=tab
snp_col=rsID
chr_col=chr
pos_col=pos
window=500000

outf=scz

python3 ${script} --dir-pops-out ${dir_pops_result} \
                --file-leadsnp ${file_leadsnp} \
                --delim-in ${delim} \
                --snp-col ${snp_col} \
                --chr-col ${chr_col} \
                --pos-col ${pos_col} \
                --window ${window} \
                --out-prefix ${outf}
