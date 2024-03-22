#!/bin/bash
script=/data1/sanghyeon/wonlab_contribute/combined/src/liftover_auto/liftover_auto.py

#####################
# Example 1
#####################
### Chromosome, position, SNP column이 따로 있는 경우.
python3 ${script} \
    --file ./data/finngen_R8_T2D.autosome.10k.subset \
    --delim tab \
    --snp-col rsids \
    --chr-col "#chrom" \
    --pos-col pos \
    --build-from 38 \
    --build-to 37 \
    --outf example1_out \
    --verbose


#####################
# Example 2
#####################
### Chromosome, position, SNP column이 따로 없고 infer 해야하는 경우.
# --do-not-save-unlifted: unlifted에 대한 정보는 따로 저장하지 않음.
# --save-mapping-file: mapping file 따로 저장.
# --rm-intermediate-file: intermediate file 다 지우기.
# --drop-pos-build-before: 이전 build 정보들 다 새로운 GWAS에서 지우기.
python3 ${script} \
    --file ./data/finngen_R8_T2D.autosome.10k.subset.infer \
    --delim tab \
    --infer-col variant CHR:POS:X1:X2 : CHR,POS \
    --build-from 38 \
    --build-to 37 \
    --chr-col-new CHR_b37 \
    --pos-col-new POS_b37 \
    --outf example2_out \
    --verbose \
    --do-not-save-unlifted \
    --save-mapping-file \
    --rm-intermediate-file \
    --drop-pos-build-before
