#!/bin/bash

# script=/data1/sanghyeon/wonlab_contribute/combined/src/compute_ld/compute_ld_single.py

######################################################
# SNP in rsID format
# European population 1000 Genome reference
# Both SNPs available in the reference
# 결과: 0.402158
######################################################
python3 ../../src/compute_ld/compute_ld_single.py \
    --snp1 rs2840528 \
    --snp2 rs7545940 \
    --reference 1kg \
    --variant-identifier rsid \
    --ancestry European

######################################################
# SNP in chr:pos format
# European population 1000 Genome reference
# Both SNPs available in the reference
# 결과: 0.402158
######################################################
python3 ../../src/compute_ld/compute_ld_single.py \
    --snp1 1:2283896 \
    --snp2 1:2299627 \
    --reference 1kg \
    --variant-identifier chr:pos \
    --ancestry European

######################################################
# SNP in chr:pos format
# East-asian population 1000 Genome reference
# Both SNPs available in the reference
# 결과: 0.419706
######################################################
python3 ../../src/compute_ld/compute_ld_single.py \
    --snp1 1:2283896 \
    --snp2 1:2299627 \
    --reference 1kg \
    --variant-identifier chr:pos \
    --ancestry East-asian


######################################################
# SNP in rsID format
# European population 1000 Genome reference
# SNP2 not available in the reference
# 결과: -9
######################################################
python3 ../../src/compute_ld/compute_ld_single.py \
    --snp1 rs2840528 \
    --snp2 rs754232320 \
    --reference 1kg \
    --variant-identifier rsid \
    --ancestry European