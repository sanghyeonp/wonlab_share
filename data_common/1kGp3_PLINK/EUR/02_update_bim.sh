#!/bin/bash
exe_plink=/home/sanghyeon/tools/plink/plink
dir_1kGp3=/data/public/1kG/phase3/plink

bed_file=${dir_1kGp3}/1kg.phase3.auto.snp.qc.chrposabc.EUR.clean.bed
fam_file=${dir_1kGp3}/1kg.phase3.auto.snp.qc.chrposabc.EUR.clean.fam

bim_file=1kg.phase3.auto.snp.qc.chrposabc.EUR.clean.rsid.new.bim
${exe_plink} --bed ${bed_file} \
            --bim ${bim_file} \
            --fam ${fam_file} \
            --make-bed \
            --out "1kg.phase3.auto.snp.qc.EUR.clean.rsid"


bim_file=1kg.phase3.auto.snp.qc.chrposabc.EUR.clean.chrpos.new.bim
${exe_plink} --bed ${bed_file} \
            --bim ${bim_file} \
            --fam ${fam_file} \
            --make-bed \
            --out "1kg.phase3.auto.snp.qc.EUR.clean.chrpos"
