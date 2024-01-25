#!/bin/bash
#SBATCH -J 01_run_independent_from_gwas
#SBATCH -p cpu
#SBATCH -o ./%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=31
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=sh.austin.park@gmail.com   # 알림 메일 보낼 주소
#SBATCH --mail-type=END,FAIL              # 알림 메일 타입 (END, FAIL, REQUEUE, TIME_LIMIT_80)


### 개인 conda 환경 쓸 때
module purge

CONDA_PATH=/data1/software/anaconda3         # 콘다 실행 경로
source $CONDA_PATH/bin/activate /home/sanghyeon/.conda/envs/R4_HPC    # 콘다 환경 실행

script_path=/data1/sanghyeon/wonlab_contribute/combined/src/independent_snp/snp_independence.exe.R

/home/sanghyeon/.conda/envs/R4_HPC/bin/Rscript ${script_path} \
    --gwas1-snplist ./data/snplist.txt \
    --delim-gwas1-snplist whitespace \
    --snp-col-gwas1-snplist SNP \
    --chr-col-gwas1-snplist Chr \
    --pos-col-gwas1-snplist bp \
    --gwas2 ./data/gwas.tsv \
    --delim-gwas2 tab \
    --snp-col-gwas2 SNP \
    --chr-col-gwas2 CHR \
    --pos-col-gwas2 POS \
    --pval-col-gwas2 Pval \
    --reference-panel "1kG" \
    --r2-threshold 0.1 \
    --window 500 \
    --pval-threshold 5e-8 \
    --n-thread 30 \
    --prefix-out SNP_independence.FG_SNPlist.HTN_GWAS \
    --delim-out comma
