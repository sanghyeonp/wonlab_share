#!/bin/bash
#SBATCH -J run_example
#SBATCH -p cpu
#SBATCH -o ./%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=sh.austin.park@gmail.com   # 알림 메일 보낼 주소
#SBATCH --mail-type=END,FAIL              # 알림 메일 타입 (END, FAIL, REQUEUE, TIME_LIMIT_80)


### 개인 conda 환경 쓸 때
module purge

CONDA_PATH=/data1/software/anaconda3         # 콘다 실행 경로
source $CONDA_PATH/bin/activate /home/sanghyeon/.conda/envs/R4_HPC    # 콘다 환경 실행


script=/data1/sanghyeon/wonlab_contribute/combined/src/liftover_auto/liftover_auto.R

#####################
# Example 1: Basic
#####################
/home/sanghyeon/.conda/envs/R4_HPC/bin/Rscript ${script} \
    --file-in ./data/finngen_R8_T2D.autosome.10k.subset \
    --delim tab \
    --snp-col rsids \
    --chr-col "#chrom" \
    --pos-col pos \
    --build-from 38 \
    --build-to 37 \
    --out-pref example1_out \
    --save-mapping-file



# ### Chromosome, position, SNP column이 따로 있는 경우.
# python3 ${script} \
#     --file ./data/finngen_R8_T2D.autosome.10k.subset \
#     --delim tab \
#     --snp-col rsids \
#     --chr-col "#chrom" \
#     --pos-col pos \
#     --build-from 38 \
#     --build-to 37 \
#     --outf example1_out \
#     --verbose


