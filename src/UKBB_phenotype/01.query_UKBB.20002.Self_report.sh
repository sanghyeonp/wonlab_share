#!/bin/bash
#SBATCH -J 01.query_UKBB.20002.Self_report
#SBATCH -p cpu
#SBATCH -o ./%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=21
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=sh.austin.park@gmail.com
#SBATCH --mail-type=END,FAIL

module purge

CONDA_PATH=/data1/software/anaconda3 
source $CONDA_PATH/bin/activate /home/sanghyeon/.conda/envs/R4_HPC 


script=/data1/sanghyeon/wonlab_contribute/combined/src/ukb/ukb_field_extract.R

# 20002: Non-cancer illness code, self-reported

/home/sanghyeon/.conda/envs/R4_HPC/bin/Rscript ${script} \
    --field 20002 \
    --retain-all-instances \
    --outf UKBB.20002.Self_report.rds \
    --delim-out rds \
    --n-cores 20