#!/bin/bash
#SBATCH -J 01.query_UKBB.41271.ICD9_main
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

# 41271: Diagnoses - ICD9

/home/sanghyeon/.conda/envs/R4_HPC/bin/Rscript ${script} \
    --field 41271 \
    --retain-all-instances \
    --outf UKBB.41271.diagnoses_main_ICD9.rds \
    --delim-out rds \
    --n-cores 20