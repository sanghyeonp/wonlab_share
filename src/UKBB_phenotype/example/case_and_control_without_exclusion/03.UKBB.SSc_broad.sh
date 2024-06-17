#!/bin/bash
#SBATCH -J 03.UKBB.SSc_broad
#SBATCH -p cpu
#SBATCH -o ./%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=sh.austin.park@gmail.com
#SBATCH --mail-type=END,FAIL

module purge

CONDA_PATH=/data1/software/anaconda3 
source $CONDA_PATH/bin/activate /home/sanghyeon/.conda/envs/R4_HPC 

/home/sanghyeon/.conda/envs/R4_HPC/bin/Rscript 03.UKBB.SSc_broad.R