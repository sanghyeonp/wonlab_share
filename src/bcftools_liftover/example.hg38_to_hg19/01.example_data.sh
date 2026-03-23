#!/bin/bash
#SBATCH -J 01.example_data
#SBATCH -p cpu
#SBATCH -o ./%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=sh.austin.park@gmail.com 
#SBATCH --mail-type=END,FAIL 

module purge


### 개인 conda 환경 쓸 때
CONDA_PATH=/data1/software/anaconda3 
source $CONDA_PATH/bin/activate /data1/software/anaconda3/envs/R_SP_4.4 

/data1/software/anaconda3/envs/R_SP_4.4/bin/Rscript 01.example_data.R