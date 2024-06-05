#!/bin/bash
#SBATCH -J example
#SBATCH -p cpu
#SBATCH -o ./%x_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=25
#SBATCH --time=14-00:00:00
#SBATCH --mail-user=sh.austin.park@gmail.com
#SBATCH --mail-type=END,FAIL

module purge

CONDA_PATH=/data1/software/anaconda3 
source $CONDA_PATH/bin/activate /home/sanghyeon/.conda/envs/R4_HPC 

src=gwas_basic_QC_check.R

/home/sanghyeon/.conda/envs/R4_HPC/bin/Rscript ${src} \
    --gwas /data1/sanghyeon/Projects/atherosclerosis_gwas/data/gwas/AAA/AAA_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz \
    --snp-col rsid \
    --chr-col '#CHR' \
    --pos-col POS \
    --a1-col REF \
    --a2-col ALT \
    --af-col all_meta_AF \
    --effect-col inv_var_meta_beta \
    --se-col inv_var_meta_sebeta \
    --additional-criteria 'inv_var_het_p<=0.05;is_diff_AF_gnomAD!="no"' \
    --additional-col 'inv_var_het_p;is_diff_AF_gnomAD' \
    --pref example \
    --save-snp-count 'csv' \
    --save-snp-list 'rds' \
    --n-thread 20