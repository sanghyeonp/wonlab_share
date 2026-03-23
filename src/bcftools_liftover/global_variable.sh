# /data1/sanghyeon/wonlab_contribute/combined/src/bcftools_liftover/global_variable.sh

# Reference and chain files
HG38_REF="/data1/ccl/KAIST/ADNI/liftover/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
HG19_REF="/data1/ccl/KAIST/ADNI/liftover/human_g1k_v37.fasta"
CHAIN_hg38ToHg19="/data1/sanghyeon/wonlab_contribute/combined/src/bcftools_liftover/chainfile/hg38ToHg19.over.chain.gz"
CHAIN_hg19ToHg38="/data1/sanghyeon/wonlab_contribute/combined/src/bcftools_liftover/chainfile/hg19ToHg38.over.chain.gz"

# Executable file
bcftools_exe="/data1/sanghyeon/wonlab_contribute/combined/software/bcftools_1.23/bcftools/bcftools"

BCFTOOLS_PLUGINS="/data1/sanghyeon/wonlab_contribute/combined/software/bcftools_1.23/bcftools/plugins"
export BCFTOOLS_PLUGINS

# Chromosome rename map
CHRMAP="/data1/sanghyeon/wonlab_contribute/combined/src/bcftools_liftover/nochr_to_chr.map"