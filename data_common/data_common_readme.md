# Data description

# GWAS summary statistics from UKB Neale lab

- variants.tsv.bgz  
    Neale lab에서 제공하는 UKB GWAS summary statistics에서 variant로 rsID를 matching 할 수 있는 파일.  
    Genome build: GRCh37
    Download: `wget https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz`


# Gene annotation

- Homo_sapiens.GRCh37.87.gtf.gz  
    Ensembl에서 제공하는 gene annotation 파일.  
    Ensembl gene ID와 gene의 transcription start site (TSS) position과 strand direction 등에 대한 정보를 얻을 수 있음.  
    Genome build: GRCh37  
    Download: `wget https://ftp.ensembl.org/pub/grch37/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz`

- modified.Homo_sapiens.GRCh37.87.gtf.gz  
    `Homo_sapiens.GRCh37.87.gtf.gz`를 modify한 파일.
    Script: ./src/modify_Homo_sapiens.GRCh_X.py  

- Homo_sapiens.GRCh38.109.gtf.gz
    Ensembl에서 제공하는 gene annotation 파일.  
    Ensembl gene ID와 gene의 transcription start site (TSS) position과 strand direction 등에 대한 정보를 얻을 수 있음.  
    Genome build: GRCh38  
    Download: `wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz`

- modified.Homo_sapiens.GRCh38.109.gtf.gz  
    `Homo_sapiens.GRCh38.109.gtf.gz`를 modify한 파일.  
    Script: ./src/modify_Homo_sapiens.GRCh_X.py  

# NCBI Entrez gene ID to gene name mapping

- gene_info.gz
    NCBI에서 제공하는 gene information 파일.  
    NCBI Entrez gene ID를 gene symbol로 mapping 할 때 사용.  
    Download: `wget https://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz`


# dbGAP Variant information mapping file

- 00_All_b37.vcf.gz
    dbGAP에서 제공하는 SNP 정보 파일. (human_9606_b151_GRCh37p13, 2018-04-27 release)
    SNP rsID만 있을 때, chromosome, base position, reference allele, alternative allele mapping 할 때 사용하면 됨.
    이건 GRCh37.
    Downlaod: `wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz -O 00_All_b37.vcf.gz`
    md5: `wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz.md5 -O 00_All_b37.vcf.gz.md5`

- 00_All_b38.vcf.gz
    dbGAP에서 제공하는 SNP 정보 파일. (human_9606_b151_GRCh38p7, 2018-04-24 release)
    SNP rsID만 있을 때, chromosome, base position, reference allele, alternative allele mapping 할 때 사용하면 됨.
    이건 GRCh38.
    Downlaod: `wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz -O 00_All_b38.vcf.gz`
    md5: `wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz.md5 -O 00_All_b38.vcf.gz.md5`