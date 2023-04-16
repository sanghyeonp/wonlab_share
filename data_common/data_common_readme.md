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


# NCBI Entrez gene ID to gene name mapping

- gene_info.gz
    NCBI에서 제공하는 gene information 파일.
    NCBI Entrez gene ID를 gene symbol로 mapping 할 때 사용.
    Download: `wget https://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz`