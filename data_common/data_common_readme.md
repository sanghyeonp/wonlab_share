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

# NCBI variant annotation

- 00_All_b37.vcf.gz
    NCBI에서 제공하는 variant information 파일.  
    Genome build: GRCh37  
    Download: `https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz -O 00_All_b37.vcf.gz`

- /00_All_b37_parquet_partitioned/part.*.parquet
    `00_All_b37.vcf.gz`에서 rsID, CHR, POS, REF, ALT dataframe을 chromosome별로 parquet format으로 저장해놓은 것.   
    Script: `/data_common/src/modify_00_All_bX.vcf.py`

- 00_All_b38.vcf.gz
    NCBI에서 제공하는 variant information 파일.  
    Genome build: GRCh38  
    Download: `https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz -O 00_All_b38.vcf.gz`

- /00_All_b38_parquet_partitioned/part.*.parquet
    `00_All_b38.vcf.gz`에서 rsID, CHR, POS, REF, ALT dataframe을 chromosome별로 parquet format으로 저장해놓은 것.   
    Script: `/data_common/src/modify_00_All_bX.vcf.py`

# Minor allele frequency (MAF) annotation

- 1000 Genome (1kG) phase 3 VCF files
    - Wonlab 위치: `/data/public/1kG/phase3/vcf`
        - `5a`와 `5b` 차이: `5b`에서 population AF 수정이 있었음.
        ```
        The super-population level allele frequency is over-estimated for chrX and chrY 
        Last update 20141017 v5
        Fixed 20150813 chrX and chrY genotype VCF files v1b version, the whole genome site VCF file v5b version now have correct super-population AF
        ```
    - Genome build: GRCh37
    - 1kG에서 제공하는 VCF 파일. 각 SNP 별로 East-Asian (EAS), European (EUR), African (AFR), American (AMR), and South-Asian (SAS)에 대한 allele frequency 정보가 INFO column에 있음.
    - Multi-allelic인 경우, ","로 AF 정보 나뉘어져 있음.

    ```
    2. Allele frequency by continental super population

    The 2504 samples in the phase3 release are from 26 populations which can be categorised into five super-populations 
    by continent (listed below).  As well as the global AF in the INFO field. We added AF for each super-population to the INFO field.

    East Asian	EAS
    South Asian	SAS
    African		AFR
    European	EUR
    American	AMR

    These allele frequences were calculated by counting the AC and AN for all the individuals from a particular super population and using that
    to calculate the AF. The info tag which represents the AFs are EAS_AF, EUR_AF, AFR_AF, AMR_AF and SAS_AF

    The super population assignment for each sample can be found in integrated_call_samples_v3.20130502.ALL.panel

    AF for multi-allelic variants are reported for each allele independently, separated by ",".
    ```
    ```
    ##INFO=<ID=EAS_AF,Number=A,Type=Float,Description="Allele frequency in the EAS populations calculated from AC and AN, in the range (0,1)">
    ##INFO=<ID=EUR_AF,Number=A,Type=Float,Description="Allele frequency in the EUR populations calculated from AC and AN, in the range (0,1)">
    ##INFO=<ID=AFR_AF,Number=A,Type=Float,Description="Allele frequency in the AFR populations calculated from AC and AN, in the range (0,1)">
    ##INFO=<ID=AMR_AF,Number=A,Type=Float,Description="Allele frequency in the AMR populations calculated from AC and AN, in the range (0,1)">
    ##INFO=<ID=SAS_AF,Number=A,Type=Float,Description="Allele frequency in the SAS populations calculated from AC and AN, in the range (0,1)">
    ```

- /1kGp3_vcf/1kgp3_chr*.tsv.gz
    - VCF에서 제공하는 #CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, 그리고 INFO에서 extract한 정보를 dataframe으로 저장함.
    - tab-delimited and gzip
