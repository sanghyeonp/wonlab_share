# Cohen's kappa coefficient for effect estimate direction concordance between GWAS

본 repository는 두 GWAS간의 effect direction concordance를 Cohen's kappa test라는 statistical test를 통해 통계적으로 concordant 한지 아닌지를 확인하는 script입니다.

저희 연구실에서는 김소연 박사님 해당 방법을 활용한 논문 (Figure 2)이 publish 되어있습니다 (Kim, 2022).

Cohen's kappa test는 R package *fmsb*의 *Kappa.test*를 활용하였습니다.

추가 되었으면 하는 기능 및 설명이 있으시면 언제든지 알려주세요!

---

### **To-do**
- 여러 GWAS pair에 대한 test 할 수 있게.
- table (csv) 결과 설명.

### **Version control**
*Version 1:*
- Cohen's kappa test를 따로 function을 만들어서, 나중에 multi-pair script에서 사용할 수 있게.

---

## How to use?
1. `git clone`을 활용해 관련 script를 다운로드 하기.
```
git clone https://github.com/sanghyeonp/cohenkappa.git
```

2. `cohen_kappa.R`를 이용하여 Cohen's kappa test 진행하기. (**An absolute path of the script MUST be passed.**)
```
Rscript <Full directory path>/cohen_kappa.R --ref_gwas <GWAS summary statistics 1> \
                    --alt_gwas <GWAS summary statistics 2>
```
&nbsp;

<ins>**아래는 현재 구현되어있는 arguments에 대한 설명입니다.**</ins>

- 아래 arguments들은 필수 입니다.  
`--ref_gwas` :  Path to GWAS summary statistics 1.  
`--alt_gwas` :  Path to GWAS summary statistics 2.  

- 아래 arguments들은 summary statistics와 관련된 **optional arguments** 입니다.  
`--ref_gwas_delim` : GWAS summary statistics 1의 delimiter. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.  
`--ref_gwas_snp` : GWAS summary statistics 1의 SNP column 이름. Default = 'SNP'.  
`--ref_gwas_beta` : GWAS summary statistics 1의 SNP column 이름. Default = 'Beta'.  
`--ref_name` : GWAS summary statistics 1의 trait 이름. Default = Summary statistics 파일 이름.  
`--alt_gwas_delim` : GWAS summary statistics 2의 delimiter. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.  
`--alt_gwas_snp` : GWAS summary statistics 2의 SNP column 이름. Default = 'SNP'.  
`--alt_gwas_beta` : GWAS summary statistics 2의 SNP column 이름. Default = 'Beta'.  
`--alt_name` : GWAS summary statistics 2의 trait 이름. Default = Summary statistics 파일 이름.  
`--alt_rev_beta` : Specify to reverse the effect direction of GWAS summary statistics 2. Default = FALSE.  

- 아래 arguments들은 특정 SNP만 Cohen's kappa test에 포함시키기 위한 **optional arguments** 입니다. 만약 특정 SNP을 specify하지 않으면, GWAS summary statistics 1에 포함된 모든 SNP들을 기준으로 test가 진행됩니다.  
`--snplist` : Path to file with SNP list. Default = 'NA'.  
`--snplist_FUMA` : FUMA에서 얻은 `leadSNPs.txt` 결과 path. Default = 'NA'.  

- 아래 arguments들은 결과 저장과 관련된 **optional arguments** 입니다.  
`--outf` : Name of the output files. Default='cohenkappa'.  
`--outd` : Path to output directory. Default=Current working directory.  

- Additional **optional arguments**.  
`--verbose` : (flag) 결과를 terminal에 print. Default = FALSE.  
`--rds` : (flag) *Kappa.test* 결과를 RDS 파일로 저장. Default = FALSE.  
`--table` : (flag) 결과를 table로 정리해서 csv 파일로 저장. Default = FALSE.  

---

## Example

`/example` 에 있는 예제를 참고하세요.

용량 제한으로 GWAS1과 GWAS2 summary statistics는 subset만 공유드립니다.

`/example/run_example1.sh` 은 `--snplist` 를 이용한 예제입니다.  
`/example/run_example2.sh` 는 `--snplist_FUMA` 를 이용한 예제입니다.


---
## Contributors


---
## Reference

1. Kim, S., Kim, K., Hwang, M. Y., Ko, H., Jung, S.-H., Shim, I., Cha, S., Lee, H., Kim, B., Yoon, J., Ha, T. H., Kim, D. K., Kim, J., Park, W.-Y., Okbay, A., Kim, B.-J., Kim, Y. J., Myung, W., & Won, H.-H. (2022). Shared genetic architectures of subjective well-being in East Asian and European ancestry populations. Nature Human Behaviour, 6(7), Article 7. https://doi.org/10.1038/s41562-022-01343-5
