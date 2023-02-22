# FUMA results to table

안녕하세요? 박상현입니다.

GWAS 연구를 수행하면서 FUMA라는 web-based software를 사용해서 SNP annotation, SNP to gene mapping, 다양한 MAGMA 분석 등을 진행하게 됩니다.  

다양한 분석이 진행되고, 결과들을 다운로드 받으면, 어디서부터 어떤 결과를 보고 해석해야 할지 막막할 수 있는데요.  

그래서 보다 FUMA 결과를 1차로 쉽게 정리할 수 있게, 코드 공유드립니다.

(제 study note 이기도 합니다. ㅎㅎ)

### To-do
- 결과를 다른 directory에 저장할 수 있게. (현재는 current working directory)
- 사용된 results 설명
- 정리된 table 설명

### Version control
Version 1:
- SNP annotation 결과 정리
- SNP-to-gene 결과 정리
- MAGMA gene-based analysis 결과 정리
- MAGMA gene-set analysis 결과 정리
- MAGMA gene-property analysis 결과 정리

## How to use?
1. `git clone`을 활용해 관련 script를 다운로드 하기.
```
git clone https://github.com/sanghyeonp/FUMAresult2table.git
```

2. `main.py`를 이용하여 table 만들기.
```
python3 main.py --result_dir <FUMA 결과 있는 directory> --run_all
```

## Example
Lagou et al,에서 다운로드 받은 fasting glucose GWAS를 예제 삼아 진행해보았습니다.



## Reference
1. Lagou, V. et al. Sex-dimorphic genetic effects and novel loci for fasting glucose and insulin variability. Nat Commun 12, 24 (2021).