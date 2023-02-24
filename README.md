# FUMA results to table

안녕하세요? 박상현입니다.

GWAS 연구를 수행하면서 FUMA (Watanabe, 2017)라는 web-based software를 사용해서 SNP annotation, SNP to gene mapping, 다양한 MAGMA 분석 등을 진행할 경우가 종종 있는데요.  

FUMA라는 하나의 tool로 다양한 분석이 진행되고, 따라서 여러 결과들을 다운로드 받으면, 어디서부터 어떤 결과를 중요하게 보아야하고, 어떻게 결과를 정리하고, 해석해야 할 지 막막할 수 있다고 생각합니다. (저도 많은 시행착오 겪었고, 현재도 겪고 있습니다.)  

그래서 FUMA 결과를 보다 쉽게 정리할 수 있는 script를 공유드립니다.

본 repository는 저의 FUMA 관련 study note 이기도 합니다.

앞으로 꾸준히 업데이트 할 예정이며, 추가 되었으면 하는 기능 및 설명이 있으시면 언제든지 알려주세요!

---

### **To-do**
- 결과를 다른 directory에 저장할 수 있게. (현재는 current working directory)
- 사용된 results 설명
- 정리된 table 설명

### **Version control**
*Version 1:*
- SNP annotation 결과 정리
- SNP-to-gene 결과 정리
- MAGMA gene-based analysis 결과 정리
- MAGMA gene-set analysis 결과 정리
- MAGMA gene-property analysis 결과 정리

---

## How to use?
1. `git clone`을 활용해 관련 script를 다운로드 하기.
```
git clone https://github.com/sanghyeonp/FUMAresult2table.git
```

2. `FUMAresult2table.py`를 이용하여 table 만들기.
```
python3 FUMAresult2table.py --result_dir <FUMA 결과 있는 directory> --run_all
```
&nbsp;

<ins>**아래는 현재 구현되어있는 arguments에 대한 설명입니다.**</ins>

- 만약 *Version 1*에 list된 특정 결과만 원하시면, `--run_all` 대신 아래 argument를 선택적으로 사용가능합니다.  
`--snp_annotation` : (flag) SNP annotation 결과 정리  
`--snp2gene_mapping` : (flag) SNP-to-gene 결과 정리  
`--magma_genebased` : (flag) MAGMA gene-based analysis 결과 정리  
`--magma_geneset` : (flag) MAGMA gene-set analysis 결과 정리  
`--magma_geneproperty_gtex_specific` : (flag) MAGMA gene-property analysis specific-tissue GTEx v8 결과 정리  
`--magma_geneproperty_gtex_general` : (flag) MAGMA gene-property analysis general-tissue GTEx v8 결과 정리  
`--magma_geneproperty_brainspan_age` : (flag) MAGMA gene-property analysis 29 ages BrainSpan 결과 정리  
`--magma_geneproperty_brainspan_dev` : (flag) MAGMA gene-property analysis 11 developmental stage 결과 정리  
&nbsp;

- 결과 정리가 되면, 결과 table과 log가 만들어지지만, terminal에서 real-time으로 보고싶으면, 아래 argument를 사용하면 됩니다.  
`--verbose` : (flag) 결과 정리된 summary 내용을 terminal에서 확인  

---

## Example

`/example` 에 있는 예제를 참고하세요.

Lagou, 2021 에서 다운로드 받은 fasting glucose GWAS를 예제 삼아 진행해보았습니다.

---

## 정리된 결과 output table column 설명

1. SNP annotation 결과 (`--snp_annotation`)
    - ㅇㅇㄹ

2. SNP-to-gene 결과 (`--snp2gene_mapping`)

3. MAGMA gene-based analysis 결과 (`--magma_genebased`)

---
## Contributors

---

## Reference
1. Lagou, V. et al. Sex-dimorphic genetic effects and novel loci for fasting glucose and insulin variability. Nat Commun 12, 24 (2021).

2. Watanabe, K., Taskesen, E., van Bochoven, A. & Posthuma, D. Functional mapping and annotation of genetic associations with FUMA. Nat Commun 8, 1826 (2017). https://doi.org:10.1038/s41467-017-01261-5