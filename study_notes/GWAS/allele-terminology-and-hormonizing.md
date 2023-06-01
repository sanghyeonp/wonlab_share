# Allele terminology에 대해서
GWAS를 하다보면 비슷한 의미를 가진 것 같은 여러 terminology 때문에 머리가 복잡해지고, 이러한 terminology들이 정확하게 정의가 되지 않는다면, 여러 GWAS summary statistics를 다룰 때, 문제가 생기게 됩니다. 

특히 GWAS summary statistics를 공유할 때, format에 대한 정해진 guideline이 없기 때문에 연구자들이 공유하는 데이터마다 column의 정보가 다를 수 있고, 또 시간이 좀 지난 데이터의 경우에도 column이 의미하는게 다를 수 있습니다.

GWAS summary statistics에 report되는 것들 중, 특히 allele의 종류와 관련한 terminology에 대해서 간략하게 정리해보았습니다.  

```
도대체 아래 allele 이름들 차이가 뭐냐?!
Reference allele, Alternative allele  
Major allele, Minor allele  
Effect allele, Other allele  
```

**1. Reference allele vs Alternative allele**

- Reference allele 이란?
    - "base found in the reference genome"
- Alternative allele 이란?
    - "any base, other than the reference, that is found at that locus"  

```
Reference allele과 alternative allele의 차이는 단순히 해당 SNP의 allele이 reference genome에 있는지 없는지로 나뉘는 것 같습니다.
```

Source: https://www.ebi.ac.uk/training/online/courses/human-genetic-variation-introduction/what-is-genetic-variation/what-are-variants-alleles-and-haplotypes/#:~:text=When%20working%20with%20genome%20scale,is%20found%20at%20that%20locus.


**2. Major allele vs Minor allele**

- Major allele 이란?
    - Population내에서 allele의 frequency가 더 많은 allele.
- Minor allele 이란?
    - Population 내에서 allele의 frequency가 더 적은 allele.

```
Major allele과 minor allele을 구분 짓는 것은 해당 allele이 population내에서 얼마만큼의 빈도 (frequency)를 보이는지에 따라 나뉘는 것 같습니다.
```
```
그렇다면 Reference allele이 major allele이랑 같은 것이냐?
>>> 아님!
>>> Reference allele은 SNP의 2가지 allele 중에 어떤 allele이 reference genome에 있었는지에 따라, reference genome에 있었던 allele이 reference allele 이라고 부르게 되고, reference genome에 있지 않았던 다른 allele을 alternative allele 이라고 부르게 됨.
>>> 반면 major allele은 reference genome에 있던 allele이든 아니든 상관없이, population내에서 allele frequency가 더 많은 allele을 major allele이라고 부름. 그렇기 때문에 reference allele이지만 population 내에서 해당 allele의 frequency가 적으면, minor allele임.
```

Source : https://www.biostars.org/p/310841/

**3. Effect allele vs Other allele**

- Effect allele 이란?
    - GWAS 분석을 할 때, phenotype과 association을 확인한 allele을 effect allele이라고 부른다고 합니다.
- Other allele 이란?
    - Effect allele의 반대.

```
GWAS 분석을 할 때, 어떤 allele을 기준으로 SNP과 phenotype 간의 association을 확인했는지에 따라, association을 확인한 allele을 effect allele 이라 부르는 것 같습니다.
```
```
Effect allele이 minor allele 과 같은 의미인가?
>>> 그렇게 생각할 수 있음.
>>> 왜냐하면 GWAS 분석으로 SNP과 phenotype간의 association을 확일할 때, 기준이 되는 allele은 less frequent allele (즉 minor allele)을 이용함. 그리고 association을 확인한 allele을 effect allele이라고 부르기 때문에, minor allele이 effect allele과 같은 의미로 생각할 수 있음.
>>> 그래서 minor allele frequency (MAF)가 effect allele의 frequency를 의미한다고 생각 할 수 있음.
```
```
Effect allele이 alternative allele과 같은 의미인가?
>>> 아님.
>>> Alternative allele이 association 분석의 주가 되었다면 effect allele이라고 부를 수 있음. 하지만 항상 그렇지는 않기 때문에 단순히 alternative allele을 effect allele이라고 생각하면 안됨.
```

Source : https://www.biostars.org/p/310841/

---

# Harmonization

### GWAS summary statistics가 있을 때, minor allele frequency (MAF)가 0.5 이상인 경우.
