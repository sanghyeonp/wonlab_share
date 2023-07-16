# Effective sample size calculation

Binary trait GWAS의 경우 GWAS마다 case:control ratio가 다른 경우가 많음.

이때 study간 sample size를 비교할 때, case + control = N 를 비교하면 comparable 하지가 않음.

그렇기 때문에 sample size를 case:control = 50:50 으로 transform 한 effective sample size (Neff) 수치를 estimate해서 비교가 가능하다고 함.

```
Neff = 4 * v * (1 - v) * N, where v is sample prevalence (v = case / (case + control))
```

또, effective sample size는 GWAS regression coefficient 중 standard error (SE)를 이용해서 estimate 할 수 있음.

SE를 이용해 추정한 sample size는 expected sample size 라고 부르기도 함.

```
MAF 10-40% 인 SNP만 이용해서,

Neff = Mean( 1/(2 * MAF * (1 - MAF)) * SE^2 )
```

Reference:
    - https://groups.google.com/g/genomic-sem-users/c/Nre4EseOet4
    - https://www.sciencedirect.com/science/article/pii/S0006322322013166?via%3Dihub
    - https://www.sciencedirect.com/science/article/pii/S2666247722000525?via%3Dihub
