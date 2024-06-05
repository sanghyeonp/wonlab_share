# Allele frequency (AF) calculation if AFcase and AFcontrol given separately

Allele 1에 대해서 allele frequency를 줬는데, case and control 각각 나눠서 주었다면, allele 1에 대한 total AF는 계산을 따로 해야함.

이때, case and control에 대한 sample size 정보도 있어야 함.

$AF_{allele 1} = ((AF_{case} * N_{case})  + (AF_{control} * N_{control})) / (N_{case} + N_{control})$

- Reference  
  - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3480678/
