
`run_tsmr_multiple_exposures_outcomes_gene.R`
- 같은 exposure면 clumping한 dat 공용으로 쓰게 따로 save

- 그리고 exposure iv 만들기 전엔 항상 공용으로 동일한 이름이 있는지 확인.

- 최종 결과에서 Pleiotropy 판단해서, 어떤 method를 해석하면 좋을지 확인할 수 있게.


`run_tsmr_multiple_exposures_outcomes_specifiedIV.R`
- IV 따로 추출해도, endophenotype 별로 p-value가 유의하지 않을 수 있어서, 이걸 고려해서 추가로 filtering 해주어야 함. (argument 추가하기.)