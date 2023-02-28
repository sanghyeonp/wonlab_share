#!/bin/bash

metasoft=../src/Metasoft.jar
pval_table=../src/HanEskinPvalueTable.txt
input=./data/example.txt

java -jar ${metasoft} -input ${input} -pvalue_table ${pval_table} -log example.log -output example.out
