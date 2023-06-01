#!/bin/bash

file=/data1/sanghyeon/Projects/MetabolicSyndrome/GWASsumstat/QC_final/CLEANED.TG_GLGC_UKB.txt
awk '{ print $1, $2 }' ${file} | shuf -n 1000 > random_1k.txt