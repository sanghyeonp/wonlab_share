#!/bin/bash

# Input file containing SNPs for all chromosomes
input_file=/data1/sanghyeon/Projects/psych_gwasbysub/src/03_modify_EasyQC_out/CLEANED.als2023.2.txt

output_file=MDD_als2023.random10k.txt

# Extract and save the header
head -n 1 "$input_file" > "$output_file"

# Loop through chromosomes 1 to 22
for chromosome in {1..22}; do
    # Generate a random sample of 1000 SNPs for the current chromosome
    awk -v chromosome="$chromosome" 'BEGIN { srand() } $1 == chromosome { print $0 }' "$input_file" | shuf -n 10000 >> "$output_file"
done
