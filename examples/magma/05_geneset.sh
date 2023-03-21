#!/bin/bash

magma=../../software/MAGMA/magma

genes_raw=AN_genebased.genes.raw
set_file=./geneset/c5.go.bp.v2023.1.Hs.entrez.gmt.modified

${magma} --gene-results ${genes_raw} --set-annot ${set_file} --out AN_geneset

