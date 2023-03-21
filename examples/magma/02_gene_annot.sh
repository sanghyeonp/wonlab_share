#!/bin/bash

magma=../../software/MAGMA/magma

snploc_file=AN.snp.loc
geneloc_file=../../software/MAGMA/aux_files/NCBI37.3/NCBI37.3.gene.loc

${magma} --annotate --snp-loc ${snploc_file} --gene-loc ${geneloc_file} --out AN
