#!/bin/bash

for dir in Sample_?;
do
    echo $dir;
    cd $dir;
    $(echo -e "EMBLID\t${dir}" > gene_cnts.txt)
    $(grep -v "N_" merged_ReadsPerGene.out.tab | cut -f 1,2 >> gene_cnts.txt)
    cd ..

done
