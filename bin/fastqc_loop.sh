#!/bin/bash

for dir
do
    echo $dir
    cd $dir
#    mkdir -p fastqc
#    $( bsub -J $dir -R 'span[hosts=1]' -n 4 fastqc --threads 4 --outdir fastqc set1.fq )

    mv set1_fastqc.html ${dir}_set1_fastqc.html
    cd ..
done


