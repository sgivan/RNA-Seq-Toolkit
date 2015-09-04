#!/bin/bash

for dir
do
    echo $dir
    cd $dir
#    mkdir -p fastqc
#    $( bsub -o fastqc.o -e fastqc.e -J ${dir} -R 'span[hosts=1]' -n 4 "fastqc --threads 4 --outdir fastqc set1.fq; cd fastqc; mv set1_fastqc.html ${dir}_set1_fastqc.html; cd ..;" )
#    mv set1_fastqc.html ${dir}_set1_fastqc.html

    $( sed "s/JOBNAME/${dir}/" ~/projects/RNA-Seq-Toolkit/bin/fastqc.bs > fastqc.bs )
    bsub < fastqc.bs
    cd ..
done

