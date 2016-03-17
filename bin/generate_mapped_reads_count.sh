#!/bin/bash

for dir
do
    echo $dir
    cd ${dir}/merged
    bsub -J $dir -o %J.o -e %J.e count_mapped_reads_bam.sh merged.bam
    cd ../..
done


