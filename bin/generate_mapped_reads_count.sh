#!/bin/bash

for dir
do
    echo $dir
    cd ${dir}/merged
    sbatch -J $dir -o %J.o -e %J.e --wrap="count_mapped_reads_bam.sh merged.bam"
    cd ../..
done


