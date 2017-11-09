#!/bin/bash

for dir
do
    echo $dir
    cd $dir/merged
    sbatch -J $dir --wrap="count_mapped_reads_bam.sh merged.bam"
    cd ../../
done


