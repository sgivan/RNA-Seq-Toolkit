#!/bin/bash

for dir
do
    echo $dir
    cd ${dir}
    #bsub -J $dir -o %J.o -e %J.e count_mapped_reads_bam.sh merged.bam
    sbatch -J $dir -o %J.o -e %J.e --wrap="count_mapped_reads_bam.sh Merged_Aligned.out.bam"
    cd ../..
done


