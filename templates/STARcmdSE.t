#!/bin/bash

#PBS -l nodes=1:ppn=$threads,mem=$mem
#PBS -q $clusterQ
#PBS -d ./

module load VARI/genomics
STAR --genomeLoad NoSharedMemory --quantMode GeneCounts --outFileNamePrefix SE_ --runMode alignReads \
    --alignIntronMin $minIntronLength --alignIntronMax $maxIntronLength --genomeDir $index --runThreadN $threads --readFilesIn read_1.1,read_2.1 
