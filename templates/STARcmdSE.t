#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$threads
#SBATCH --mem=$mem

STAR --genomeLoad NoSharedMemory --quantMode GeneCounts --outFileNamePrefix SE_ --runMode alignReads \
    --alignIntronMin $minIntronLength --alignIntronMax $maxIntronLength --genomeDir $index --runThreadN $threads --readFilesIn read_1.1,read_2.1 
