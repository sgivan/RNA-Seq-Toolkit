#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$threads
#SBATCH --mem=$mem

STAR --genomeLoad LoadAndRemove --quantMode GeneCounts --outFileNamePrefix PE_ --runMode alignReads \
    --genomeDir ../../index --runThreadN $threads --readFilesIn ../read_1 ../read_2
