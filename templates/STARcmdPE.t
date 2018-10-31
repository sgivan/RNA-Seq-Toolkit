#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$threads
#SBATCH --mem=$mem
#SBATCH -time=0

STAR --genomeLoad NoSharedMemory --quantMode GeneCounts --outFileNamePrefix PE_ --runMode alignReads \
     --alignIntronMin $minIntronLength --alignIntronMax $maxIntronLength --genomeDir $index --runThreadN $threads \
     --alignMatesGapMax $alignGapMax --readFilesIn read_1 read_2
