#!/bin/bash
# run_test.sh
# Tests non-paired sequence reads.
# Runs initial mapping to reference genome.
# Creates a merged transcript gtf file.
# Then, re-runs HISAT2 pipepline using the merged GTF for all samples.

HISAT_VERSION_BEING_TESTED=2.0.4

module load Python-shared
module load R-3.3.0-sharedlib
module load bowtie2-2.3.2
module load stringtie-1.3.0
module load HISAT2-$HISAT_VERSION_BEING_TESTED

echo "running using hisat2 version '`hisat2 --version | head -n1`' (expected $HISAT_VERSION_BEING_TESTED)"

# Sanity check (Are we where we want to be and have the files been copied?):
sanity_check() {
    echo ""
    echo "---------------------------------------------------------"
    echo "------------- Files in `pwd` and subfolders -------------"
    ls *
    echo "---------------------------------------------------------"
    echo "---------------------------------------------------------"
    echo ""
}

sanity_check

#EXPLANATION OF FLAGS USED (listed alphabetically)
# -a use "transcripts.gtf" for gene models and skip preprocessing
# -e Single-end reads
# -f --full (will run full analysis, including short read preprocessing)
# -H number of threads to use
# -t generate gtf file of empirical transcripts
# -Y Phred quality values encoded as Phred + 64 
../../bin/RNAseq_process_data.sh -e -H 4 -t -f -Y s_?

sanity_check

../../bin/RNAseq_process_data.sh -e -H 4 -a -Y  s_?

sanity_check

# run ballgown
bin/ballgown_setup.pl --cont "s_[56]" --exp "s_[78]"

sanity_check

../../R/ballgown.Rscript 

sanity_check

# Test result files
t/test.pl 2.0.4


