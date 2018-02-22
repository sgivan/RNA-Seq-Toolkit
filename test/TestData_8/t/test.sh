#!/bin/bash
# run_test.sh
# Tests non-paired sequence reads.
# Runs initial mapping to reference genome.
# Creates a merged transcript gtf file.
# Then, re-runs HISAT2 pipepline using the merged GTF for all samples.

DEBUG=0

# Sanity check (Are we where we want to be with the files we expect?)
debug_check_pwd_and_files() {
    if [ $DEBUG -eq 1 ]; then 
        echo "------------- Files in `pwd` and subfolders -------------"
        ls -A *
        echo "---------------------------------------------------------"
    fi
}

debug_check_pwd_and_files

#EXPLANATION OF FLAGS USED (listed alphabetically)
# -a use "transcripts.gtf" for gene models and skip preprocessing
# -e Single-end reads
# -f --full (will run full analysis, including short read preprocessing)
# -H number of threads to use
# -t generate gtf file of empirical transcripts
# -Y Phred quality values encoded as Phred + 64 
RNAseq_process_data.sh -e -H 4 -t -f -Y s_?

debug_check_pwd_and_files

RNAseq_process_data.sh -e -H 4 -a -Y  s_?

debug_check_pwd_and_files

# run ballgown
../../bin/ballgown_setup.pl --cont "s_[56]" --exp "s_[78]"

debug_check_pwd_and_files

ballgown.Rscript 

debug_check_pwd_and_files

# Test result files
t/test.t $HISAT_VERSION_BEING_TESTED
