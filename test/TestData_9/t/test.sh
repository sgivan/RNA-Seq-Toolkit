#!/bin/bash
# run_test.sh
# Tests non-paired sequence reads.
# Starts with reference GTF file
# Creates a merged transcript gtf file.
# Then, re-runs HISAT2 pipepline using the merged GTF for all samples.
#TODO: Are all those notes above correct for this test?


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

echo 'running RNAseq_process_data.sh -e -H 4 -O -C s_?'
#EXPLANATION OF FLAGS USED (listed alphabetically)
# -a use "transcripts.gtf" for gene models and skip preprocessing
# -C leave temporary files on file system
# -e Single-end reads
# -H number of threads to use
# -O preprocess data
# -X Phred quality values encoded as Phred + 33
RNAseq_process_data.sh -e -H 4 -O -C s_?

debug_check_pwd_and_files

echo 'running RNAseq_process_data.sh -e -H 4 -a -X s_?'
RNAseq_process_data.sh -e -H 4 -a -X s_?

debug_check_pwd_and_files

# run ballgown
ballgown_setup.pl --cont "s_[13]" --exp "s_[24]"

debug_check_pwd_and_files

ballgown.Rscript

debug_check_pwd_and_files

# Test result files
t/test.t
