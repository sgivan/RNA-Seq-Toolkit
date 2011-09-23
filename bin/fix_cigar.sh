#!/bin/bash
file=$1

if [ $# != 1 ]
then
        echo "what file?"
        exit
fi
grep -v N536870...M $file | samtools view -b -S - > accepted_hits.bam

