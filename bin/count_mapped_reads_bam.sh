#!/bin/bash
file=$1
#samtools view merged.bam | cut -f 1 | psort --parallel=4 --temporary-directory=/scratch | uniq | wc -l > merged.cnt
#samtools view $file | cut -f 1 | psort --parallel=4 --temporary-directory=/scratch | uniq | wc -l > ${file}.cnt
# use following line on BioCluster
#samtools view $file | grep -v '^@' | awk '{ if ($2 != 4) { print $0; } }' | cut -f 1 | psort --parallel=4 --temporary-directory=/scratch | uniq | wc -l > ${file}.cnt
# for most machines, don't use psort and /scratch, b/c they may not be available
# psort is just a link to a newer version of gnu sort that is threaded
# https://www.gnu.org/software/coreutils/manual/html_node/sort-invocation.html
#
samtools view $file | grep -v '^@' | awk '{ if ($2 != 4) { print $0; } }' | cut -f 1 | sort --temporary-directory /scratch | uniq | wc -l > ${file}.cnt
