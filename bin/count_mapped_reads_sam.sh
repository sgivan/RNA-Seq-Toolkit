#!/bin/bash
file=$1
# Only use following line on IRCF BioCluster
# It uses a newer version of Gnu sort that is threaded
# https://www.gnu.org/software/coreutils/manual/html_node/sort-invocation.html
# It also specifies a temporary directory that probably isn't available on most machines.
#grep -v '^@' $file | awk 'BEGIN { FS = "\t" } { if ($2 != 4) { print $0; } }' | cut -f 1 | psort --parallel=4 --temporary-directory=/scratch | uniq | wc -l > ${file}.cnt
grep -v '^@' $file | awk 'BEGIN { FS = "\t" } { if ($2 != 4) { print $0; } }' | cut -f 1 | sort | uniq | wc -l > ${file}.cnt
