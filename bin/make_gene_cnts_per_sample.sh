#!/bin/bash

TEMP=`getopt -o p -- "$@"`

paired=0

eval set -- "$TEMP"
while true ; do
    case "$1" in
        -p) paired=1 ; shift ;;
        --) shift ; break ;;
        *) break ;;
    esac
done

for dir in Sample_*;
do
    echo $dir;
    cd $dir;
#    $(echo -e "EMBLID\t${dir}" > gene_cnts.txt)
    $(echo -e "GENEID\t${dir}" > gene_cnts.txt)

    $(grep -v "N_" merged_ReadsPerGene.out.tab | cut -f 1,4 >> gene_cnts.txt)

    cd ..

done
