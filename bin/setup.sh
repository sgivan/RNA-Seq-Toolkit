#!/bin/bash

for dir in $(ls -1d ../data/current/Sample_*)
do
    echo $dir
    newdir=$( echo $dir | sed -E 's/.+\///' )
    echo "newdir: "$newdir
    mkdir -p ${newdir}
    cd ${newdir}
    file=$( ls -1 ../$dir/*.fastq | head -n 1 )
    ln -sf $file ./set1_00.fq
    ln -sf set1_00.fq set1.fq
    file=$( ls -1 ../$dir/*.fastq | tail -n 1 )
    ln -sf $file ./set2_00.fq
    ln -sf set2_00.fq set2.fq
    cd ..
done

ln -s index hisat_index
