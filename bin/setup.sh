#!/bin/bash

#for dir in $(ls -1d ../data/current/Sample_*)
for dir in $(ls -1d ${1}/Sample_*)
do
    echo $dir
    newdir=$( echo $dir | sed -E 's/.+\///' )
    echo "newdir: "$newdir
    mkdir -p ${newdir}
    cd ${newdir}
    file=$( ls -1 ../$dir/*.fq | head -n 1 )
    ln -sf $file ./set1_00.fq
    ln -sf set1_00.fq set1.fq
    ln -s set1.fq read_1
    file=$( ls -1 ../$dir/*.fq | tail -n 1 )
    ln -sf $file ./set2_00.fq
    ln -sf set2_00.fq set2.fq
    ln -s set2.fq read_2
    cd ..
done

#ln -s index hisat_index
