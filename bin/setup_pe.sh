#!/bin/bash
#
# 

if [[ $# -eq 0 ]]; then
    echo "

    To use this script, there needs to be
    a directory ../data/current with directories
    named like 'Sample_*'. Or, you can pass
    the name of the data directory as an
    argument. Also, there needs
    to be an 'index' directory (or symlink)
    in the current directory.

    "
    datadir='../data/current'

elif [[ $# -eq 1 ]]; then
    datadir=$1
fi


#for dir in $(ls -1d ../data/current/Sample_*)
for dir in $(ls -1d ${datadir}/Sample_*)
do
    echo $dir
    newdir=$( echo $dir | sed -E 's/.+\///' )
    echo "newdir: "$newdir
    mkdir -p ${newdir}
    cd ${newdir}
    file=$( ls -1 ../$dir/*.fastq | head -n 1 )
    ln -sf $file ./set1_00.fq
    ln -sf set1_00.fq set1.fq
#    cd ..
#    cd ${newdir}-2
    file=$( ls -1 ../$dir/*.fastq | tail -n 1 )
    ln -sf $file ./set2_00.fq
    ln -sf set2_00.fq set2.fq
    cd ..
done

ln -s index hisat_index
