#!/bin/bash

TEMP=`getopt -o q -- "$@"`

quiet=0

eval set -- "$TEMP"
while true ; do
    case "$1" in
        -q) quiet=1 ; shift ;;
        --) shift ; break ;;
        *) break ;;
    esac
done

for dir
do
#    echo $dir
    cd $dir
    cd preprocess
    qct=`grep 'aligned 0 times' set1_qt_qf_bwt.log | awk '{ print $1; }'`
    cd ../merged
#    mapped=`cat merged.cnt`
    mapped=`cat merged.bam.cnt`
    rtn=`echo "scale=4; ($mapped/$qct)*100" | bc | sed 's/00\$//'`
    if [[ $quiet -eq 0 ]]
    then
        echo $dir
        echo "qc reads: $qct"
        echo "mapped reads: $mapped"
    fi
    cd ../../
    echo "|$dir|$qct|$mapped|$rtn|"
done


