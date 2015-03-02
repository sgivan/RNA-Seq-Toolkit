#!/bin/bash

TEMP=`getopt -o qhc -- "$@"`

quiet=0
header=0
comma=0

eval set -- "$TEMP"
while true ; do
    case "$1" in
        -q) quiet=1 ; shift ;;
        -h) header=1 ; shift ;;
        -c) comma=1 ; shift ;;
        --) shift ; break ;;
        *) break ;;
    esac
done

if [[ $header == 1 ]]
then 
    if [[ $comma == 1 ]]
    then
        echo "Sample,QC Reads,Mapped Reads,% Mapped"
    else
        echo "|Sample|QC Reads|Mapped Reads|% Mapped|"
    fi
fi

for dir
do
#    echo $dir
    cd $dir
    cd preprocess
    qct=`grep 'failed to align' set1_qt_qf_bwt.log |  awk '{ print $7; }'`
    cd ../merged
    mapped=`cat merged.bam.cnt`
    rtn=`echo "scale=4; ($mapped/$qct)*100" | bc | sed 's/00\$//'`
    if [[ $quiet -eq 0 ]]
    then
        echo $dir
        echo "qc reads: $qct"
        echo "mapped reads: $mapped"
    fi
    cd ../../
    if [[ $comma == 1 ]]
    then
        echo "$dir,$qct,$mapped,$rtn"
    else
        echo "|$dir|$qct|$mapped|$rtn|"
    fi
done


