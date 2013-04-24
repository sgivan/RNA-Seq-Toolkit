#!/bin/bash

for dir
do
    echo $dir
    cd $dir
    cd preprocess
    qct=`grep 'aligned 0 times' set1_qt_qf_bwt.log | awk '{ print $1; }'`
    cd ../merged
    mapped=`cat merged.cnt`
    rtn=`echo "scale=4; ($mapped/$qct)*100" | bc | sed 's/00\$//'`
    echo "qc reads: $qct"
    echo "mapped reads: $mapped"
    cd ../../
    echo "|$qct|$mapped|$rtn|"
done


