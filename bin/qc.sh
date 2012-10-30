#!/bin/bash

for dir
do
    echo $dir
    cd $dir
#    total=`fqseqs set1_00.fq`
#    echo "raw seqs $total"
    cd preprocess
#    cat *.log
    raw=`grep Input set1_qt.log | cut -f 2 -d " "`
    qt=`grep discarded set1_qt.log | awk '{ print $2; }'`
    qf=`grep discarded set1_qt_qf.log | awk '{ print $2; }'`
    qm=`grep Reported set1_qt_qf_bwt.log | cut -f 2 -d " "`
    echo "trimmed: $qt"
    echo "filtred: $qf"
    echo "matched: $qm"
    cd ../../
    echo "|$dir|$raw|$qt|$qf|$qm|"
    echo
done


