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
#    qm=`grep Reported set1_qt_qf_bwt.log | cut -f 2 -d " "`# doesn't work with bowtie2
    qct=`grep 'aligned 0 times' set1_qt_qf_bwt.log | awk '{ print $1; }'`
    qcfi=`grep reads set1_qt_qf_bwt.log | cut -f 1 -d " "`
    qm=`expr $qcfi - $qct`
    rtn=`echo "scale=4; ($qct/$raw)*100" | bc | sed 's/00\$//'`
    echo "trimmed: $qt"
    echo "filtered: $qf"
    echo "matched: $qm"
    echo "qc reads: $qct"
    echo "% retained: $rtn"
    cd ../../
    echo "|$dir|$raw|$qt|$qf|$qm|$qct|$rtn|"
    echo
done


