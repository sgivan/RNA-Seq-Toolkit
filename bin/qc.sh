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
    raw=`grep Input set1_qt.log | cut -f 2 -d " "`
    qt=`grep discarded set1_qt.log | awk '{ print $2; }'`
    qf=`grep discarded set1_qt_qf.log | awk '{ print $2; }'`
    qct=`grep 'aligned 0 times' set1_qt_qf_bwt.log | awk '{ print $1; }'`
    qcfi=`grep reads set1_qt_qf_bwt.log | cut -f 1 -d " "`
    qm=`expr $qcfi - $qct`
    rtn=`echo "scale=4; ($qct/$raw)*100" | bc | sed 's/00\$//'`
    if [[ $quiet -eq 0 ]]
    then
        echo $dir
        echo "trimmed: $qt"
        echo "filtered: $qf"
        echo "matched: $qm"
        echo "qc reads: $qct"
        echo "% retained: $rtn"
    fi
    cd ../../
    echo "|$dir|$raw|$qt|$qf|$qm|$qct|$rtn|"
#    echo
done


