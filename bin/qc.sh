#!/bin/bash

TEMP=`getopt -o qhcu -- "$@"`

quiet=0
header=0
comma=0
cumulative=0

eval set -- "$TEMP"
while true ; do
    case "$1" in
        -q) quiet=1 ; shift ;;
        -h) header=1 ; shift ;;
        -c) comma=1 ; shift ;;
        -u) cumulative=1 ; shift ;;
        --) shift ; break ;;
        *) break ;;
    esac
done

if [[ $header == 1 ]]
then 
    if [[ $comma == 1 ]]
    then
        echo "Sample,Raw,Trimmed,Filtered,Matched,Retained,% Retained"
    else
        echo "|Sample|Raw|Trimmed|Filtered|Matched|Retained|% Retained|"
    fi
fi

for dir
do
    cd $dir
    cd preprocess
    raw=`grep Input set1_qt.log | cut -f 2 -d " "`
    qt=`grep discarded set1_qt.log | awk '{ print $2; }'`
    qf=`grep discarded set1_qt_qf.log | awk '{ print $2; }'`
    qct=`grep 'failed to align' set1_qt_qf_bwt.log |  awk '{ print $7; }'`
    qcfi=`grep -m 1 reads set1_qt_qf_bwt.log | cut -f 4 -d " "`
    qm=`expr $qcfi - $qct`
    rtn=`echo "scale=4; ($qct/$raw)*100" | bc | sed 's/00\$//'`

    if [[ $cumulative -eq 1 ]]
    then
        qf=`expr $raw - $qt - $qf`
        qt=`expr $raw - $qt`
    fi

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
    if [[ $comma == 1 ]]
    then
        echo "$dir,$raw,$qt,$qf,$qm,$qct,$rtn"
    else
        echo "|$dir|$raw|$qt|$qf|$qm|$qct|$rtn|"
    fi
#    echo
done


