#!/bin/bash

TEMP=`getopt -o qhcu -- "$@"`

quiet=0
header=0
comma=0
cumulative=0
paired=0

eval set -- "$TEMP"
while true ; do
    case "$1" in
        -q) quiet=1 ; shift ;;
        -h) header=1 ; shift ;;
        -c) comma=1 ; shift ;;
        -u) cumulative=1 ; shift ;;
        -p) paired=1 ; shift ;;
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
        #echo "^Sample^Raw<sup>1</sup>^Trimmed<sup>2</sup>^Filtered<sup>3</sup>^Matched<sup>4</sup>^Retained<sup>5</sup>^% Retained^"
        echo "||= Sample =||= Raw^1^ =||= Trimmed^2^ =||= Filtered^3^ =||= Matched^4^ =||= Retained^5^ =||= % Retained =||"
    fi
fi

for dir
do
    cd $dir
    cd preprocess

    if [[ -e 'set2_qt.log' ]]
    then
        let paired=1
#        echo "paired = '"$paired"'"
    fi

    raw=$(grep Input set1_qt.log | cut -f 2 -d " ")
    qt=$(grep discarded set1_qt.log | awk '{ print $2; }')
    qf=$(grep discarded set1_qt_qf.log | awk '{ print $2; }')
    qct=$(grep 'failed to align' set1_qt_qf_bwt.log |  awk '{ print $7; }')
    qcfi=$(grep -m 1 reads set1_qt_qf_bwt.log | cut -f 4 -d " ")
    qm=$(expr $qcfi - $qct)
    rtn=$(echo "scale=4; ($qct/$raw)*100" | bc | sed 's/00\b//')

    if [[ $cumulative -eq 1 ]]
    then
        qf=$(expr $raw - $qt - $qf)
        qt=$(expr $raw - $qt)
    fi

    if [[ $quiet -eq 0 ]]
    then
        echo $dir
        echo "Read set 1"
        echo "trimmed: $qt"
        echo "filtered: $qf"
        echo "matched: $qm"
        echo "qc reads: $qct"
        echo "% retained: $rtn"
    fi
#    cd ../../

    if [[ $comma == 1 ]]
    then
        echo "$dir Read1,$raw,$qt,$qf,$qm,$qct,$rtn"
    else
        echo "||$dir Read1||$raw||$qt||$qf||$qm||$qct||$rtn||"
    fi
#    echo

#    if (($paired > 0));
    if [[ $paired -gt 0 ]];
    then

        raw2=$(grep Input set2_qt.log | cut -f 2 -d " ")
        qt2=$(grep discarded set2_qt.log | awk '{ print $2; }')
        qf2=$(grep discarded set2_qt_qf.log | awk '{ print $2; }')
        qct2=$(grep 'failed to align' set2_qt_qf_bwt.log |  awk '{ print $7; }')
        qcfi2=$(grep -m 1 reads set2_qt_qf_bwt.log | cut -f 4 -d " ")
        qm2=$(expr $qcfi2 - $qct2)
        rtn2=$(echo "scale=4; ($qct2/$raw2)*100" | bc | sed 's/00\b//')

        if [[ $cumulative -eq 1 ]]
        then
            qf2=$(expr $raw2 - $qt2 - $qf2)
            qt2=$(expr $raw2 - $qt2)
        fi

        if [[ $quiet -eq 0 ]]
        then
            echo $dir
            echo "Read set 2"
            echo "trimmed: $qt2"
            echo "filtered: $qf2"
            echo "matched: $qm2"
            echo "qc reads: $qct2"
            echo "% retained: $rtn2"
        fi

        if [[ $comma == 1 ]]
        then
            echo "$dir Read2,$raw2,$qt2,$qf2,$qm2,$qct2,$rtn2"
        else
            echo "||$dir Read2||$raw2||$qt2||$qf2||$qm2||$qct2||$rtn2||"
        fi
#        cd ../../

    fi
    cd ../../
done


if [[ $header == 1 ]]
then 
    echo ""
    echo "1. Number of raw reads at beginning of QC regimen"
    echo "2. Number of reads that had low quality base calls trimmed"
    echo "3. Number of reads filtered due to overall low quality"
    echo "4. Number of reads that match PhiX, rRNA genes or relevant repeat sequences"
    echo "5. Number of reads remaining after all quality control steps"
fi
