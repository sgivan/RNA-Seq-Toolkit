#!/bin/bash

# Don't use this script
# Use qc.sh for a QC summary.
# use map_efficiency.sh for a mapping summary

TEMP=`getopt -o qHhcu -- "$@"`

quiet=0
header=0
comma=0
cumulative=0
hlp=0

eval set -- "$TEMP"
while true ; do
    case "$1" in
        -q) quiet=1 ; shift ;;
        -H) header=1 ; shift ;;
        -c) comma=1 ; shift ;;
        -u) cumulative=1 ; shift ;;
        -h) hlp=1 ; shift ;;
        --) shift ; break ;;
        *) break ;;
    esac
done

if [[ $hlp == 1 ]]
then
    echo "
    -q don't print verbose information
    -H print column headers
    -c print comma-separated data (default = | separated)
    -u print cumulative data
    -h print this help menu
    "
    exit
fi

if [[ $header == 1 ]]
then 
    if [[ $comma == 1 ]]
    then
        echo "Sample,Raw,Trimmed,Filtered,Matched,Retained,% Retained,PE Input,PE Mapped,%PE Mapped,SE Input,SE Mapped,%SE Mapped,All Mapped, % All Mapped"
    else
        echo "||= Sample =||= Raw^1^ =||= Trimmed^2^ =||= Filtered^3^ =||= Matched^4^ =||= Retained^5^ =||= % Retained =||= PE Input^6^ =||= PE Mapped^7^ =||= % PE Mapped  =||= SE Input^8^ =||= SE Mapped^9^ =||= % SE Mapped =||= All Mapped =||= % All Mapped^10^ =||"
    fi
fi

for dir
do
    cd $dir
    cd preprocess
    raw=$(grep Input set1_qt.log | cut -f 2 -d " ")
    qt=$(grep discarded set1_qt.log | awk '{ print $2; }')
    qf=$(grep discarded set1_qt_qf.log | awk '{ print $2; }')
    qct=$(grep 'failed to align' set1_qt_qf_bwt.log |  awk '{ print $7; }')
    qcfi=$(grep -m 1 reads set1_qt_qf_bwt.log | cut -f 4 -d " ")
    qm=$(expr $qcfi - $qct)
    rtn=$(echo "scale=4; ($qct/$raw)*100" | bc | sed 's/00\b//')

#    cd ../merged
#    mapped=$(cat merged.bam.cnt)
#    mrtn=$(echo "scale=4; ($mapped/$qct)*100" | bc | sed 's/00\b//')
#    cd ../../
#
    # strings to grep
    # Number of input reads
    # Uniquely mapped reads number
    # Number of reads mapped to multiple loci
    #
    # file1: PE_Log.final.out
    # file2: SE_Log.final.out

    # Paired-End (PE)
    # PE_input_reads is actually *pairs* of reaads
    PE_input_reads=$(grep 'Number of input reads' PE_Log.final.out | cut -f 2)
    #echo "PE_input_reads: '"$PE_input_reads"'"
    PE_unique_mapped=$(grep 'Uniquely mapped reads number' PE_Log.final.out | cut -f 2)
    #echo "PE_unique_reads: '"$PE_unique_mapped"'"
    PE_multi_mapped=$(grep 'Number of reads mapped to multiple loci' PE_Log.final.out | cut -f 2)
    #echo "PE_multi_mapped: '"$PE_multi_mapped"'"
    PE_mapped=$(expr $PE_unique_mapped + $PE_multi_mapped)

    # Single-End (SE)
    SE_input_reads=$(grep 'Number of input reads' SE_Log.final.out | cut -f 2)
    #echo "SE_input_reads: '"$SE_input_reads"'"
    SE_unique_mapped=$(grep 'Uniquely mapped reads number' SE_Log.final.out | cut -f 2)
    #echo "SE_unique_reads: '"$SE_unique_mapped"'"
    SE_multi_mapped=$(grep 'Number of reads mapped to multiple loci' SE_Log.final.out | cut -f 2)
    #echo "SE_multi_mapped: '"$SE_multi_mapped"'"
    SE_mapped=$(expr $SE_unique_mapped + $SE_multi_mapped)

    if [[ $cumulative -eq 1 ]]
    then
        qf=$(expr $raw - $qt - $qf)
        qt=$(expr $raw - $qt)
    fi

    if [[ $quiet -eq 0 ]]
    then
        echo $dir
        echo "trimmed: $qt"
        echo "filtered: $qf"
        echo "matched: $qm"
        echo "qc reads: $qct"
        echo "% retained: $rtn"
        echo "mapped: $mapped"
        echo "% mapped: $mrtn"
    fi

    if [[ $comma == 1 ]]
    then
        echo "$dir,$raw,$qt,$qf,$qm,$qct,$rtn,$mapped,$mrtn"
    else
        echo "|$dir|$raw|$qt|$qf|$qm|$qct|$rtn|$mapped|$mrtn|"
    fi

done

if [[ $header == 1 ]]
then 
    echo ""
    echo "  - Number of raw reads at beginning of QC regimen"
    echo "  - Number of reads that had low quality base calls trimmed"
    echo "  - Number of reads filtered due to overall low quality"
    echo "  - Number of reads that match PhiX, rRNA genes or relevant repeat sequences"
    echo "  - Number of reads remaining after all quality control steps"
    echo "  - Number of QC reads mapped to reference genome"
fi
