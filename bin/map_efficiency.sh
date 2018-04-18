#!/bin/bash

TEMP=`getopt -o qhcp -- "$@"`

quiet=0
header=0
comma=0
paired=0

eval set -- "$TEMP"
while true ; do
    case "$1" in
        -q) quiet=1 ; shift ;;
        -h) header=1 ; shift ;;
        -c) comma=1 ; shift ;;
        -p) paired=1 ; shift ;;
        --) shift ; break ;;
        *) break ;;
    esac
done

if [[ $paired == 1 ]]
then

    if [[ $header == 1 ]]
    then 
        if [[ $comma == 1 ]]
        then
            echo "Sample,QC Reads,PE Reads,PE Mapped,%PE Mapped,SE Reads,SE Mapped,%SE Mapped,All Mapped, % All Mapped"
        else
            echo "||= Sample =||= QC Reads^1^ =||= PE QC^2^ =||= PE Mapped^3^ =||= % PE Mapped  =||= SE QC^4^ =||= SE Mapped^5^ =||= % SE Mapped =||= All Mapped^6^ =||= % All Mapped =||"
        fi
    fi

    for dir
    do
        cd $dir
        cd preprocess
        qct1=`grep 'failed to align' set1_qt_qf_bwt.log |  awk '{ print $7; }'`
        qct2=`grep 'failed to align' set2_qt_qf_bwt.log |  awk '{ print $7; }'`
        qctAll=$(expr $qct1 + $qct2)
        cd ..
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

        mapped=$(expr 2 \* $PE_mapped + $SE_mapped)


        ppe_mapped=`echo "scale=4; ($PE_mapped/$PE_input_reads)*100" | bc | sed 's/00\$//'`
        pse_mapped=`echo "scale=4; ($SE_mapped/$SE_input_reads)*100" | bc | sed 's/00\$//'`
        rtn=`echo "scale=4; ($mapped/$qctAll)*100" | bc | sed 's/00\$//'`
        if [[ $quiet -eq 0 ]]
        then
            echo $dir
            echo "qc reads: $qctAll"
            echo "mapped reads: $mapped"
        fi

        if [[ $comma == 1 ]]
        then

#            echo "Sample,QC Reads,PE Reads,PE Mapped,%PE Mapped,SE Reads,SE Mapped,%SE Mapped,All Mapped, % All Mapped"
            echo "$dir,$qctAll,$PE_input_reads,$PE_mapped,$ppe_mapped,$SE_input_reads,$SE_mapped,$pse_mapped,$mapped,$rtn"
        else
            echo "||$dir||$qctAll||$PE_input_reads||$PE_mapped||$ppe_mapped||$SE_input_reads||$SE_mapped||$pse_mapped||$mapped||$rtn||"
        fi
        cd ..
    done
else

    if [[ $header == 1 ]]
    then 
        if [[ $comma == 1 ]]
        then
            echo "Sample,QC Reads,Mapped Reads,% Mapped"
        else
            echo "||Sample||QC Reads||Mapped Reads||% Mapped||"
        fi
    fi

    for dir
    do
    #    echo $dir
        cd $dir
        cd preprocess
        qct=`grep 'failed to align' set1_qt_qf_bwt.log |  awk '{ print $7; }'`
        cd ..

        # Single-End (SE)
        SE_input_reads=$(grep 'Number of input reads' SE_Log.final.out | cut -f 2)
        #echo "SE_input_reads: '"$SE_input_reads"'"
        SE_unique_mapped=$(grep 'Uniquely mapped reads number' SE_Log.final.out | cut -f 2)
        #echo "SE_unique_reads: '"$SE_unique_mapped"'"
        SE_multi_mapped=$(grep 'Number of reads mapped to multiple loci' SE_Log.final.out | cut -f 2)
        #echo "SE_multi_mapped: '"$SE_multi_mapped"'"
        mapped=$(expr $SE_unique_mapped + $SE_multi_mapped)

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
            echo "||$dir||$qct||$mapped||$rtn||"
        fi
    done
fi

if [[ $header == 1 ]]
then 
    if [[ $paired == 1 ]]
    then
#       echo "||= Sample =||= QC Reads^1^ =||= PE Input^2^ =||= PE Mapped^3^ =||= % PE Mapped  =||= SE Input^4^ =||= SE Mapped^5^ =||= % SE Mapped =||= All Mapped^6^ =||= % All Mapped =||"
        echo "1. During the QC phase, some reads are removed from the experiment based on a variety of factors (low quality, below a minimum length, matches a repeat element, etc.). In some cases, only a single member of a paired set of reads is removed. The reads that are '''not''' removed are subsequently referred to as Single End (SE) reads. We collect and map SE reads independently of the PE reads. The number in this column represents the Single End (SE) + Paired End (PE) reads remaining after all quality control (QC) steps. Note that this number represents individual reads."
        echo "2. Number of  Paired End (PE) reads after QC. Note that this number represents '''pairs''' of reads. So, 2 X (# of PE reads) = # of individual PE reads."
        echo "3. Number of PE reads that can be mapped to the reference genome. Note that this number represents '''pairs''' of reads."
        echo "4. Number of Single End (SE) reads after QC. Note that this number represents individual reads."
        echo "5. Number of SE reads that can be mapped to the reference genome. Note that this number represents individual reads."
        echo "6. Number of SE + PE reads that can be mapped to the reference genome. This number represents individual reads."

    else

#       echo "||Sample||QC Reads^1^||Mapped Reads^2^||% Mapped||"
        echo ""
        echo "1. Number of reads remaining after all quality control steps"
        echo "2. Number of reads mapped to the reference genome."
    fi
fi
