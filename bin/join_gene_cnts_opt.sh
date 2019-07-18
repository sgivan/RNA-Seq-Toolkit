#!/bin/bash

TEMP=`getopt -o hc:e:E: -- "$@"`

help=0

cntrl=3
exp1=4
exp2=6

eval set -- "$TEMP"
while true ; do
    case "$1" in
        -h) help=1 ; shift ;;
        -c) cntrl=$2 ; shift 2 ;;
        -e) exp1=$2 ; shift 2 ;;
        -E) exp2=$2 ; shift 2 ;;
        --) shift ; break ;;
        *) break ;;
    esac
done

echo "last control is Sample_${cntrl}"
echo "first exp is Sample_${exp1}"
echo "last exp is Sample_${exp2}"

#cp Sample_1/gene_cnts.txt joined.txt
#for file in `ls -v Sample_{2..${cntrl}}/gene_cnts.txt`; do echo $file; join --header joined.txt $file > joined2.txt; mv joined2.txt joined.txt; done
#for file in `ls -v Sample_{${exp1}..${exp2}}/gene_cnts.txt`; do echo $file; join --header joined.txt $file > joined2.txt; mv joined2.txt joined.txt; done
#mv joined.txt C_v_E.txt

echo "controls"
echo "Sample_1/gene_cnts.txt"

cp Sample_1/gene_cnts.txt joined.txt

for (( dir=2; dir<=$cntrl; dir++ ));
do
    file="Sample_${dir}/gene_cnts.txt"
    ls $file
    join --header joined.txt $file > joined2.txt
    mv joined2.txt joined.txt;
done

echo "experimentals"

for (( dir=$exp1; dir<=$exp2; dir++ ));
do
    file="Sample_${dir}/gene_cnts.txt"
    ls $file
    join --header joined.txt $file > joined2.txt
    mv joined2.txt joined.txt;
done

mv joined.txt C_v_E.txt

echo ""
echo "final file is C_v_E.txt"

#cp Sample_1/gene_cnts.txt joined.txt
#for file in $(ls -v Sample_{2..$cntrl}/gene_cnts.txt); do echo $file; join --header joined.txt $file > joined2.txt; mv joined2.txt joined.txt; done
#for file in `ls -v Sample_{4..6}/gene_cnts.txt`; do echo $file; join --header joined.txt $file > joined2.txt; mv joined2.txt joined.txt; done
#mv joined.txt C_v_E.txt

#for dir
#do
#    echo $dir
#done


