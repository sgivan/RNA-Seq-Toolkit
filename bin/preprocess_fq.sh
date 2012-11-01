#!/bin/bash
# copyright Scott Givan, The University of Missouri, July 6, 2012
#
#    This file is part of the RNA-seq Toolkit, or RST.
#
#    RST is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RST is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RST.  If not, see <http://www.gnu.org/licenses/>.
#
#export PATH="~/projects/RNAseq/bin/:$PATH"
wd=`pwd`
osname=`uname -s`
export PATH=".:$wd:$wd/bin:$HOME/bin:$PATH"
#echo $PATH
#exit
#
# initialize variables with default values
# which can be superseded with clo's
#
min_qual=13
min_length=32
percent_high_quality=90
bowtie_threads=8
BOWTIE_INDEXES='index'
BOWTIE2_INDEXES='index'
bowtie_cmd='bowtie2'
filter='filter'
leave_temp=0
qualscores='NULL'
seonly=0
#
function help_messg {
    echo "invoke script like this:"
    echo "preprocess_fq.sh min_qual min_length %high_quality #bowtie_threads"
    echo "or, to accept the default values:"
    echo "preprocess_fq.sh"
    echo "example with the default values:"
    echo "preprocess.sh --min_qual 13 --min_length 32 --percent_high_quality 90 --bowtie_threads 8"
    echo ""
}

# command line option parsing adpated from /usr/share/doc/util-linux-2.13/getopt-parse.bash
#
case "$osname" in

    Linux)
        TEMP=`getopt -o q:l:p:t:hi:f:cQ:L:H:se --long min_qual:,min_length:,percent_high_quality:,bowtie_threads:,indexpath:,filter:,leave_temp,min_qual:,min_length:,percent_high_quality:,solexa,seonly -- "$@"`
        ;;

    Darwin)
        TEMP=`getopt q:l:p:t:hi:f:cQ:L:H:se $*`
        ;;

    *)
        TEMP=`getopt -o q:l:p:t:hi:f:cQ:L:H:se --long min_qual:,min_length:,percent_high_quality:,bowtie_threads:,indexpath:,filter:,leave_temp,min_qual:,min_length:,percent_high_quality:,solexa,seonly -- "$@"`
        ;;
esac

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

while true ; do
    case "$1" in
        -q|--min_qual) min_qual=$2 ; shift 2 ;;
        -l|--min_length) min_length=$2 ; shift 2 ;;
        -p|--percent_high_quality) percent_high_quality=$2 ; shift 2 ;; 
        -t|--bowtie_threads) bowtie_threads=$2 ; shift 2 ;;
        -i|--indexpath) BOWTIE_INDEXES=$2 ; shift 2 ;;
        -f|--filter) filter=$2 ; shift 2 ;;
        -Q|--min_qual) min_qual=$2 ; shift 2 ;;
        -L|--min_length) min_length=$2 ; shift 2 ;;
        -H|--percent_high_quality) percent_high_quality=$2 ; shift 2 ;;
        -s|--solexa) qualscores=1 ; shift ;;
        -e|--seonly) seonly=1 ; shift ;;
        -h|--help) help_messg ; exit ;;
        -c|--leave_temp) leave_temp=1; shift ;;
#        --) shift ; break ;;
        *) break ;;
    esac
done

export BOWTIE_INDEXES
export BOWTIE2_INDEXES
echo 'BOWTIE2_INDEXES: '$BOWTIE2_INDEXES

echo "min_qual min_length percent_high_quality bowtie_threads indexpath = $min_qual $min_length $percent_high_quality $bowtie_threads $BOWTIE2_INDEXES"
#
echo "creating preprocess directory in $wd"
mkdir -p preprocess
#mv set1.fq set1_input.fq
#mv set2.fq set2_input.fq
cd preprocess
ln -sf ../set1.fq ./set1.fq
ln -sf ../../index
#exit
if [[ $seonly -ne 1 ]]
then
    ln -sf ../set2.fq ./set2.fq
fi

echo "quality trimming and filtering"
trimmer_flags="-t $min_qual -l $min_length -v"
filter_flags="-p $percent_high_quality -q $min_qual"
bowtie_flags="-q --threads $bowtie_threads"
#echo trimmer_flags="-t $min_qual -l $min_length -v"
#echo filter_flags="-p $percent_high_quality -q $min_qual"
#echo bowtie_flags="-q --threads $bowtie_threads" 
#exit
if [[ $qualscores != 'NULL' ]]
then
#    trimmer_flags="$trimmer_flags -Q 33"
#    filter_flags="$filter_flags -Q 33"
    #bowtie_flags="$bowtie_flags --solexa-quals"
    bowtie_flags="$bowtie_flags --phred64-quals"
else
#    bowtie_flags="$bowtie_flags --solexa1.3-quals"
    trimmer_flags="$trimmer_flags -Q 33"
    filter_flags="$filter_flags -Q 33"
fi
echo trimmer_flags="-t $min_qual -l $min_length -v"
echo filter_flags="-p $percent_high_quality -q $min_qual"
echo bowtie_flags="-q --threads $bowtie_threads" 
#exit
#echo "fastq_quality_trimmer -i set1.fq -t $min_qual -l $min_length -v 2> set1_qt.log | fastq_quality_filter -p $percent_high_quality -q $min_qual -o set1_qt_qf.fq -v | tee set1_qt_qf.log" 
echo "fastq_quality_trimmer -i set1.fq $trimmer_flags 2> set1_qt.log | fastq_quality_filter $filter_flags -o set1_qt_qf.fq -v | tee set1_qt_qf.log" 
#fastq_quality_trimmer -i set1.fq $trimmer_flags 2> set1_qt.log | fastq_quality_filter $filter_flags -o set1_qt_qf.fq -v | tee set1_qt_qf.log 
eval fastq_quality_trimmer -i set1.fq $trimmer_flags 2> set1_qt.log | fastq_quality_filter $filter_flags -o set1_qt_qf.fq -v | tee set1_qt_qf.log 
if [[ $seonly -ne 1 ]]
then
    #echo "fastq_quality_trimmer -i set2.fq -t $min_qual -l $min_length -v 2> set2_qt.log | fastq_quality_filter -p $percent_high_quality -q $min_qual -o set2_qt_qf.fq -v | tee set2_qt_qf.log"
    echo "fastq_quality_trimmer -i set2.fq $trimmer_flags 2> set2_qt.log | fastq_quality_filter $filter_flags -o set2_qt_qf.fq -v | tee set2_qt_qf.log"
    #fastq_quality_trimmer -i set2.fq $trimmer_flags 2> set2_qt.log | fastq_quality_filter $filter_flags -o set2_qt_qf.fq -v | tee set2_qt_qf.log
    eval fastq_quality_trimmer -i set2.fq $trimmer_flags 2> set2_qt.log | fastq_quality_filter $filter_flags -o set2_qt_qf.fq -v | tee set2_qt_qf.log
fi
#exit
echo "sequence similarity filtering using $bowtie_cmd"
echo "reads that fail to align are retained, reads that align are filtered out of data"
echo "first data file ..."
#echo "bowtie -q --solexa1.3-quals --un set1_qt_qf_sf.fq --threads $bowtie_threads $filter set1_qt_qf.fq > set1_qt_qf_filter_matched.sam 2> set1_qt_qf_bwt.log"
echo "$bowtie_cmd $bowtie_flags --un set1_qt_qf_sf.fq $filter set1_qt_qf.fq > set1_qt_qf_filter_matched.sam 2> set1_qt_qf_bwt.log"
$bowtie_cmd $bowtie_flags --un set1_qt_qf_sf.fq $filter set1_qt_qf.fq > set1_qt_qf_filter_matched.sam 2> set1_qt_qf_bwt.log
cat set1_qt_qf_bwt.log

if [[ -e set2_qt_qf.fq ]]
then
    echo "second data file"
    echo "$bowtie_cmd $bowtie_flags --un set2_qt_qf_sf.fq $filter set2_qt_qf.fq > set2_qt_qf_filter_matched.sam 2> set2_qt_qf_bwt.log"
    $bowtie_cmd $bowtie_flags --un set2_qt_qf_sf.fq $filter set2_qt_qf.fq > set2_qt_qf_filter_matched.sam 2> set2_qt_qf_bwt.log
    cat set2_qt_qf_bwt.log
fi

if [[ leave_temp -eq 1 ]]
then
    echo "removing temporary files"
    rm *.sam
    #rm set?_qt_qf_sf.fq
    rm set?_qt_qf.fq
fi

echo "creating symlinks to final files"
#
# why are the following lines commented?
# uncommenting because adapter_trim.pl will use set1.fq and set2.fq as input files
# however, now adapter_trim.pl runs before preprocess_fq.sh, so re-comment
#mv set1.fq set1_input.fq
#mv set2.fq set2_input.fq
ln -sf set1_qt_qf_sf.fq set1.fq
if [[ $seonly -ne 1 ]]
then
    ln -sf set2_qt_qf_sf.fq set2.fq
fi
cd ..
#ln -sf preprocess/set1_qt_qf_sf.fq set1.fq
#ln -sf preprocess/set2_qt_qf_sf.fq set2.fq
ln -sf preprocess/set1.fq ./
if [[ $seonly -ne 1 ]]
then
    ln -sf preprocess/set2.fq ./
fi

