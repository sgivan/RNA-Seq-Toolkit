#!/bin/bash
#
# run this script from directory containing flowcell directories; ie, the one containing all the FCXXX directories
# then, invoke like FGMG_process_data.sh FC???/s_? > & FGMG_process_agg_data.log &
export PATH=~/projects/RNAseq/bin:$PATH
#
# initialize variables using default values
#run_type='full'
run_type='partial'
rd=`pwd`

function help_messg () {
    echo "invoke script with one of the following options:"
    echo "--full [will run full analysis, including short read preprocessing]"
    echo "--partial [will skip preprocessing steps]"
    echo "--aggregate [will use an aggregate junctions file and skips preprocessing]"
    echo "-h [print this help message]"
    echo ""
}

function mk_agg_txpts () {
    cd $rd
    echo "generating aggregate transcripts file"
    mkdir -p merged_aggregates
    cd merged_aggregates
    samtools merge all_merged.bam ../*/merged/merged.bam
    echo "cufflinks -p $threads -N --library-type $library_type -I 25000 -L allmerge -r ../$refseq all_merged.bam"
    cufflinks -p $threads -N --library-type $library_type -I 25000 -L allmerge -r ../$refseq all_merged.bam
    cd ..
    ln -s merged_aggregates/transcripts.gtf ./
#    ln -s merged_aggregates/aggregate_junctions.txt ./
}

#function mk_agg_jncts {
#
#}

mate_inner_distance=165
min_intron_length=50
##max_intron_length_I=10000
max_intron_length=25000
aggregate_junctions=0
aggregate_transcripts=0
refseq='refseq'
threads=8
library_type='fr-unstranded'
use_aggregates=0
seonly=0
adapter='NULL'
indexpath="$rd/index/"

# command line option parsing adpated from /usr/share/doc/util-linux-2.13/getopt-parse.bash
#
TEMP=`getopt -o pafhr:i:I:jts:H:l:ueA:P: --long full,aggregate,partial,mate_inner_distance:,min_intron_length:,max_intron_length:,agg_junctions,agg_transcripts,refseq:,threads:,library_type:,use_aggregates,seonly,adapter:,indexpath: -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!

eval set -- "$TEMP"
while true ; do
    case "$1" in
        -f|--full) run_type='full' ; shift ;;
        -p|--partial) run_type='partial' ; shift ;;
        -a|--aggregate) run_type='aggregate' ; shift ;;
        -r|--mate_inner_distance) mate_inner_distance=$2 ; shift 2 ;;
        -i|--min_intron_length) min_intron_length=$2 ; shift 2 ;;
        -I|--max_intron_length) max_intron_length=$2 ; shift 2 ;;
        -j|--agg_junctions) aggregate_junctions=1 ; shift ;;
        -t|--agg_transcripts) aggregate_transcripts=1 ; shift ;;
        -s|--refseq) refseq=$2 ; shift 2 ;;
        -H|--threads) threads=$2 ; shift 2 ;;
        -l|--library_type) library_type=$2 ; shift 2 ;;
        -u|--use_aggregates) use_aggregates=1 ; shift ;;
        -e|--seonly) seonly=1 ; shift ;;
        -A|--adapter) adapter=$2 ; shift 2 ;;
        -P|--indexpath) indexpath=$2 ; shift 2 ;;
        -h) help_messg ; exit ;;
        --) shift ; break ;;
        *) break ;;
    esac
done
#for arg do echo '--> '"\`$arg'" ; done

echo "run type is '$run_type'"
flags="--mate_inner_distance $mate_inner_distance --min_intron_length $min_intron_length --max_intron_length $max_intron_length --procs $threads --librarytype $library_type --indexpath $indexpath/"
echo "flags: $flags"
#exit

if [[ $seonly -eq 1 ]]
then
    flags="$flags --seonly"
fi

if [[ $adapter != "NULL" ]]
then
    flags="$flags --adapter_seq $adapter"    
fi

#echo "flags: '$flags'"
#exit

for dir 
do
	echo `date`
	#echo $rd/$dir
	#cd $rd/$dir
    echo $dir
    cd $dir

    case "$run_type" in

        partial)
	
            # for partial runs (when you don't need to run preprocessing steps
            echo "RNAseq.sh $flags --partial" ;
            RNAseq.sh $flags --partial ;;

        full)

            # for full runs
            echo "RNAseq.sh $flags --full" ;
            RNAseq.sh $flags --full ;;

        aggregate)

            # for aggregate junction runs
            # must have already run script with --agg_transcripts flag
            mkdir non-aggregate
            echo "moving old output files to 'non-aggregate'"
            mv pe_tophat* singles_tophat* non-aggregate/
            echo "creating symbolic link to transcript.gtf"
            ln -sf ../transcripts.gtf ./
            more_flags=""
            echo "running tophat"
            echo "RNAseq.sh --aggregate $flags $more_flags" ;
            RNAseq.sh --aggregate $flags $more_flags ;;

    esac
    cd $rd
done

if [ $aggregate_junctions = 1 ]
    then
        echo "generating aggregate junctions file"
        mkdir -p merged_aggregates
        cd merged_aggregates
        cat ../*/*/junctions.bed | awk '{ if ($1 != "track") {split($11,len,","); split($12,blstrt,","); printf "%s\t%i\t%i\t%s\n", $1, $2 + len[1] - 1, $2 + blstrt[2], $6; }}' | sort -k 1,1 -gk 2,2 | uniq > aggregate_junctions.txt
        cd ..
fi

if [ $aggregate_transcripts = 1 ]
    then
        mk_agg_txpts
fi

if [ $use_aggregates = 1 ]
    then
        echo "re-running pipeline using aggregate files"
fi

