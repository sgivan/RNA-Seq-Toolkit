#!/bin/bash
#
# run this script from directory containing flowcell directories; ie, the one containing all the FCXXX directories
# then, invoke like FGMG_process_data.sh FC???/s_? 
wd=`pwd`
osname=`uname -s`
#export PATH=~givans/projects/RNAseq/bin:$PATH
#export PATH=/ircf/sgivan/projects/RNAseq/bin:$PATH
export PATH=".:$wd:$wd/bin:$HOME/bin:$PATH"
#

function help_messg () {
    case "$osname" in

        Linux) 
            echo ""
            echo "invoke script with the following options [default value]:"
            echo "you can use either - or -- flags"
            echo "-f | --full (will run full analysis, including short read preprocessing)"
            echo "-p | --partial (will skip preprocessing steps)"
            echo "-a | --transcripts (use transcripts.gtf for gene models and skip preprocessing)"
            echo "-r | --mate_inner_distance [165] (expected mean inner distance between mate pairs (PE only))"
            echo "-i | --min_intron_length [50] (minimum intron length)"
            echo "-I | --max_intron_length [25000] (maximum intron length)"
            echo "-t | --agg_transcripts (generate gtf file of empirical transcripts)"
            echo "-s | --refseq [refseq] (name of file containing reference DNA seqeunce)"
            echo "-H | --threads [4] (number of threads to use)"
            echo "-l | --library_type [fr-unstranded] (library type as defined in TopHat manual)"
            #echo "--use_aggregates"
            echo "-e | --seonly (use if NOT working with paired-end sequence data)"
            echo "-A | --adapter (provide the adapter sequence to remove)"
            echo "-P | --indexpath [index] (provide the path to the directory containing the bowtie indexes for refseq and filter)"
            echo "-T | --toolpath [.] (provide path to directory containing RNAseq tools)"
            echo "-h | --help [print this help message]"
            echo "" ;;

        Darwin)
            echo ""
            echo "invoke script with the following options [default value]:"
            echo "-f (will run full analysis, including short read preprocessing)"
            echo "-p (will skip preprocessing steps)"
            echo "-a (use transcripts.gtf for gene models and skip preprocessing)"
            echo "-r [165] (expected mean inner distance between mate pairs (PE only))"
            echo "-i [50] (minimum intron length)"
            echo "-I [25000] (maximum intron length)"
            echo "-t (generate gtf file of empirical transcripts)"
            echo "-s [refseq] (name of file containing reference DNA seqeunce)"
            echo "-H [4] (number of threads to use)"
            echo "-l [fr-unstranded] (library type as defined in TopHat manual)"
            #echo "--use_aggregates"
            echo "-e (use if NOT working with paired-end sequence data)"
            echo "-A (provide the adapter sequence to remove)"
            echo "-P [index] (provide the path to the directory containing the bowtie indexes for refseq and filter)"
            echo "-T [.] (provide path to directory containing RNAseq tools)"
            echo "-h [print this help message]"
            echo "" ;;

        *)
            echo ""
            echo "invoke script with the following options [default value]:"
            echo "you can use either - or -- flags"
            echo "-f | --full (will run full analysis, including short read preprocessing)"
            echo "-p | --partial (will skip preprocessing steps)"
            echo "-a | --transcripts (use transcripts.gtf for gene models and skip preprocessing)"
            echo "-r | --mate_inner_distance [165] (expected mean inner distance between mate pairs (PE only))"
            echo "-i | --min_intron_length [50] (minimum intron length)"
            echo "-I | --max_intron_length [25000] (maximum intron length)"
            echo "-t | --agg_transcripts (generate gtf file of empirical transcripts)"
            echo "-s | --refseq [refseq] (name of file containing reference DNA seqeunce)"
            echo "-H | --threads [4] (number of threads to use)"
            echo "-l | --library_type [fr-unstranded] (library type as defined in TopHat manual)"
            #echo "--use_aggregates"
            echo "-e | --seonly (use if NOT working with paired-end sequence data)"
            echo "-A | --adapter (provide the adapter sequence to remove)"
            echo "-P | --indexpath [index] (provide the path to the directory containing the bowtie indexes for refseq and filter)"
            echo "-T | --toolpath [.] (provide path to directory containing RNAseq tools)"
            echo "-h | --help [print this help message]"
            echo "" ;;

    esac
}

function mk_agg_txpts () {
    cd $wd
    echo "generating transcripts file"
    mkdir -p transcripts
    cd transcripts
    echo "cuffcompare -s ../index/$refseq.fa ../*/cufflinks/transcripts.gtf"
    cuffcompare -s ../index/$refseq.fa ../*/cufflinks/transcripts.gtf
#    samtools merge all_merged.bam ../*/merged/merged.bam
#    echo "cufflinks -p $threads -N --library-type $library_type -I 25000 -L allmerge -r ../index/$refseq.fa all_merged.bam"
#    cufflinks -p $threads -N --library-type $library_type -I 25000 -L allmerge -r ../index/$refseq.fa all_merged.bam
    cd ..
    ln -sf transcripts/stdout.combined.gtf ./transcripts.gtf
}

#function mk_agg_txpts () {
#    cd $wd
#    echo "generating aggregate transcripts file"
#    mkdir -p merged_aggregates
#    cd merged_aggregates
#    samtools merge all_merged.bam ../*/merged/merged.bam
#    echo "cufflinks -p $threads -N --library-type $library_type -I 25000 -L allmerge -r ../index/$refseq.fa all_merged.bam"
#    cufflinks -p $threads -N --library-type $library_type -I 25000 -L allmerge -r ../index/$refseq.fa all_merged.bam
#    cd ..
#    ln -sf merged_aggregates/transcripts.gtf ./
#}

#function mk_agg_jncts {
#
#}

run_type='NULL'
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
indexpath="$wd/index/"

# edit this variable to be the path to RNAseq toolkit an you won't need to use the --toolpath command line flag
toolpath='.'

# command line option parsing adpated from /usr/share/doc/util-linux-2.13/getopt-parse.bash
#
case "$osname" in

    Linux)
        TEMP=`getopt -o pafhr:i:I:jts:H:l:ueA:P:T: --long help,full,transcripts,partial,mate_inner_distance:,min_intron_length:,max_intron_length:,agg_junctions,agg_transcripts,refseq:,threads:,library_type:,use_aggregates,seonly,adapter:,indexpath:,toolpath: -- "$@"`
        ;;

    Darwin)
        TEMP=`getopt pafhr:i:I:jts:H:l:ueA:P:T: $*`
        ;;

    *)
        echo "something else"
        ;;
esac

if [ $? != 0 ] ; then echo "Terminating..." ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!

eval set -- "$TEMP"
while true ; do
    case "$1" in
        -f|--full) run_type='full' ; shift ;;
        -p|--partial) run_type='partial' ; shift ;;
        -a|--transcripts) run_type='transcripts' ; shift ;;
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
        -T|--toolpath) toolpath=$2 ; shift 2 ;;
        -h|--help) help_messg ; exit ;;
        --) shift ; break ;;
        *) break ;;
    esac
done
#for arg do echo '--> '"\`$arg'" ; done

echo "run type is '$run_type'"
#flags="--refseq $refseq --mate_inner_distance $mate_inner_distance --min_intron_length $min_intron_length --max_intron_length $max_intron_length --procs $threads --librarytype $library_type --indexpath $indexpath/" # change to work with MacOSX, below
flags="-s $refseq -r $mate_inner_distance -i $min_intron_length -I $max_intron_length -t $threads -l $library_type -P $indexpath/"
#echo "flags: $flags"
#exit

if [[ $seonly -eq 1 ]]
then
#    flags="$flags --seonly"
    flags="$flags -e"
fi

if [[ $adapter != "NULL" ]]
then
#    flags="$flags --adapter_seq $adapter"    
    flags="$flags -A $adapter"    
fi

#echo "flags: '$flags'"
#exit

for dir 
do
	echo `date`
	#echo $wd/$dir
	#cd $wd/$dir
    echo $dir
    cd $dir

    case "$run_type" in

        partial)
	
            # for partial runs (when you don't need to run preprocessing steps
#            echo "RNAseq.sh $flags --partial" ;
            echo "RNAseq.sh $flags -p" ;
#            RNAseq.sh $flags --partial ;;
            RNAseq.sh $flags -p ;;

        full)

            # for full runs
#            echo "RNAseq.sh $flags --full" ;
            echo "RNAseq.sh $flags -f" ;
#            RNAseq.sh $flags --full ;;
            RNAseq.sh $flags -f ;;

        transcripts)

            # for transcripts runs
            # must have already run script with --agg_transcripts flag
            mkdir -p non-aggregate
            echo "moving old output files to 'non-aggregate'"
            mv -f merged cufflinks pe_tophat* singles_tophat* non-aggregate/
            echo "creating symbolic link to transcript.gtf"
#            ln -sf ../transcripts.gtf ./
            ln -sf $wd/transcripts.gtf ./
            more_flags=""
            echo "running tophat"
#            echo "RNAseq.sh --transcripts $flags $more_flags" ;
            echo "RNAseq.sh -a $flags $more_flags" ;
#            RNAseq.sh --transcripts $flags $more_flags ;;
            RNAseq.sh -a $flags $more_flags ;;

        *)
            echo "RNAeq.sh $flags"
            RNAseq.sh $flags ;;

    esac
    cd $wd
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

