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
            echo "-m | --splice_mismatches [0] (max number of mismatches in anchor region of spliced alignment)"
            echo "-c | --min_anchor_length [8] (minimum number of reads on each side of splice junction)"
            echo "-S | --mate_std_dev [20] (std dev for inner distances between mate pairs)"
            echo "-F | --min-isoform-fraction [0.15] (min isoform fraction)"
            echo "-g | --max_multihits [20] (max multihits)"
            echo "-v | --coverage_search (enables coverage search)"
            echo "-b | --butterfly_search (enables butterfly search)"
            echo "-L | --segment-length [25] (segment length)"
            echo "-M | --segment-mismatches [2] (segment mismatches [0-3])"
            echo "-t | --agg_transcripts (generate gtf file of empirical transcripts)"
            echo "-s | --refseq [refseq] (name of file containing reference DNA seqeunce)"
            echo "-H | --threads [8] (number of threads to use)"
            echo "-l | --library_type [fr-unstranded] (library type as defined in TopHat manual)"
            #echo "--use_aggregates"
            echo "-e | --seonly (use if NOT working with paired-end sequence data)"
            echo "-A | --adapter (provide the adapter sequence to remove)"
            echo "-P | --indexpath [index] (provide the path to the directory containing the bowtie indexes for refseq and filter)"
            echo "-T | --toolpath [.] (provide path to directory containing RNAseq tools)"
            echo "-R | --preprocess (use if you want to ensure running sequence preprocessing routines)"
            echo "-q | --min_qual [13] during preprocessing, minimum quality of base to avoid trimming"
            echo "-n | --min_length [32] during preprocessing, minimum acceptable length after trimming"
            echo "-E | --percent_high_quality [90] during preprocessing, minimum percentage of bases >= min_qual"
            echo "-Q | --solexa Solexa quality scores. Current default is Illumina (Phred-scaled base-64). Illumina switched back to Solexa in 2011."
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
            echo "-m [0] (max number of mismatches in anchor region of spliced alignment)"
            echo "-c [8] (minimum number of reads on each side of splice junction)"
            echo "-S [20] (std dev for inner distances between mate pairs)"
            echo "-F [0.15] (min isoform fraction)"
            echo "-g [20] (max multihits)"
            echo "-v (enables coverage search)"
            echo "-b (enables butterfly search)"
            echo "-L [25] (segment length)"
            echo "-M [2] (segment mismatches [0-3])"
            echo "-t (generate gtf file of empirical transcripts)"
            echo "-s [refseq] (name of file containing reference DNA seqeunce)"
            echo "-H [8] (number of threads to use)"
            echo "-l [fr-unstranded] (library type as defined in TopHat manual)"
            #echo "--use_aggregates"
            echo "-e (use if NOT working with paired-end sequence data)"
            echo "-A (provide the adapter sequence to remove)"
            echo "-P [index] (provide the path to the directory containing the bowtie indexes for refseq and filter)"
            echo "-T [.] (provide path to directory containing RNAseq tools)"
            echo "-R (use if you want to ensure running sequence preprocessing routines)"
            echo "-q [13] during preprocessing, minimum quality of base to avoid trimming"
            echo "-n [32] during preprocessing, minimum acceptable length after trimming"
            echo "-E [90] during preprocessing, minimum percentage of bases >= min_qual"
            echo "-Q Solexa quality scores. Current default is Illumina (Phred-scaled base-64). Illumina switched back to Solexa in 2011."
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
            echo "-m | --splice_mismatches [0] (max number of mismatches in anchor region of spliced alignment)"
            echo "-c | --min_anchor_length [8] (minimum number of reads on each side of splice junction)"
            echo "-S | --mate_std_dev [20] (std dev for inner distances between mate pairs)"
            echo "-g | --max_multihits [20] (max multihits)"
            echo "-v | --coverage_search (enables coverage search)"
            echo "-b | --butterfly_search (enables butterfly search)"
            echo "-L | --segment-length [25] (segment length)"
            echo "-M | --segment-mismatches [2] (segment mismatches [0-3])"
            echo "-F | --min-isoform-fraction [0.15] (min isoform fraction)"
            echo "-t | --agg_transcripts (generate gtf file of empirical transcripts)"
            echo "-s | --refseq [refseq] (name of file containing reference DNA seqeunce)"
            echo "-H | --threads [8] (number of threads to use)"
            echo "-l | --library_type [fr-unstranded] (library type as defined in TopHat manual)"
            #echo "--use_aggregates"
            echo "-e | --seonly (use if NOT working with paired-end sequence data)"
            echo "-A | --adapter (provide the adapter sequence to remove)"
            echo "-P | --indexpath [index] (provide the path to the directory containing the bowtie indexes for refseq and filter)"
            echo "-T | --toolpath [.] (provide path to directory containing RNAseq tools)"
            echo "-R | --preprocess (use if you want to ensure running sequence preprocessing routines)"
            echo "-q | --min_qual [13] during preprocessing, minimum quality of base to avoid trimming"
            echo "-n | --min_length [32] during preprocessing, minimum acceptable length after trimming"
            echo "-E | --percent_high_quality [90] during preprocessing, minimum percentage of bases >= min_qual"
            echo "-Q | --solexa Solexa quality scores. Current default is Illumina (Phred-scaled base-64). Illumina switched back to Solexa in 2011."
            echo "-h | --help [print this help message]"
            echo "" ;;

    esac
}

function mk_agg_txpts () {
    cd $wd
    echo "generating transcripts file"
    mkdir -p transcripts
    cd transcripts
    echo "cuffcompare -s $wd/index/$refseq.fa ../*/cufflinks/transcripts.gtf"
    cuffcompare -s $wd/index/$refseq.fa ../*/cufflinks/transcripts.gtf
#    samtools merge all_merged.bam ../*/merged/merged.bam
#    echo "cufflinks -p $threads -N --library-type $library_type -I 25000 -L allmerge -r ../index/$refseq.fa all_merged.bam"
#    cufflinks -p $threads -N --library-type $library_type -I 25000 -L allmerge -r ../index/$refseq.fa all_merged.bam
    cd ..
    #ln -sf transcripts/stdout.combined.gtf ./transcripts.gtf
    # name change means the link above doesn't work
    ln -sf transcripts/cuffcmp.combined.gtf ./transcripts.gtf
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
preprocess=0
splice_mismatches=0
min_anchor_length=8
mate_std_dev=20
min_isoform_frac=0.15
max_multihits=20
coverage_search=0
butterfly_search=0
segment_length=25
segment_mismatches=2
min_qual=13
min_length=32
percent_high_quality=90
qualscores='NULL'

# edit this variable to be the path to RNAseq toolkit an you won't need to use the --toolpath command line flag
toolpath='.'

# command line option parsing adpated from /usr/share/doc/util-linux-2.13/getopt-parse.bash
#
case "$osname" in

    Linux)
        TEMP=`getopt -o pafhr:i:I:jts:H:l:ueA:P:T:Rm:c:S:F:g:vbL:M:q:n:E:Q --long help,full,transcripts,partial,mate_inner_distance:,min_intron_length:,max_intron_length:,agg_junctions,agg_transcripts,refseq:,threads:,library_type:,use_aggregates,seonly,adapter:,indexpath:,toolpath:,preprocess,splice_mismatches:,min_anchor_length:,mate_std_dev:,min_isoform_fraction:,max_multihits:,coverage_search,butterfly_search,segment_length:,segment_mismatches:,min_qual:,min_length:,percent_high_quality:,solexa -- "$@"`
        ;;

    Darwin)
        TEMP=`getopt pafhr:i:I:jts:H:l:ueA:P:T:Rm:c:S:F:g:vbL:M:q:n:E:Q $*`
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
        -m|--splice_mismatches) splice_mismatches=$2 ; shift 2 ;;
        -c|--min_anchor_length) min_anchor_length=$2 ; shift 2 ;;
        -S|--mate_std_dev) mate_std_dev=$2 ; shift 2 ;;
        -F|--min_isoform_fraction) min_isoform_frac=$2 ; shift 2 ;;
        -g|--max_multihits) max_multihits=$2 ; shift 2 ;;
        -v|--coverage_search) coverage_search=1 ; shift ;;
        -b|--butterfly_search) butterfly_search=1 ; shift ;;
        -L|--segment_length) segment_length=$2 ; shift 2 ;;
        -M|--segment_mismatches) segment_mismatches=$2 ; shift 2 ;; 
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
        -R|--preprocess) preprocess=1 ; shift ;;
        -q|--min_qual) min_qual=$2 ; shift 2 ;;
        -n|--min_length) min_length=$2 ; shift 2 ;;
        -E|--percent_high_quality) percent_high_quality=$2 ; shift 2 ;;
        -Q|--solexa) qualscores=1 ; shift ;;
        -h|--help) help_messg ; exit ;;
        --) shift ; break ;;
        *) break ;;
    esac
done
#for arg do echo '--> '"\`$arg'" ; done

echo "run type is '$run_type'"
#flags="-s $refseq -r $mate_inner_distance -i $min_intron_length -I $max_intron_length -t $threads -l $library_type -P $indexpath/ -m $splice_mismatches -c $min_anchor_length -S $mate_std_dev -F $min_isoform_frac -g $max_multihits -L $segment_length -M $segment_mismatches -q $min_qual -n $min_length -E $percent_high_quality"
flags="-s $refseq -r $mate_inner_distance -i $min_intron_length -I $max_intron_length -t $threads -l $library_type -P $indexpath/ -m $splice_mismatches -c $min_anchor_length -S $mate_std_dev -F $min_isoform_frac -g $max_multihits -L $segment_length -M $segment_mismatches -q $min_qual -n $min_length -E $percent_high_quality"
#echo "flags: $flags"
#exit

if [[ $run_type = "full" ]]
then
    preprocess=1
fi

if [[ $seonly -eq 1 ]]
then
#    flags="$flags --seonly"
    flags="$flags -e"
fi

if [[ $adapter != "NULL" ]]
then
    flags="$flags -A $adapter"    
fi

if [[ $qualscores != "NULL" ]]
then
    if [[ $qualscores -eq 1 ]]
    then
        flags="$flags -Q"
    fi
fi

if [[ $coverage_search -ne 0 ]]
then
    flags="$flags -v"
fi
if [[ $butterfly_search -ne 0 ]]
then
    flags="$flags -b"
fi
echo "flags: '$flags'"
#exit

#echo "preprocess = " $preprocess

if [[ $preprocess -ne 0 ]]
then
    more_flags="-R"
else
    more_flags=""
fi
#echo "preprocess = " $preprocess
#echo "flags = " $flags
#echo "more_flags = " $more_flags
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
            RNAseq.sh $flags -p $more_flags ;;

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

#            if [[ preprocess != "NULL" ]]
#            then
#                more_flags="-R"
#            else
#                more_flags=""
#            fi

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

