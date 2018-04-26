#!/bin/bash
#
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
wd=`pwd`
osname=`uname -s`
#export PATH=# <-- make sure RNAseq scripts are in your path
export PATH=".:$wd:$wd/bin:$HOME/bin:$PATH"
script=`which RNAseq.sh`
echo "script \"$script\""
#

function help_messg () {
    case "$osname" in

        Linux) 
            echo "invoke as RNAseq_process_data.sh [options] set1.fq set2.fq"
            echo "invoke script with the following options [default value]:"
            echo "you can use either - or -- flags"
            echo "-f | --full (will run full analysis, including short read preprocessing)"
            echo "-p | --partial (will skip preprocessing steps)"
            echo "-O | --preprocess_only (only run preprocessing routines)"
            echo "-B | --submit submit to Slurm job manager" 
            echo "-D | --partition Slurm partition [CLUSTER]"
            echo "-w | --wait wait for first job submitted to Slurm job manager to finish before running other jobs"
            echo "-i | --min_intron_length [20] (minimum intron length)"
            echo "-I | --max_intron_length [500000] (maximum intron length)"
#            echo "-M | --segment-mismatches [2] (segment mismatches [0-3])"
            echo "-H | --threads [8] (number of threads to use)"
            echo "-l | --library_type [fr-unstranded] (library type as defined in TopHat manual)"
            #echo "--use_aggregates"
            echo "-e | --seonly (use if NOT working with paired-end sequence data)"
            echo "-A | --adapter (provide the adapter sequence to remove)"
            echo "-P | --indexpath [index] (provide the path to the directory containing the STAR indexes for ref seq and bowtie index for filter.fa)"
            echo "-T | --toolpath [.] (provide path to directory containing RNAseq tools)"
            echo "-R | --preprocess (use if you want to ensure running sequence preprocessing routines)"
            echo "-q | --min_qual [13] during preprocessing, minimum quality of base to avoid trimming"
            echo "-n | --min_length [32] during preprocessing, minimum acceptable length after trimming"
            echo "-E | --percent_high_quality [90] during preprocessing, minimum percentage of bases >= min_qual"
            echo "-Q | --solexa Solexa quality scores (Phred-33)"
            echo "-X | --phred33 Phred quality values encoded as Phred + 33"
            echo "-Y | --phred64 Phred quality values encoded as Phred + 64"
            echo "-C | --leave_temp leave temporary files on file system"
            echo "-F | --nofilter do not do sequence similarity filtering"
            echo "-x | --memory amount of memory to allocate [50G]"
#            echo "-J | --ignore_single_exons ignore single exon transfrags (& reference transcripts) when combining from multiple GTF files"
            echo "-h | --help [print this help message]"
            echo "" ;;

        Darwin)
            echo ""
            echo "invoke script with the following options [default value]:"
            echo "-f (will run full analysis, including short read preprocessing)"
            echo "-p (will skip preprocessing steps)"
            echo "-B submit to Slurm job manager" 
            echo "-D Slurm Partition [CLUSTER]"
            echo "-w wait for first job submitted to Slurm job manager to finish before running other jobs"
            echo "-r [165] (expected mean inner distance between mate pairs (PE only))"
            echo "-i [50] (minimum intron length)"
            echo "-I [25000] (maximum intron length)"
            echo "-m [0] (max number of mismatches in anchor region of spliced alignment)"
            echo "-c [8] (minimum number of reads on each side of splice junction)"
            echo "-S [20] (std dev for inner distances between mate pairs)"
            echo "-g [20] (max multihits)"
            echo "-G [2] (max mismatches)"
            echo "-v (enables coverage search)"
            echo "-M [2] (segment mismatches [0-3])"
            echo "-t (generate gtf file of empirical transcripts)"
            echo "-K generate output files compatible with cufflinks"
            echo "-H [8] (number of threads to use)"
            echo "-l [fr-unstranded] (library type as defined in TopHat manual)"
            #echo "--use_aggregates"
            echo "-e (use if NOT working with paired-end sequence data)"
            echo "-A (provide the adapter sequence to remove)"
            echo "-P [index] (provide the path to the directory containing the bowtie indexes for refseq.fa and filter.fa)"
            echo "-T [.] (provide path to directory containing RNAseq tools)"
            echo "-R (use if you want to ensure running sequence preprocessing routines)"
            echo "-q [13] during preprocessing, minimum quality of base to avoid trimming"
            echo "-n [32] during preprocessing, minimum acceptable length after trimming"
            echo "-E [90] during preprocessing, minimum percentage of bases >= min_qual"
            echo "-Q solexa Solexa quality scores (Phred-33)"
            echo "-X Phred quality values encoded as Phred + 33"
            echo "-Y Phred quality values encoded as Phred + 64"
            echo "-C | --leave_temp leave temporary files on file system"
            echo "-F | --nofilter do not do sequence similarity filtering"
            echo "-x | --memory amount of memory to allocate [50G]"
#            echo "-J | --ignore_single_exons ignore single exon transfrags (& reference transcripts) when combining from multiple GTF files"
            echo "-h [print this help message]"
            echo "" ;;

        *)
            echo ""
            echo "invoke script with the following options [default value]:"
            echo "you can use either - or -- flags"
            echo "-f | --full (will run full analysis, including short read preprocessing)"
            echo "-p | --partial (will skip preprocessing steps)"
            echo "-B | --submit submit to Slurm job manager" 
            echo "-D | --partition Slurm partition [CLUSTER]"
            echo "-w | --wait wait for first job submitted to Slurm job manager to finish before running other jobs"
            echo "-i | --min_intron_length [20] (minimum intron length)"
            echo "-I | --max_intron_length [500000] (maximum intron length)"
#            echo "-M | --segment-mismatches [2] (segment mismatches [0-3])"
            echo "-H | --threads [8] (number of threads to use)"
            echo "-l | --library_type [fr-unstranded] (library type as defined in TopHat manual)"
            #echo "--use_aggregates"
            echo "-e | --seonly (use if NOT working with paired-end sequence data)"
            echo "-A | --adapter (provide the adapter sequence to remove)"
            echo "-P | --indexpath [index] (provide the path to the directory containing the STAR indexes for ref seq and bowtie index for filter.fa)"
            echo "-T | --toolpath [.] (provide path to directory containing RNAseq tools)"
            echo "-R | --preprocess (use if you want to ensure running sequence preprocessing routines)"
            echo "-q | --min_qual [13] during preprocessing, minimum quality of base to avoid trimming"
            echo "-n | --min_length [32] during preprocessing, minimum acceptable length after trimming"
            echo "-E | --percent_high_quality [90] during preprocessing, minimum percentage of bases >= min_qual"
            echo "-Q | --solexa Solexa quality scores (Phred-33)"
            echo "-X | --phred33 Phred quality values encoded as Phred + 33"
            echo "-Y | --phred64 Phred quality values encoded as Phred + 64"
            echo "-C | --leave_temp leave temporary files on file system"
            echo "-F | --nofilter do not do sequence similarity filtering"
            echo "-x | --memory amount of memory to allocate [50G]"
#            echo "-J | --ignore_single_exons ignore single exon transfrags (& reference transcripts) when combining from multiple GTF files"
            echo "-h | --help [print this help message]"
            echo "" ;;

    esac
}

run_type='NULL'
min_intron_length=20
max_intron_length=500000
aggregate_junctions=0
aggregate_transcripts=0
threads=8
library_type='NULL'
use_aggregates=0
seonly=0
adapter='NULL'
indexpath="$wd/index/"
preprocess=0
preprocess_only=0
min_qual=13
min_length=32
percent_high_quality=90
# default qualscores should be Phred+33
qualscores='NULL'
dev=0
oldid=0
RNAseq_script='NULL'
bsub=0
queue='CLUSTER'
no_new_txpts='NULL'
leave_temp=0
ignore_single_exons=0
wait4first=0
run_STAR=1
cufflinks_compatible=0
nofilter=0
memory='50G'

# edit this variable to be the path to RNAseq toolkit an you won't need to use the --toolpath command line flag
toolpath='.'

# command line option parsing adpated from /usr/share/doc/util-linux-2.13/getopt-parse.bash
#
case "$osname" in

    Linux)
        TEMP=`getopt -o pfhr:i:I:jH:l:ueA:P:T:ROm:c:S:g:vbM:q:n:E:QdCNXx:YoG:B:wJF --long help,full,partial,min_intron_length:,max_intron_length:,agg_junctions,threads:,library_type:,use_aggregates,seonly,adapter:,indexpath:,toolpath:,preprocess,preprocess_only,min_qual:,min_length:,percent_high_quality:,solexa,dev,leave_temp,oldid,phred33,phred64,submit,partition:,wait,ignore_single_exons,nofilter,memory: -- "$@"`
        ;;

    Darwin)
        TEMP=`getopt pfhr:i:I:jH:l:ueA:P:T:Rm:c:S:g:vbM:q:n:E:QdCNXx:YoG:B:wJF $*`
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
        -i|--min_intron_length) min_intron_length=$2 ; shift 2 ;;
        -I|--max_intron_length) max_intron_length=$2 ; shift 2 ;;
        -j|--agg_junctions) aggregate_junctions=1 ; shift ;;
        -H|--threads) threads=$2 ; shift 2 ;;
        -l|--library_type) library_type=$2 ; shift 2 ;;
        -u|--use_aggregates) use_aggregates=1 ; shift ;;
        -e|--seonly) seonly=1 ; shift ;;
        -A|--adapter) adapter=$2 ; shift 2 ;;
        -P|--indexpath) indexpath=$2 ; shift 2 ;;
        -T|--toolpath) toolpath=$2 ; shift 2 ;;
        -R|--preprocess) preprocess=1 ; shift ;;
        -O|--preprocess_only) run_type='preprocess' ; shift ;;
        -q|--min_qual) min_qual=$2 ; shift 2 ;;
        -n|--min_length) min_length=$2 ; shift 2 ;;
        -E|--percent_high_quality) percent_high_quality=$2 ; shift 2 ;;
        -Q|--solexa) qualscores=3 ; shift ;;
        -X|--phred33) qualscores=1 ; shift ;;
        -Y|--phred64) qualscores=2 ; shift ;;
        -h|--help) help_messg ; exit ;;
        -d|--dev) dev=1 ; shift ;;
        -C|--leave_temp) leave_temp=1 ; shift ;;
        -J|--ignore_single_exons) ignore_single_exons=1 ; shift ;;
        -N|--oldid) oldid=1 ; shift ;;
        -B|--submit) bsub=1 ; shift ;;
        -D|--partition) queue=$2 ; shift 2 ;;
        -w|--wait) wait4first=1 ; shift ;;
        -F|--nofilter) nofilter=1 ; shift ;;
        -x|--memory) memory=$2 ; shift 2 ;;
        --) shift ; break ;;
        *) break ;;
    esac
done
for arg do echo '--> '"\`$arg'" ; done

echo "run type is '$run_type'"

flags="-i $min_intron_length -I $max_intron_length -t $threads -P $indexpath/ -q $min_qual -n $min_length -E $percent_high_quality" 
#echo "flags: $flags"
#exit

if [[ $library_type -ne 'NULL' ]]
then
    # Illumina TruSeq stranded should be fr-firststrand
    flags="$flags --rna-strandness $library_type"
fi

if [[ $run_type = "full" ]]
then
    preprocess=1
fi

if [[ $seonly -eq 1 ]]
then
#    flags="$flags --seonly"
    flags="$flags -e"
fi

echo "RNAseq_script = '$RNAseq_script'"

if [[ $adapter != "NULL" ]]
then
    flags="$flags -A $adapter"    
fi

# default quality scores should be Phred+33
if [[ $qualscores != "NULL" ]]
then
    if [[ $qualscores -eq 2 ]]
    then
        flags="$flags -Q"
    elif [[ $qualscores -eq 3 ]]
    then
        flags="$flags -B"
    fi
fi

if [[ $oldid -ne 0 ]]
then
    flags="$flags -N"
fi

if [[ $leave_temp -ne 0 ]]
then
    flags="$flags -C"
fi

if [[ $nofilter -ne 0 ]]
then
    flags="$flags -F"
fi

#echo "flags: '$flags'"
#exit

#echo "preprocess = " $preprocess

more_flags=""
if [[ $preprocess -ne 0 ]]
then
    more_flags="-R"
#else
#    more_flags=""
fi

if [[ $preprocess_only -ne 0 ]]
then
    preprocess=1
    more_flags="-R -O"
fi

if [[ $dev -ne 0 ]]
then
    # this works in-house at IRCF
    export PATH="/ircf/ircfapps/dev/bin:$PATH"
    export PATH="/home/sgivan/projects/RNA-Seq-Toolkit/bin:$PATH"
    flags="$flags -d"
fi

#echo "preprocess = " $preprocess
echo "flags = " $flags
#echo "more_flags = " $more_flags
#exit

let "cnt = 0"
let "slurm_id = 0"

echo "entering loop"
for dir 
do
    let "++cnt"
    echo "cnt = '$cnt'"
	echo `date`
	#echo $wd/$dir
	#cd $wd/$dir
    echo $dir
    cd $dir






    RNAseq_script="$script"
    echo "RNAseq_script set to'${script}'"

    batch=0
    case "$run_type" in

        partial)
            echo "run type still 'partial'"	
            # for partial runs (when you don't need to run preprocessing steps
#            echo "RNAseq.sh $flags --partial" ;
            echo "$RNAseq_script $flags -p $more_flags" ;
            $($RNAseq_script $flags -p $more_flags) ;
#            echo "$RNAseq_script $flags -p $more_flags " > cmd ;
            batch=1;;
            #eval $RNAseq_script $flags -p $more_flags > RNAseq.log 2>&1 ;;

        full)

            # for full runs
#            echo "RNAseq.sh $flags --full" ;
            echo "$RNAseq_script $flags -f" ;
#            RNAseq.sh $flags --full ;;
            echo "$RNAseq_script $flags -f " > cmd ;
            batch=1;;
            #eval $RNAseq_script $flags -f > RNAseq.log 2>&1 ;;

        preprocess)

            echo "running preprocessing routines only" ;
            echo "$RNAseq_script -R -O $flags $more_flags" ;
            echo "$RNAseq_script -R -O $flags $more_flags " > cmd ;;

        *)
            echo "RNAeq.sh $flags"
            echo "$RNAseq_script $flags $more_flags " > cmd ;;

    esac

    if [[ $bsub != 0 ]]
    then
        if [[ $wait4first -eq 1 && $cnt -ne 1 ]]
        then
            echo "subsequent jobs will wait for job JOB ID: '"$slurm_id"'"
            if [[ $batch -eq 1 ]]
            then
                OUTPUT="$(sbatch -J $dir cmd)"
            else
                OUTPUT="$(sbatch --mem=${memory} --depend=afterok:${slurm_id} -o ./${dir}.o -e ./${dir}.e -J $dir  --partition $queue --ntasks=1 --cpus-per-task $threads   --wrap='sh cmd')"
            fi

        else
            if [[ $batch -eq 1 ]]
            then
                OUTPUT="$(sbatch -J ${dir} cmd)"
                slurm_id=$(echo $OUTPUT | sed 's/Submitted batch job //')
                #echo "subsequent jobs will wait for job JOB ID: '"$slurm_id"'"
            else

                OUTPUT="$(sbatch --mem=${memory} -o ./${dir}.o -e ./${dir}.e -J $dir  --partition $queue --ntasks=1 --cpus-per-task $threads --wrap='sh cmd')"
                slurm_id=$(echo $OUTPUT | sed 's/Submitted batch job //')
                #echo "subsequent jobs will wait for job JOB ID: '"$slurm_id"'"
            fi
        fi
    fi
    echo "OUTPUT: "$OUTPUT 
#    slurm_id=$(echo $OUTPUT | sed 's/Submitted batch job //')
#    echo "JOB ID: '"$slurm_id"'"

    RETVAL=$?
    if [[ $RETVAL -eq 0 ]]
    then
        echo "success!"
    else
        echo "non-zero exit value"
        exit
    fi

    cd $wd
done

