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
            echo "-a | --transcripts (use transcripts.gtf for gene models and skip preprocessing)"
            echo "-k | --nonewtranscripts (if using --transcripts, only align to known transcripts)"
            echo "-B | --slurm submit to slurm queueing system" 
            echo "-D | --queue LSF queue [normal]"
            echo "-w | --wait wait for first job submitted to LSF queue finish before running other jobs"
            echo "-i | --min_intron_length [20] (minimum intron length)"
            echo "-I | --max_intron_length [500000] (maximum intron length)"
            echo "-F | --min-isoform-fraction [0.15] (min isoform fraction)"
            echo "-L | --segment-length [25] (segment length)"
#            echo "-M | --segment-mismatches [2] (segment mismatches [0-3])"
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
            echo "-Q | --solexa Solexa quality scores (Phred-33)"
            echo "-X | --phred33 Phred quality values encoded as Phred + 33"
            echo "-Y | --phred64 Phred quality values encoded as Phred + 64"
            echo "-N | --oldid fastq files use old sequence/quality ID structure"
            echo "-C | --leave_temp leave temporary files on file system"
            echo "-J | --ignore_single_exons ignore single exon transfrags (& reference transcripts) when combining from multiple GTF files"
            echo "-h | --help [print this help message]"
            echo "--no_hisat don't run hisat"
            echo "" ;;

        Darwin)
            echo ""
            echo "invoke script with the following options [default value]:"
            echo "-f (will run full analysis, including short read preprocessing)"
            echo "-p (will skip preprocessing steps)"
            echo "-a (use transcripts.gtf for gene models and skip preprocessing)"
            echo "-B submit to LSF queueing system" 
            echo "-D LSF queue [normal]"
            echo "-w wait for first job submitted to LSF queue finish before running other jobs"
            echo "-r [165] (expected mean inner distance between mate pairs (PE only))"
            echo "-i [50] (minimum intron length)"
            echo "-I [25000] (maximum intron length)"
            echo "-m [0] (max number of mismatches in anchor region of spliced alignment)"
            echo "-c [8] (minimum number of reads on each side of splice junction)"
            echo "-S [20] (std dev for inner distances between mate pairs)"
            echo "-F [0.15] (min isoform fraction)"
            echo "-g [20] (max multihits)"
            echo "-G [2] (max mismatches)"
            echo "-v (enables coverage search)"
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
            echo "-Q solexa Solexa quality scores (Phred-33)"
            echo "-X Phred quality values encoded as Phred + 33"
            echo "-Y Phred quality values encoded as Phred + 64"
            echo "-C | --leave_temp leave temporary files on file system"
            echo "-J | --ignore_single_exons ignore single exon transfrags (& reference transcripts) when combining from multiple GTF files"
            echo "-h [print this help message]"
            echo "" ;;

        *)
            echo ""
            echo "invoke script with the following options [default value]:"
            echo "you can use either - or -- flags"
            echo "-f | --full (will run full analysis, including short read preprocessing)"
            echo "-p | --partial (will skip preprocessing steps)"
            echo "-a | --transcripts (use transcripts.gtf for gene models and skip preprocessing)"
            echo "-k | --nonewtranscripts (if using --transcripts, only align to known transcripts)"
            echo "-B | --slurm submit to slurm queueing system" 
            echo "-D | --queue LSF queue [normal]"
            echo "-w | --wait wait for first job submitted to LSF queue finish before running other jobs"
            echo "-i | --min_intron_length [20] (minimum intron length)"
            echo "-I | --max_intron_length [500000] (maximum intron length)"
            echo "-L | --segment-length [25] (segment length)"
#            echo "-M | --segment-mismatches [2] (segment mismatches [0-3])"
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
            echo "-Q | --solexa Solexa quality scores (Phred-33)"
            echo "-X | --phred33 Phred quality values encoded as Phred + 33"
            echo "-Y | --phred64 Phred quality values encoded as Phred + 64"
            echo "-N | --oldid fastq files use old sequence/quality ID structure"
            echo "-C | --leave_temp leave temporary files on file system"
            echo "-J | --ignore_single_exons ignore single exon transfrags (& reference transcripts) when combining from multiple GTF files"
            echo "-h | --help [print this help message]"
            echo "--no_hisat don't run hisat"
            echo "" ;;

    esac
}

function mk_agg_txpts () {
    cd $wd
    echo "generating transcripts file"
    if [[ -e "transcripts" ]] || mkdir -p transcripts
    then
        cd transcripts
        cuffcompare_flags="-s $wd/index/$refseq.fa ../*/cufflinks/transcripts.gtf"
        if [[ $ignore_single_exons -eq 1 ]]
        then
            cuffcompare_flags="-M $cuffcompare_flags"
        fi
#        echo "cuffcompare -s $wd/index/$refseq.fa ../*/cufflinks/transcripts.gtf"
        echo "cuffcompare $cuffcompare_flags"
#        cuffcompare -s $wd/index/$refseq.fa ../*/cufflinks/transcripts.gtf
        cuffcompare $cuffcompare_flags

        cd ..
        #ln -sf transcripts/stdout.combined.gtf ./transcripts.gtf
        # name change means the link above doesn't work
        ln -sf transcripts/cuffcmp.combined.gtf ./transcripts.gtf
    else
        echo "can't create transcripts directory"
    fi
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
min_intron_length=20
max_intron_length=500000
aggregate_junctions=0
aggregate_transcripts=0
refseq='refseq'
threads=8
library_type='NULL'
use_aggregates=0
seonly=0
adapter='NULL'
indexpath="$wd/hisat_index/"
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
slurm=0
queue='normal'
no_new_txpts='NULL'
leave_temp=0
ignore_single_exons=0
wait4first=0
run_hisat=1

# edit this variable to be the path to RNAseq toolkit an you won't need to use the --toolpath command line flag
toolpath='.'

# command line option parsing adpated from /usr/share/doc/util-linux-2.13/getopt-parse.bash
#
case "$osname" in

    Linux)
        TEMP=`getopt -o pafhr:i:I:jts:H:l:ueA:P:T:ROm:c:S:F:g:vbL:M:q:n:E:QdC:NXYoG:BD:kwJ --long help,full,transcripts,partial,min_intron_length:,max_intron_length:,agg_junctions,agg_transcripts,refseq:,threads:,library_type:,use_aggregates,seonly,adapter:,indexpath:,toolpath:,preprocess,preprocess_only,min_qual:,min_length:,percent_high_quality:,solexa,dev,leave_temp,oldid,phred33,phred64,slurm,queue:,nonewtranscripts,wait,ignore_single_exons,no_hisat -- "$@"`
        ;;

    Darwin)
        TEMP=`getopt pafhr:i:I:jts:H:l:ueA:P:T:Rm:c:S:F:g:vbL:M:q:n:E:QdC:NXYoG:BD:kwJ $*`
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
        -k|--nonewtranscripts) no_new_txpts=1 ; shift ;;
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
        -R|--preprocess) preprocess=1 ; shift ;;
#        -O|--preprocess_only) preprocess_only=1 ; shift ;;
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
        -B|--slurm) slurm=1 ; shift ;;
        -D|--queue) queue=$2 ; shift 2 ;;
        -w|--wait) wait4first=1 ; shift ;;
        --no_hisat) run_hisat=0 ; shift ;;
        --) shift ; break ;;
        *) break ;;
    esac
done
#for arg do echo '--> '"\`$arg'" ; done

echo "run type is '$run_type'"
flags="-s $refseq -i $min_intron_length -I $max_intron_length -t $threads -P $indexpath/ -q $min_qual -n $min_length -E $percent_high_quality" 
#echo "flags: $flags"
#exit

if [[ $library_type -ne 'NULL' ]]
then
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

if [[ $run_hisat -eq 0 ]]
then
    flags="$flags --no_hisat"
fi

if [[ $oldid -ne 0 ]]
then
    flags="$flags -N"
fi

if [[ $leave_temp -ne 0 ]]
then
    flags="$flags -C"
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
    flags="$flags -d"
fi

#echo "preprocess = " $preprocess
echo "flags = " $flags
#echo "more_flags = " $more_flags
#exit

let "cnt = 0"
for dir 
do
    let "++cnt"
    echo "cnt = '$cnt'"
	echo `date`
	#echo $wd/$dir
	#cd $wd/$dir
    echo $dir
    cd $dir

    if [[ $slurm != 0 ]]
    then
        if [[ $wait4first -eq 1 && $cnt -eq 1 ]]
        then
            RNAseq_script="sbatch --wait"
        else
            RNAseq_script="sbatch"
        fi

        #RNAseq_script="$RNAseq_script -R \"rusage[mem=1000] span[hosts=1]\" -o ./${dir}.o -e ./${dir}.e -J $dir -q $queue -n $threads $script"
        RNAseq_script="$RNAseq_script -o ./${dir}.%j.o -e ./${dir}.%j.e --mem-per-cpu 12000 --job-name $dir --cpus-per-task $threads --wrap=\"$script"
        echo "6 RNAseq_script: '$RNAseq_script\"'"
        #continue
        #break
        #exit

    else
        RNAseq_script="$script"
    fi

    case "$run_type" in

        partial)
	
            # for partial runs (when you don't need to run preprocessing steps
#            echo "RNAseq.sh $flags --partial" ;
            echo "1 $RNAseq_script $flags -p $more_flags \"" ;
            eval $RNAseq_script $flags -p $more_flags \" ;;
            #eval $RNAseq_script $flags -p $more_flags > RNAseq.log 2>&1 ;;

        full)

            # for full runs
            echo "2 $RNAseq_script $flags -f" ;
#            RNAseq.sh $flags --full ;;
            eval $RNAseq_script $flags -f \" ;;
            #eval $RNAseq_script $flags -f > RNAseq.log 2>&1 ;;

        transcripts)

            # for transcripts runs
            # must have already run script with --agg_transcripts flag
            if [[ $run_hisat -eq 1 ]]
            then
                if [[ -e "non-aggregate" ]] || mkdir -p non-aggregate
                then
                    echo "moving old output files to 'non-aggregate'"
                    mv -f merged cufflinks pe_hisat* singles_hisat* non-aggregate/
                else
                    echo "can't create non-aggregate directory"
                fi
            fi

            # add '-k' flag if needed and not already included
            if [[ $no_new_txps != "NULL" ]] && [[ "$more_flags" != *" -k"* ]]
            then
                more_flags="$more_flags -k"
            fi

            echo "creating symbolic link to transcript.gtf"
            ln -sf $wd/transcripts.gtf ./

            echo "running hisat"
#            echo "RNAseq.sh --transcripts $flags $more_flags" ;
            echo "3 $RNAseq_script -a $flags $more_flags" ;
#            RNAseq.sh --transcripts $flags $more_flags ;;
            eval $RNAseq_script -a $flags $more_flags \" ;;
            #eval $RNAseq_script -a $flags $more_flags > RNAseq.log 2>&1 ;;

        preprocess)

            echo "running preprocessing routines only" ;
            echo "4 $RNAseq_script -R -O $flags $more_flags" ;
            eval $RNAseq_script -R -O $flags $more_flags \" ;;
            #eval $RNAseq_script -R -O $flags $more_flags >> RNAseq.log 2>&1 ;;

        *)
            echo "5 RNAeq.sh $flags"
            eval $RNAseq_script $flags $more_flags \" ;;

    esac

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

if [ $aggregate_junctions = 1 ]
    then
        echo "generating aggregate junctions file"
        if [[ -e "merged_aggregates" ]] || mkdir -p merged_aggregates
        then
            cd merged_aggregates
            cat ../*/*/junctions.bed | awk '{ if ($1 != "track") {split($11,len,","); split($12,blstrt,","); printf "%s\t%i\t%i\t%s\n", $1, $2 + len[1] - 1, $2 + blstrt[2], $6; }}' | sort -k 1,1 -gk 2,2 | uniq > aggregate_junctions.txt
            cd ..
        else
            echo "can't create merged_aggregates directory"
        fi
fi

if [ $aggregate_transcripts = 1 ]
    then
        mk_agg_txpts
fi

if [ $use_aggregates = 1 ]
    then
        echo "re-running pipeline using aggregate files"
fi

