#!/bin/bash

# set path to find RNAseq scripts
wd=`pwd`
#export PATH=# <-- make sure the RST scripts are in your path
export PATH=".:$wd:$wd/bin:$HOME/bin:$PATH"
#
# some variables to set ...
#

cutadapt=$(which cutadapt)
use_cutadapt=1

# HISAT_INDEXES contains the path to the directory containing the reference index files
# built by the hisat-build program. Be sure to end with a "/".
#
HISAT_INDEXES='hisat_index'

#
# some other variables that are typically required by tophat ...
# initialize with reasonable values
#

seonly=0 # for single-end sequence data only -- no paired-end 
procs=8 # number of processors to use
run_type='NULL'
adapter_seq='NULL'
preprocess=0
min_qual=13
min_length=32
percent_high_quality=90
#qualscores='--solexa1.3-quals'
qualscores='NULL'
dev=0
leave_temp=0
oldid=0

osname=`uname -s`
echo "osname '$osname'"

case "$osname" in

    Linux)
            TEMP=`getopt -o et:pafhr:i:I:P:l:as:A:ROm:c:S:F:g:vbL:M:q:n:E:QdNoBG:kCK --long full,transcripts,partial,mate_inner_distance:,min_intron_length:,max_intron_length:,procs:,librarytype:,indexpath:,refseq:,seonly,adapter_seq:,preprocess,preprocess_only,min_qual:,min_length:,percent_high_quality:,solexa,dev,oldid,solexa_p13,nonewtranscripts,leave_temp,no_hisat,cufflinks_compatible -- "$@"`
            ;;

    Darwin)
            TEMP=`getopt et:pafhr:i:I:P:l:as:A:ROm:c:S:F:g:vbL:M:q:n:E:QdNoBG:kCK $*`
            ;;

        *)
            TEMP=`getopt -o et:pafhr:i:I:P:l:as:A:ROm:c:S:F:g:q:n:E:QdNoBG:kCK --long full,transcripts,partial,mate_inner_distance:,min_intron_length:,max_intron_length:,procs:,librarytype:,indexpath:,refseq:,seonly,adapter_seq:,preprocess,preprocess_only,min_qual:,min_length:,percent_high_quality:,solexa,dev,oldid,solexa_p13,nonewtranscripts,leave_temp,no_hisat,cufflinks_compatible -- "$@"`
            ;;
esac

if [ $? != 0 ] ; then echo "Terminating..." ; exit 1 ; fi

function help_messg () {

    echo "

        -f|--full) run_type='full' ; shift ;;
        -p|--partial) run_type='partial' ; shift ;;
        -a|--transcripts) run_type='transcripts' ; shift ;;
        -k|--nonewtranscripts) no_new_txpts=1 ; shift ;;
        -r|--mate_inner_distance) mate_inner_distance_r=$2 ; shift 2 ;;
        -i|--min_intron_length) min_intron_length_i=$2 ; shift 2 ;;
        -I|--max_intron_length) max_intron_length_I=$2 ; shift 2 ;;
        -t|--procs) procs=$2 ; shift 2 ;;
        -l|--librarytype) librarytype=$2 ; shift 2 ;;
        -e|--seonly) seonly=1 ; shift ;;
        -A|--adapter_seq) adapter_seq=$2 ; shift 2 ;;
        -P|--indexpath) HISAT_INDEXES=$2 ; shift 2 ;;
        -s|--refseq) fasta_file=$2 ; shift 2 ;;
        -R|--preprocess) preprocess=1 ; shift ;;
        -O|--preprocess_only) preprocess_only=1; shift ;;
        -q|--min_qual) min_qual=$2 ; shift 2 ;;
        -n|--min_length) min_length=$2 ; shift 2 ;;
        -E|--percent_high_quality) percent_high_quality=$2 ; shift 2 ;;
        -Q|--solexa) qualscores=1 ; shift ;;
        -B|--solexa_p13) qualscores=2 ; shift ;;
        -h) help_messg ; exit ;;
        -d|--dev) dev=1 ; shift ;;
        -C|--leave_temp) leave_temp=1 ; shift ;;
        -N|--oldid) oldid=1 ; shift ;;
        --no_hisat) run_hisat=0 ; shift ;;
        -K|--cufflinks_compatible) generate output files compatible with cufflinks

        "
}

# Note the quotes around `$TEMP': they are essential!

eval set -- "$TEMP"
while true ; do
    case "$1" in
        -f|--full) run_type='full' ; shift ;;
        -p|--partial) run_type='partial' ; shift ;;
        -a|--transcripts) run_type='transcripts' ; shift ;;
        -k|--nonewtranscripts) no_new_txpts=1 ; shift ;;
        -r|--mate_inner_distance) mate_inner_distance_r=$2 ; shift 2 ;;
        -i|--min_intron_length) min_intron_length_i=$2 ; shift 2 ;;
        -I|--max_intron_length) max_intron_length_I=$2 ; shift 2 ;;
        -t|--procs) procs=$2 ; shift 2 ;;
        -l|--librarytype) librarytype=$2 ; shift 2 ;;
        -e|--seonly) seonly=1 ; shift ;;
        -A|--adapter_seq) adapter_seq=$2 ; shift 2 ;;
        -P|--indexpath) HISAT_INDEXES=$2 ; shift 2 ;;
        -s|--refseq) fasta_file=$2 ; shift 2 ;;
        -R|--preprocess) preprocess=1 ; shift ;;
        -O|--preprocess_only) preprocess_only=1; shift ;;
        -q|--min_qual) min_qual=$2 ; shift 2 ;;
        -n|--min_length) min_length=$2 ; shift 2 ;;
        -E|--percent_high_quality) percent_high_quality=$2 ; shift 2 ;;
        -Q|--solexa) qualscores=1 ; shift ;;
        -B|--solexa_p13) qualscores=2 ; shift ;;
        -h) help_messg ; exit ;;
        -d|--dev) dev=1 ; shift ;;
        -C|--leave_temp) leave_temp=1 ; shift ;;
        -N|--oldid) oldid=1 ; shift ;;
        --no_hisat) run_hisat=0 ; shift ;;
        -K|--cufflinks_compatible) cufflinks_compatible=1 ; shift ;;
        --) shift ; break ;;
        *) break ;;
    esac
done
echo "extra arguments: $@"
echo "run type is " $run_type

preprocess_flags="-i $HISAT_INDEXES -t $procs -Q $min_qual -L $min_length -H $percent_high_quality"

if [[ $seonly -eq 1 ]]
then
    preprocess_flags="$preprocess_flags -e"
fi

if [[ $leave_temp -ne 0 ]]
then
    preprocess_flags="$preprocess_flags -C"
fi


if [[ $preprocess -ne 0 ]]
then
    echo "passing reads through preprocessing routines"
    
    if [[ $seonly -eq 0 ]] # then these are paired-end data
    then
        echo "working with paired-end seuqence data"
        if [[ $adapter_seq != 'NULL' ]] # then we want to remove adapter sequence
        then
            echo "removing adapter sequence '$adapter_seq'"
            if [[ $use_cutadapt == 'NULL' ]]
            then
                echo "use adapter_trim.pl"
                echo "adapter_trim.pl --prefix --fastq --infile set1.fq --outfile set1_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters"
                adapter_trim.pl --prefix --fastq --infile set1.fq --outfile set1_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters
                echo "adapter_trim.pl --prefix --fastq --infile set2.fq --outfile set2_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters"
                adapter_trim.pl --prefix --fastq --infile set2.fq --outfile set2_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters
                echo "linking to noadapter files"
                #ln -sf set1_noadapters.fq set1.fq
                #ln -sf set2_noadapters.fq set2.fq
                touch adapter_trim_finished
            else
                echo "use cutadapt"
                echo "$cutadapt --anywhere=$adapter_seq --output=set1_noadapters.fq set1.fq 2> cutadapt_set1.log"
                $($cutadapt --anywhere=$adapter_seq --output=set1_noadapters.fq set1.fq 2> cutadapt_set1.err 1> cutadapt_set1.log)
                echo "$cutadapt --anywhere=$adapter_seq --output=set2_noadapters.fq set2.fq 2> cutadapt_set2.log"
                $($cutadapt --anywhere=$adapter_seq --output=set2_noadapters.fq set2.fq 2> cutadapt_set2.err 1> cutadapt_set1.log)
                touch cutadapt_finished
            fi
            ln -sf set1_noadapters.fq set1.fq
            ln -sf set2_noadapters.fq set2.fq
        fi
#
        echo "preprocess_fq.sh $preprocess_flags"
#        preprocess_flags="-i $HISAT_INDEXES -t $procs -Q $min_qual -L $min_length -H $percent_high_quality"
#
#        if [[ $qualscores -eq 1 ]]
#        then
#            preprocess_flags="$preprocess_flags -s"
#        fi

        preprocess_fq.sh $preprocess_flags
#
        if [[ $oldid -ne 1 ]]
        then
            echo "fastq_pe_matchup.pl --read_1 set1.fq --read_2 set2.fq --nomaxN --newid"
            fastq_pe_matchup.pl --read_1 set1.fq --read_2 set2.fq --nomaxN --newid
        else
            echo "fastq_pe_matchup.pl --read_1 set1.fq --read_2 set2.fq --nomaxN "
            fastq_pe_matchup.pl --read_1 set1.fq --read_2 set2.fq --nomaxN 
        fi
#
        echo "linking new files"
        ln -fs set1.fq.matched.fq read_1
        ln -fs set2.fq.matched.fq read_2
        ln -fs set1.fq.nomate.fq read_1.1
        ln -fs set2.fq.nomate.fq read_2.1
    else # then these are single-end data
        echo "working with single-end seuqence data"
        if [[ $adapter_seq != 'NULL' ]] # then we want to remove adapter sequence
        then
            echo "removing adapter sequence '$adapter_seq'"
            if [[ $use_cutadapt == 'NULL' ]]
            then
                echo "use adapter_trim.pl"
                echo "adapter_trim.pl --prefix --fastq --infile set1.fq --outfile set1_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters"
                adapter_trim.pl --prefix --fastq --infile set1.fq --outfile set1_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters
                touch adapter_trim_finished
            else
                echo "use cutadapt"
                echo "$cutadapt --anywhere=$adapter_seq --output=set1_noadapters.fq set1.fq 2> cutadapt_set1.log"
                $($cutadapt --anywhere=$adapter_seq --output=set1_noadapters.fq set1.fq 2> cutadapt_set1.err 1> cutadapt_set1.log)
                touch cutadapt_finished
            fi
            ln -sf set1_noadapters.fq set1.fq
        fi 
        echo "preprocess_fq.sh $preprocess_flags"
        preprocess_fq.sh $preprocess_flags
        touch preprocess_fq finished
        ln -sf set1.fq read_1.1
    fi
fi # end of preprocessing




for dir
do
    echo $dir
done


