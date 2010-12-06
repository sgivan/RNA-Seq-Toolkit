#!/bin/bash
#
# This bash script runs the program "tophat" on a series of fastq
# files of short read sequences
#
# copyright Scott Givan, The University of Missouri, November 11, 2011.
#
# set path to find RNAseq scripts

export PATH="~givans/projects/RNAseq/bin:~givans/bin:$PATH"

#
# some variables to set ...
#

#tophat=/home/cgrb/givans/bin/tophat # full path to tophat
tophat=tophat # because it's in my $PATH
#cufflinks=/home/cgrb/givans/bin/cufflinks # full path to cufflinks
cufflinks=cufflinks 

#
# librarytype contains the "library type", as defined in the tophat documentation:
# http://tophat.cbcb.umd.edu/manual.html
# For Illumina, it is usually fr-unstranded
#

librarytype=fr-unstranded 

#
# fasta_file contains the root of the fasta file that contains the reference sequence that
# the index files are derived from. It should be in the same location as the bowtie index
# files made with bowtie-build, which is specified by BOWTIE_INDEXES, below
#

fasta_file=maize # actual file should be have .fa suffix

#
# BOWTIE_INDEXES contains the path to the directory containing the reference index files
# built by the bowtie-build program. Be sure to end with a "/".
#
#BOWTIE_INDEXES=/dbase/genomes/maize/4a.53v2/repeat_masked/
#BOWTIE_INDEXES=/dbase/genomes/maize/bowtie/

#BOWTIE_INDEXES=/dbase/genomes/maize/4a.53v2/chrv2/bowtie/
#BOWTIE_INDEXES=/dbase/genomes/maize/4a.53v2/chrv2/bowtie_shortchroms/
BOWTIE_INDEXES='index'

#
# some other variables that are typically required by tophat ...
# initialize with reasonable values
#

seonly=0 # for single-end sequence data only -- no paired-end 
procs=8 # number of processors to use
run_type='partial'
mate_inner_distance_r=165
min_intron_length_i=50
max_intron_length_I=25000
adapter_seq='NULL'

#
# command line option parsing adpated from /usr/share/doc/util-linux-2.13/getopt-parse.bash
#

TEMP=`getopt -o pafhr:i:I:P:l:sA: --long full,aggregate,partial,mate_inner_distance:,min_intron_length:,max_intron_length:,procs:,librarytype:,indexpath:,refseq:,seonly,adapter_seq: -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!

eval set -- "$TEMP"
while true ; do
    case "$1" in
        -f|--full) run_type='full' ; shift ;;
        -p|--partial) run_type='partial' ; shift ;;
        -a|--aggregate) run_type='aggregate' ; shift ;;
        -r|--mate_inner_distance) mate_inner_distance_r=$2 ; shift 2 ;;
        -i|--min_intron_length) min_intron_length_i=$2 ; shift 2 ;;
        -I|--max_intron_length) max_intron_length_I=$2 ; shift 2 ;;
        -P|--procs) procs=$2 ; shift 2 ;;
        -l|--librarytype) librarytype=$2 ; shift 2 ;;
        -s|--seonly) seonly=1 ; shift ;;
        -A|--adapter_seq) adapter_seq=$2 ; shift 2 ;;
        --indexpath) BOWTIE_INDEXES=$2 ; shift 2 ;;
        --refseq) fasta_file=$2 ; shift 2 ;;
        -h) help_messg ; exit ;;
        --) shift ; break ;;
        *) break ;;
    esac
done
#echo "extra arguments: $@"

#
# tophatcmd, pe_extra_cmd and singles_extra_cmd contain the actual invocation of
# tophat, so it contains arguments and options passed to the program
#

#tophatcmd="$tophat --library-type $librarytype -p $procs -r $mate_inner_distance_r -i $min_intron_length_i -I $max_intron_length_I --solexa1.3-quals"
tophatcmd="$tophat --library-type $librarytype -p $procs -i $min_intron_length_i -I $max_intron_length_I --solexa1.3-quals"
pe_extra_cmd="-r $mate_inner_distance_r -o pe_tophat_out $fasta_file read_1 read_2 "
singles_extra_cmd="-o singles_tophat_out $fasta_file read_1.1,read_2.1 "
# 
#cufflinkscmd="$cufflinks -m $mate_inner_distance_r -I $max_intron_length_I "
cufflinkscmd="$cufflinks -I $max_intron_length_I --library-type $librarytype -r ../../refseq -p $procs"

#
# END OF USER-DEFINED VARIABLES
# 
# Typically, it shouldn't be necessay to change anything below this line
#
#

if [ $run_type = full ]
then
    echo "preprocess_fq.sh"
    preprocess_fq.sh --indexpath $BOWTIE_INDEXES
    
    if [[ $seonly -eq 0 ]] # then these are paired-end data
    then
        echo "working with paired-end seuqence data"
        if [[ $adapter_seq != 'NULL' ]] # then we want to remove adapter sequence
        then
            echo "removing adapter sequence '$adapter_seq'"
            cd preprocess
            echo "adapter_trim.pl --fastq --infile set1.fq --outfile set1_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters"
            adapter_trim.pl --fastq --infile set1.fq --outfile set1_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters
            echo "adapter_trim.pl --fastq --infile set2.fq --outfile set2_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters"
            adapter_trim.pl --fastq --infile set2.fq --outfile set2_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters
            echo "linking new files"
            ln -sf set1_noadapters.fq set1.fq
            ln -sf set2_noadapters.fq set2.fq
            cd ..
        fi
        echo "fastq_pe_matchup.pl --read_1 set1.fq --read_2 set2.fq --nomaxN"
        fastq_pe_matchup.pl --read_1 set1.fq --read_2 set2.fq --nomaxN
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
            cd preprocess
            echo "adapter_trim.pl --fastq --infile set1.fq --outfile set1_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters"
            adapter_trim.pl --fastq --infile set1.fq --outfile set1_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters
            ln -sf set1_noadapters.fq set1.fq
            cd ..
            ln -sf preprocess/set1.fq ./
        fi 
        ln -sf set1.fq read_1.1
    fi
fi

export BOWTIE_INDEXES # this is the directory containing the index files created with bowtie-build
if [ $run_type = aggregate ]
then
    if [[ $seonly -eq 0 ]]
    then
        echo "running tophat using precomputed annotation file"
        echo $tophatcmd --GTF transcripts.gtf $pe_extra_cmd $@
        $tophatcmd --GTF transcripts.gtf $pe_extra_cmd $@ > pe_tophat.stdout 2>&1
        echo $tophatcmd --GTF transcripts.gtf $singles_extra_cmd $@
        $tophatcmd --GTF transcripts.gtf $singles_extra_cmd  $@ > singles_tophat.stdout 2>&1
    else
        singles_extra_cmd=`echo $singles_extra_cmd | sed 's/,read_2.1//'`
        echo $tophatcmd --GTF transcripts.gtf $singles_extra_cmd $@
        $tophatcmd --GTF transcripts.gtf $singles_extra_cmd $@ > singles_tophat.stdout 2>&1
    fi

else # maybe this should be a separate if clause
    echo "running tophat without precomputed annotations"
    if [[ $seonly -eq 0 ]]
    then
        echo $tophatcmd $pe_extra_cmd $@
        $tophatcmd $pe_extra_cmd $@ > pe_tophat.stdout 2>&1
        echo $tophatcmd $singles_extra_cmd $@
        $tophatcmd $singles_extra_cmd $@ > singles_tophat.stdout 2>&1
    else
        singles_extra_cmd=`echo $singles_extra_cmd | sed 's/,read_2.1//'`
        echo $tophatcmd $singles_extra_cmd $@
        $tophatcmd $singles_extra_cmd $@ > singles_tophat.stdout 2>&1
    fi
        
fi

if [[ $run_type -eq full ]]
then

    mkdir -p merged
    cd merged

    if [[ $seonly -eq 0 ]]
    then
        echo "now merging PE and SE alignment data"
        
        #echo "creating aggregate_junctions.txt"
        #cat ../*/junctions.bed | awk '{ if ($1 != "track") {split($11,len,","); split($12,blstrt,","); printf "%s\t%i\t%i\t%s\n", $1, $2 + len[1] - 1, $2 + blstrt[2], $6; }}' | sort -k 1,1 -gk 2,2 | uniq > aggregate_junctions.txt
        echo "merging tophat bam files"
        samtools merge merged.bam ../*/accepted_hits.bam
    else
        ln -fs ../singles_tophat_out/accepted_hits.bam ./merged.bam
    fi
    cd ..
fi

if [ $run_type = aggregate ]
then
    echo "running cufflinks using aggregate transcript file"
    #samtools merge - ../*/accepted_hits.bam | samtools view -o - - | cuff_sam_to_gff.pl --infile - --outfile all_sorted.gff --source GA2 --type $bioclass$lane
    
    bioclass=`pwd | sed -r 's/.+\/(.+)\/(.+)/\1/'`
    lane=`pwd | sed -r 's/.+\/(.+)\/(.+)/\2/' | sed 's/_//'`
    
    echo "creating new directory: merged_aggregate"
    mkdir -p merged_aggregate
    cd merged_aggregate
    echo "creating bam file from pe and single runs"
    if [[ $seonly -eq 0 ]]
    then
        samtools merge merged.bam ../*/accepted_hits.bam 
    else
        ln -s ../singles_tophat_out/accepted_hits.bam ./merged.bam
    fi
    echo "running cufflinks"
    #cufflinks_extra_cmd="--GTF ../transcripts.gtf -L $bioclass$lane merged.bam > cufflinks.stdout 2>&1"
    cufflinks_extra_cmd="--GTF ../transcripts.gtf -L $bioclass$lane merged.bam"
    echo "$cufflinks $cufflinks_extra_cmd"
    $cufflinks $cufflinks_extra_cmd
    #cd ..
    #ln -s merged/aggregate_junctions.txt junctions.txt
    #ln -s merged/transcripts.gtf ./
    
    #echo "rerunning tophat with merged data"
    
fi
