#!/bin/bash
#
# This bash script runs the program "tophat" on a series of fastq
# files of short read sequences
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
#
# set path to find RNAseq scripts
wd=`pwd`
#export PATH=# <-- make sure the RST scripts are in your path
export PATH=".:$wd:$wd/bin:$HOME/bin:$PATH"
#
# some variables to set ...
#

#tophat= # full path to tophat
tophat=tophat # because it's in my $PATH
#cufflinks= # full path to cufflinks
cufflinks=cufflinks 

osname=`uname -s`
echo "osname '$osname'"

case "$osname" in

    Linux)
            bioclass=`pwd | sed -r 's/.+\/(.+)\/(.+)/\1/'` ;
            lane=`pwd | sed -r 's/.+\/(.+)\/(.+)/\2/' | sed 's/_//'` ;;

    Darwin)
            bioclass=`pwd | sed -E 's/.+\/(.+)\/(.+)/\1/'` ;
            lane=`pwd | sed -E 's/.+\/(.+)\/(.+)/\2/' | sed 's/_//'` ;;

    *)
            bioclass=`pwd | sed -r 's/.+\/(.+)\/(.+)/\1/'` ;
            lane=`pwd | sed -r 's/.+\/(.+)\/(.+)/\2/' | sed 's/_//'` ;;

esac
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

fasta_file='refseq' # actual file should be have .fa suffix

#
# BOWTIE_INDEXES contains the path to the directory containing the reference index files
# built by the bowtie-build program. Be sure to end with a "/".
#
BOWTIE_INDEXES='index'

#
# some other variables that are typically required by tophat ...
# initialize with reasonable values
#

seonly=0 # for single-end sequence data only -- no paired-end 
procs=8 # number of processors to use
run_type='NULL'
mate_inner_distance_r=165
min_intron_length_i=50
max_intron_length_I=25000
adapter_seq='NULL'
preprocess=0
preprocess_only=0
splice_mismatches_m=0
min_anchor_length_a=8
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
#qualscores='--solexa1.3-quals'
qualscores='NULL'
dev=0
initial_read_mismatches=2
newid=0

#
# command line option parsing adpated from /usr/share/doc/util-linux-2.13/getopt-parse.bash
#
case "$osname" in

    Linux)
            TEMP=`getopt -o et:pafhr:i:I:P:l:as:A:ROm:c:S:F:g:vbL:M:q:n:E:QdC:N --long full,transcripts,partial,mate_inner_distance:,min_intron_length:,max_intron_length:,procs:,librarytype:,indexpath:,refseq:,seonly,adapter_seq:,preprocess,preprocess_only,splice_mismatches:,min_anchor_length:,mate_std_dev:,min_isoform_fraction:,max_multihits:,coverage_search,butterfly_search,segment_length:,segment_mismatches:,min_qual:,min_length:,percent_high_quality:,solexa,dev,initial_read_mismatches:,newid -- "$@"`
            ;;

    Darwin)
            TEMP=`getopt et:pafhr:i:I:P:l:as:A:ROm:c:S:F:g:vbL:M:q:n:E:QdC:N $*`
            ;;

        *)
            TEMP=`getopt -o et:pafhr:i:I:P:l:as:A:ROm:c:S:F:g:q:n:E:QdC:N --long full,transcripts,partial,mate_inner_distance:,min_intron_length:,max_intron_length:,procs:,librarytype:,indexpath:,refseq:,seonly,adapter_seq:,preprocess,preprocess_only,splice_mismatches:,min_anchor_length:,mate_std_dev:,min_isoform_fraction:,max_multihits:,coverage_search,butterfly_search,segment_length:,segment_mismatches:,min_qual:,min_length:,percent_high_quality:,solexa,dev,initial_read_mismatches:,newid -- "$@"`
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
        -r|--mate_inner_distance) mate_inner_distance_r=$2 ; shift 2 ;;
        -i|--min_intron_length) min_intron_length_i=$2 ; shift 2 ;;
        -I|--max_intron_length) max_intron_length_I=$2 ; shift 2 ;;
        -m|--splice_mismatches) splice_mismatches_m=$2 ; shift 2 ;;
        -c|--min_anchor_length) min_anchor_length_a=$2 ; shift 2 ;;
        -S|--mate_std_dev) mate_std_dev=$2 ; shift 2 ;;
        -F|--min_isoform_fraction) min_isoform_frac=$2 ; shift 2 ;;
        -g|--max_multihits) max_multihits=$2 ; shift 2 ;;
        -v|--coverage_search) coverage_search=1 ; shift ;;
        -b|--butterfly_search) butterfly_search=1 ; shift ;;
        -L|--segment_length) segment_length=$2 ; shift 2 ;;
        -M|--segment_mismatches) segment_mismatches=$2 ; shift 2 ;;
        -t|--procs) procs=$2 ; shift 2 ;;
        -l|--librarytype) librarytype=$2 ; shift 2 ;;
        -e|--seonly) seonly=1 ; shift ;;
        -A|--adapter_seq) adapter_seq=$2 ; shift 2 ;;
        -P|--indexpath) BOWTIE_INDEXES=$2 ; shift 2 ;;
        -s|--refseq) fasta_file=$2 ; shift 2 ;;
        -R|--preprocess) preprocess=1 ; shift ;;
        -O|--preprocess_only) preprocess_only=1; shift ;;
        -q|--min_qual) min_qual=$2 ; shift 2 ;;
        -n|--min_length) min_length=$2 ; shift 2 ;;
        -E|--percent_high_quality) percent_high_quality=$2 ; shift 2 ;;
        -Q|--solexa) qualscores=1 ; shift ;;
        -h) help_messg ; exit ;;
        -d|--dev) dev=1 ; shift ;;
        -C|--initial_read_mismatches) initial_read_mismatches=$2 ; shift 2 ;;
        -N|--newid) newid=1 ; shift ;;
        --) shift ; break ;;
        *) break ;;
    esac
done
echo "extra arguments: $@"
echo "run type is " $run_type
#exit

#
# tophatcmd, pe_extra_cmd and singles_extra_cmd contain the actual invocation of
# tophat, so it contains arguments and options passed to the program
#

#tophatcmd="$tophat --library-type $librarytype -p $procs -i $min_intron_length_i -I $max_intron_length_I --solexa1.3-quals -m $splice_mismatches_m -a $min_anchor_length_a --min-isoform-fraction $min_isoform_frac --max-multihits $max_multihits --segment-length $segment_length --segment-mismatches $segment_mismatches"
#tophatcmd="$tophat --library-type $librarytype -p $procs -i $min_intron_length_i -I $max_intron_length_I -m $splice_mismatches_m -a $min_anchor_length_a --min-isoform-fraction $min_isoform_frac --max-multihits $max_multihits --segment-length $segment_length --segment-mismatches $segment_mismatches"
tophatcmd="$tophat --library-type $librarytype -p $procs -i $min_intron_length_i -I $max_intron_length_I -m $splice_mismatches_m -a $min_anchor_length_a --min-isoform-fraction $min_isoform_frac --max-multihits $max_multihits --segment-length $segment_length --segment-mismatches $segment_mismatches --initial-read-mismatches $initial_read_mismatches"
pe_extra_cmd="-r $mate_inner_distance_r --mate-std-dev $mate_std_dev -o pe_tophat_out $BOWTIE_INDEXES/$fasta_file read_1 read_2 "
singles_extra_cmd="-o singles_tophat_out $BOWTIE_INDEXES/$fasta_file read_1.1,read_2.1 "
# 
#cufflinkscmd="$cufflinks -m $mate_inner_distance_r -I $max_intron_length_I "
#cufflinksflgs="-I $max_intron_length_I --library-type $librarytype -r ../index/$fasta_file.fa -p $procs -o cufflinks -L $bioclass$lane --min-intron-length $min_intron_length_i */accepted_hits.bam"
#cufflinksflgs="-I $max_intron_length_I --library-type $librarytype -r ../index/$fasta_file.fa -p $procs -o cufflinks -L $bioclass$lane --min-intron-length $min_intron_length_i"
#cufflinksflgs="-I $max_intron_length_I --library-type $librarytype -r $wd/index/$fasta_file.fa -p $procs -o cufflinks -L $bioclass$lane --min-intron-length $min_intron_length_i"
#cufflinksflgs="-I $max_intron_length_I --library-type $librarytype -r $BOWTIE_INDEXES/$fasta_file.fa -p $procs -o cufflinks -L $bioclass$lane --min-intron-length $min_intron_length_i"
# cufflinks version >= 1.0.0 no longer uses -r flag. Now uses -b flag
cufflinksflgs="-I $max_intron_length_I --library-type $librarytype -b $BOWTIE_INDEXES/$fasta_file.fa -p $procs -o cufflinks -L $bioclass$lane --min-intron-length $min_intron_length_i"

#
# END OF USER-DEFINED VARIABLES
# 
# Typically, it shouldn't be necessay to change anything below this line
#
#
if [[ $dev -ne 0 ]]
then
    export $PATH="/ircf/ircfapps/dev/bin:$PATH"
fi

if [[ $run_type = "full" ]]
then
    echo "setting preprocess to TRUE"
    preprocess=1
fi

if [[ $coverage_search -ne 0 ]]
then
    tophatcmd="$tophatcmd --coverage-search"
fi
if [[ $butterfly_search -ne 0 ]]
then
    tophatcmd="$tophatcmd --butterfly-search"
fi

if [[ $qualscores != 'NULL' ]]
then
    if [[ $qualscores -eq 1 ]]
    then
        tophatcmd="$tophatcmd --solexa-quals"
    fi
else
    tophatcmd="$tophatcmd --solexa1.3-quals"
fi

preprocess_flags="-i $BOWTIE_INDEXES -t $procs -Q $min_qual -L $min_length -H $percent_high_quality"
if [[ $qualscores -eq 1 ]]
then
    preprocess_flags="$preprocess_flags -s"
fi
if [[ $seonly -eq 1 ]]
then
    preprocess_flags="$preprocess_flags -e"
fi

#echo ""
#echo "tophatcmd= $tophatcmd"
#echo ""
#exit

if [[ $preprocess -ne 0 ]]
then
    echo "passing reads through preprocessing routines"
    
    if [[ $seonly -eq 0 ]] # then these are paired-end data
    then
        echo "working with paired-end seuqence data"
        if [[ $adapter_seq != 'NULL' ]] # then we want to remove adapter sequence
        then
            echo "removing adapter sequence '$adapter_seq'"
#            cd preprocess
            echo "adapter_trim.pl --fastq --infile set1.fq --outfile set1_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters"
            adapter_trim.pl --fastq --infile set1.fq --outfile set1_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters
            echo "adapter_trim.pl --fastq --infile set2.fq --outfile set2_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters"
            adapter_trim.pl --fastq --infile set2.fq --outfile set2_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters
            echo "linking to noadapter files"
            ln -sf set1_noadapters.fq set1.fq
            ln -sf set2_noadapters.fq set2.fq
#            cd ..
        fi
#
        echo "preprocess_fq.sh $preprocess_flags"
#        preprocess_flags="-i $BOWTIE_INDEXES -t $procs -Q $min_qual -L $min_length -H $percent_high_quality"
#
#        if [[ $qualscores -eq 1 ]]
#        then
#            preprocess_flags="$preprocess_flags -s"
#        fi

        preprocess_fq.sh $preprocess_flags
#
        if [[ $newid = 1 ]]
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
#            cd preprocess
            echo "adapter_trim.pl --fastq --infile set1.fq --outfile set1_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters"
            adapter_trim.pl --fastq --infile set1.fq --outfile set1_noadapters.fq -adapterseq $adapter_seq --overwrite --printall --notwoadapters
            ln -sf set1_noadapters.fq set1.fq
#            cd ..
#            ln -sf preprocess/set1.fq ./
        fi 
        echo "preprocess_fq.sh $preprocess_flags"
        preprocess_fq.sh $preprocess_flags
        ln -sf set1.fq read_1.1
    fi
fi # end of preprocessing

if [[ $preprocess_only -eq 1 ]]
then
    echo "preprocess only"
    exit
fi

export BOWTIE_INDEXES=$BOWTIE_INDEXES # this is the directory containing the index files created with bowtie-build
if [ $run_type = transcripts ]
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

    if [[ $run_type != 'NULL' ]]
    then
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
fi
            
if [[ run_type != 'NULL' ]]
then

    mkdir -p merged
    cd merged

    if [[ $seonly -eq 0 ]]
    then

        echo "merging PE and SE bam files"
        samtools merge merged.bam ../*/accepted_hits.bam

    else

        ln -s ../singles_tophat_out/accepted_hits.bam ./merged.bam
    fi

    cd ..

fi

#if [[ $run_type = full ]] # not sure why this is just for full runs
if [[ $run_type = full ]] || [[ $run_type = partial ]]
then

#    mkdir -p merged
#    cd merged

    if [[ $seonly -eq 0 ]]
    then
#        echo "now merging PE and SE alignment data"
        echo "using cufflinks to build gene models with PE and SE alignment data"
        
        #echo "creating aggregate_junctions.txt"
        #cat ../*/junctions.bed | awk '{ if ($1 != "track") {split($11,len,","); split($12,blstrt,","); printf "%s\t%i\t%i\t%s\n", $1, $2 + len[1] - 1, $2 + blstrt[2], $6; }}' | sort -k 1,1 -gk 2,2 | uniq > aggregate_junctions.txt
#        echo "merging tophat bam files"
#        samtools merge merged.bam ../*/accepted_hits.bam
#        echo $cufflinks $cufflinksflgs */accepted_hits.bam
#        $cufflinks $cufflinksflgs */accepted_hits.bam
        echo $cufflinks $cufflinksflgs merged/merged.bam
        $cufflinks $cufflinksflgs merged/merged.bam
    else
        echo "using cufflinks to build gene models with SE alignment data"
#        ln -fs ../singles_tophat_out/accepted_hits.bam ./merged.bam
#        echo $cufflinks $cufflinksflgs */accepted_hits.bam
#        $cufflinks $cufflinksflgs */accepted_hits.bam
        echo $cufflinks $cufflinksflgs merged/merged.bam
        $cufflinks $cufflinksflgs merged/merged.bam
    fi
#    cd ..
fi

if [ $run_type = transcripts ]
then
    echo "running cufflinks using transcript file"
    #samtools merge - ../*/accepted_hits.bam | samtools view -o - - | cuff_sam_to_gff.pl --infile - --outfile all_sorted.gff --source GA2 --type $bioclass$lane
    
    echo "creating new directory: merged_aggregate"
#    mkdir -p merged_aggregate
#    cd merged_aggregate
#    echo "creating bam file from pe and single runs"
#    if [[ $seonly -eq 0 ]]
#    then
#        samtools merge merged.bam ../*/accepted_hits.bam 
#    else
#        ln -s ../singles_tophat_out/accepted_hits.bam ./merged.bam
#    fi
    echo "running cufflinks"
    #cufflinks_extra_cmd="--GTF ../transcripts.gtf -L $bioclass$lane merged.bam"
    #cufflinks_extra_cmd="--GTF ../transcripts.gtf"
    cufflinks_extra_cmd="--GTF transcripts.gtf"
    echo $cufflinks $cufflinksflgs $cufflinks_extra_cmd */accepted_hits.bam
    $cufflinks $cufflinksflgs $cufflinks_extra_cmd */accepted_hits.bam
    #cd ..
    #ln -s merged/aggregate_junctions.txt junctions.txt
    #ln -s merged/transcripts.gtf ./
    
    #echo "rerunning tophat with merged data"
    
fi

#mkdir -p merged
#cd merged
#
#if [[ $seonly -eq 0 ]]
#then
#
#    echo "merging PE and SE bam files"
#    samtools merge merged.bam ../*/accepted_hits.bam
#
#else
#
#    ln -s ../singles_tophat_out/accepted_hits.bam ./merged.bam
#fi
#
#cd ..

echo "finished"
echo ""

