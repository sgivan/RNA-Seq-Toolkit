#!/bin/bash
#
# This bash script runs the program "hisat" on a series of fastq
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

hisat=$(which hisat2) # because it's in my $PATH
stringtie=$(which stringtie)
cutadapt=$(which cutadapt)
use_cutadapt=1
samfile="accepted_hits.sam"
use_stringtie=1

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
# librarytype contains the "library type", as defined in the hisat documentation:
# http://ccb.jhu.edu/software/hisat2/manual.shtml#running-hisat2 
# the default is unstranded, fr_firstrand = R, fr_secondstrand = F
#

librarytype='NULL'

#
# fasta_file contains the root of the fasta file that contains the reference sequence that
# the index files are derived from. It should be in the same location as the hisat index
# files made with hisat-build, which is specified by HISAT_INDEXES, below
#

fasta_file='refseq.fa' # actual file should be have .fa suffix

#
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
mate_inner_distance_r=165
min_intron_length_i=20
max_intron_length_I=500000
adapter_seq='NULL'
preprocess=0
preprocess_only=0
min_qual=13
min_length=32
percent_high_quality=90
#qualscores='--solexa1.3-quals'
qualscores='NULL'
dev=0
leave_temp=0
oldid=0
no_new_txpts='NULL'
run_hisat=1
cufflinks_compatible=0
#
# command line option parsing adpated from /usr/share/doc/util-linux-2.13/getopt-parse.bash
#
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
#exit

#
# hisatmcmd, pe_extra_cmd and singles_extra_cmd contain the actual invocation of
# hisat, so they contain arguments and options passed to the program
#
# retain cufflinks compatibility with the --dta-cufflinks option
#hisatcmd="$hisat --dta --no-unal -p $procs --min-intronlen $min_intron_length_i --max-intronlen $max_intron_length_I"
hisatcmd="$hisat --no-unal -p $procs --min-intronlen $min_intron_length_i --max-intronlen $max_intron_length_I"

if [[ cufflinks_compatible -ne 0 ]]
then
    hisatcmd="$hisatcmd --dta-cufflinks"
else
    hisatcmd="$hisatcmd --dta"
fi

pe_extra_cmd="--met-file hisat_metrics_pe.txt -S pe_hisat_out/$samfile --un-gz pe_hisat_out/unaligned.fa.gz --un-conc-gz --un-gz pe_hisat_out/pe_unaligned.fa.gz -x $HISAT_INDEXES/$fasta_file -1 read_1 -2 read_2 "
singles_extra_cmd="--met-file hisat_metrics_se.txt -S singles_hisat_out/$samfile --un-gz singles_hisat_out/unaligned.fa.gz -x $HISAT_INDEXES/$fasta_file -U read_1.1,read_2.1 "

if [[ cufflinks_compatible -ne 0 ]]
then
    pe_extra_cmd="$pe_extra_cmd --dta-cufflinks"
fi
# 
#cufflinksflgs="-u --max_intron_length $max_intron_length_I -b $HISAT_INDEXES/$fasta_file.fa -p $procs -o cufflinks -L $bioclass$lane --min_intron_length $min_intron_length_i"
#stringtieflgs="-o stringtie/string_transcripts.gtf -B -p $procs"
#stringtieflgs="-o stringtie/string_transcripts.gtf -B -p $procs -l $bioclass$lane"
#stringtieflgs="-o ../ballgown/$bioclass$lane/"$bioclass$lane"_transcripts.gtf -B -p $procs -l $bioclass$lane"
#stringtieflgs="-o ../ballgown/$bioclass$lane/transcripts.gtf -B -p $procs -l $bioclass$lane"
# not sure why I use the -B flag, but it is causing an error witout guide the GFF/GTF
# answer: the -B flag makes stringtie generate a file for ballgown, the DE program
# So, don't use the -B flag at the merge stage, but use it at the DE stage
#stringtieflgs="-o ../ballgown/$bioclass$lane/transcripts.gtf -p $procs -l $bioclass$lane"
# change the above to be more consistent with past directory structures.
# this will likely require changing the ballgown Rscript
stringtieflgs="-o ballgown/transcripts.gtf -p $procs -l $bioclass$lane"

#
# END OF USER-DEFINED VARIABLES
# 
# Typically, it shouldn't be necessay to change anything below this line
#
#
if [[ $dev -ne 0 ]]
then
    export PATH="/ircf/ircfapps/dev/bin:$PATH"
fi

if [[ $run_type = "full" ]]
then
    echo "setting preprocess to TRUE"
    preprocess=1
fi

# default quality scores are Phred+33
if [[ $qualscores != 'NULL' ]]
then
    if [[ $qualscores -eq 1 ]]
    then
        # quality scores are Phred+64
        #hisatcmd="$hisatcmd --solexa-quals"
        hisatcmd="$hisatcmd --solexa1.3-quals"
    elif [[ $qualscores -eq 2 ]]
    then
        hisatcmd="$hisatcmd --solexa-quals"
    fi
else
    #hisatcmd="$hisatcmd --solexa1.3-quals"
    hisatcmd=$hisatcmd
fi

if [[ $librarytype -ne 'NULL' ]]
then
    hisatcmd="$hisatcmd --rna-strandness $librarytype"
fi

preprocess_flags="-i $HISAT_INDEXES -t $procs -Q $min_qual -L $min_length -H $percent_high_quality"
if [[ $qualscores -eq 1 ]]
then
    preprocess_flags="$preprocess_flags -s"
    #preprocess_flags=$preprocess_flags
fi

if [[ $seonly -eq 1 ]]
then
    preprocess_flags="$preprocess_flags -e"
fi

if [[ $leave_temp -ne 0 ]]
then
    preprocess_flags="$preprocess_flags -C"
fi

#echo ""
#echo "hisatcmd= $hisatcmd"
#echo ""
#exit

# start running pipeline


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

if [[ $preprocess_only -eq 1 ]]
then
    echo "preprocess only"
    exit
fi

# run hisat

if [[ $run_hisat -eq 1 ]]
then

    $(mkdir -p singles_hisat_out)
    if [[ $seonly -eq 0 ]]
    then
        $(mkdir -p pe_hisat_out)
    fi


    export HISAT_INDEXES=$HISAT_INDEXES # this is the directory containing the index files created with hisat-build
    if [ $run_type = transcripts ]
    then
        echo "run type is " $run_type
        if [[ $no_new_txpts != "NULL" ]]
        then
            pe_extra_cmd=" --transcriptome-mapping-only $pe_extra_cmd"
            singles_extra_cmd=" --transcriptome-mapping-only $singles_extra_cmd"
        fi

        if [[ $seonly -eq 0 ]]
        then
            echo "running hisat with these options: "$hisatcmd $pe_extra_cmd $@
            $hisatcmd $pe_extra_cmd $@ > pe_hisat.stdout 2>&1
            echo "running hisat with these options: "$hisatcmd $singles_extra_cmd $@
            $hisatcmd $singles_extra_cmd  $@ > singles_hisat.stdout 2>&1
        else
            singles_extra_cmd=`echo $singles_extra_cmd | sed 's/,read_2.1//'`
            echo "running hisat with these options: "$hisatcmd $singles_extra_cmd $@
            pwd=`pwd`
            echo "working directory: "$pwd
            $hisatcmd $singles_extra_cmd $@ > singles_hisat.stdout 2>&1
        fi

    else # maybe this should be a separate if clause

        if [[ $run_type != 'NULL' ]]
        then 
            echo "run type is " $run_type
            echo "running tophat without precomputed annotations"
            pwd=`pwd`
            echo "working directory: "$pwd
            if [[ $seonly -eq 0 ]]
            then
                echo "running HISAT2 with these options: "$hisatcmd $pe_extra_cmd $@
                $hisatcmd $pe_extra_cmd $@ > pe_hisat.stdout 2>&1
                echo "running HISAT2 with these options: "$hisatcmd $singles_extra_cmd $@
                $hisatcmd $singles_extra_cmd $@ > singles_hisat.stdout 2>&1
            else
                singles_extra_cmd=`echo $singles_extra_cmd | sed 's/,read_2.1//'`
                echo "running HISAT2 with these options: "$hisatcmd $singles_extra_cmd $@
                $hisatcmd $singles_extra_cmd $@ > singles_hisat.stdout 2>&1
            fi
        fi
    fi
                
    if [[ run_type != 'NULL' ]]
    then

        if [[ -e "merged" ]] || mkdir -p merged
        then
            cd merged

            if [[ $seonly -eq 0 ]]
            then

                echo "merging PE and SE sam files"
                # I don't think this will work bc samfiles aren't sorted
                $(samtools merge merged.bam ../*/$samfile)

            else
                ln -s ../singles_hisat_out/$samfile ./merged.sam
            fi

            cd ..
        else
            echo "can't create merged directory"
        fi

    fi
fi

# run stringtie

#if [[ $run_type = full ]] # not sure why this is just for full runs
if [[ $run_type = full ]] || [[ $run_type = partial ]] 
then

    cd merged
    # test to see if merged_sorted.bam exists first to save time
    if [ ! -e merged_sorted.bam ]
    then
        $(samtools view -Su merged.sam | samtools sort -o merged_sorted.bam -)
        ln -s merged_sorted.bam merged.bam
    fi
    cd ..

    if [[ $seonly -eq 0 ]]
    then
#        echo "now merging PE and SE alignment data"
        echo "using stringtie to build gene models with PE and SE alignment data"
        
        echo $stringtie $stringtieflgs merged/merged.bam
        $($stringtie $stringtieflgs merged/merged.bam > stringtie.log 2>&1)
    else
        echo "using stringtie to build gene models with SE alignment data"
        echo $stringtie $stringtieflgs merged/merged.bam
        $($stringtie $stringtieflgs merged/merged.bam > stringtie.log 2>&1)
    fi
#    cd ..
fi

if [ $run_type = transcripts ]
then
    cd merged
    echo 'creating sorted bam file for input to stringtie'
    if [ ! -e merged_sorted.bam ]
    then
        $(samtools view -Su merged.sam | samtools sort -o merged_sorted.bam -)
        ln -s merged_sorted.bam merged.bam
    fi
    cd ..

    echo "running stringtie using transcript file"
    
    echo "running stringtie"
    if [[ $no_new_txpts != "NULL" ]]
    then
        stringtie_extra_cmd=" -B -G transcripts.gtf -e"
    else
        stringtie_extra_cmd=" -B -G transcripts.gtf"
    fi
    echo $stringtie $stringtieflgs $stringtie_extra_cmd merged/merged.bam 
    $($stringtie $stringtieflgs $stringtie_extra_cmd merged/merged.bam > stringtie.log 2>&1)
    
fi

echo "finished"
echo ""

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

        "
}
