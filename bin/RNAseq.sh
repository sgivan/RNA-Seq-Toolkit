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
echo "#######################"
echo "##                   ##"
echo "## Running RNAseq.sh ##"
echo "##                   ##"
echo "#######################"

# set path to find RNAseq scripts
wd=`pwd`
#export PATH=# <-- make sure the RST scripts are in your path

pth=$(echo $0 | sed -r 's/(.+)\/.+/\1/')

export PATH=".:$wd:$wd/bin:$HOME/bin:$PATH:${pth}"
#
# some variables to set ...
#

use_cutadapt=1
samfile="accepted_hits.sam"

#osname=`uname -s`
osname='Linux'
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
#
# INDEXES contains the path to the directory containing the reference index files.
# Be sure to end with a "/".
#
INDEXES='index'

#
# some other variables that are typically required by tophat ...
# initialize with reasonable values
#

seonly=0 # for single-end sequence data only -- no paired-end 
procs=8 # number of processors to use
run_type='NULL'
mate_inner_distance_r=0
min_intron_length_i=21
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
run_STAR=1
STARcmd="STAR"
nofilter=0
notrim=0
queue='shortQ'
gzip=0
#
# command line option parsing adpated from /usr/share/doc/util-linux-2.13/getopt-parse.bash
#
case "$osname" in

    Linux)
            TEMP=`getopt -o et:pfhr:i:I:P:aA:ROm:c:S:g:vbL:M:q:n:E:QdNoBCKFGw:z --long full,partial,mate_inner_distance:,min_intron_length:,max_intron_length:,procs:,indexpath:,refseq:,seonly,adapter_seq:,preprocess,preprocess_only,min_qual:,min_length:,percent_high_quality:,solexa,dev,oldid,solexa_p13,leave_temp,nofilter,notrim,queue:,gzip -- "$@"`
            ;;

    Darwin)
            TEMP=`getopt et:pfhr:i:I:P:aA:ROm:c:S:g:vbL:M:q:n:E:QdNoBCKFGw:z $*`
            ;;

        *)
            TEMP=`getopt -o et:pfhr:i:I:P:aA:ROm:c:S:g:q:n:E:QdNoBCKFGw:z --long full,partial,mate_inner_distance:,min_intron_length:,max_intron_length:,procs:,indexpath:,refseq:,seonly,adapter_seq:,preprocess,preprocess_only,min_qual:,min_length:,percent_high_quality:,solexa,dev,oldid,solexa_p13,leave_temp,nofilter,notrim,queue:,gzip -- "$@"`
            ;;
esac

if [ $? != 0 ] ; then echo "Terminating..." ; exit 1 ; fi

function help_messg () {

    echo "

        -f|--full) run_type='full' ; shift ;;
        -p|--partial) run_type='partial' ; shift ;;
        -r|--mate_inner_distance) mate_inner_distance_r=$2 ; shift 2 ;;
        -i|--min_intron_length) min_intron_length_i=$2 ; shift 2 ;;
        -I|--max_intron_length) max_intron_length_I=$2 ; shift 2 ;;
        -t|--procs) procs=$2 ; shift 2 ;;
        -e|--seonly) seonly=1 ; shift ;;
        -A|--adapter_seq) adapter_seq=$2 ; shift 2 ;;
        -P|--indexpath) INDEXES=$2 ; shift 2 ;;
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
        -F|--nofilter) nofilter=1 ; shift ;;
        -G|--notrim) notrim=1 ; shift ;;
        -w|--queue) queue=$2 ; shift 2 ;;
        -z|--gzip) gzip=1 ; shift ;;

        "
}

# Note the quotes around `$TEMP': they are essential!

eval set -- "$TEMP"
while true ; do
    case "$1" in
        -f|--full) run_type='full' ; shift ;;
        -p|--partial) run_type='partial' ; shift ;;
        -r|--mate_inner_distance) mate_inner_distance_r=$2 ; shift 2 ;;
        -i|--min_intron_length) min_intron_length_i=$2 ; shift 2 ;;
        -I|--max_intron_length) max_intron_length_I=$2 ; shift 2 ;;
        -t|--procs) procs=$2 ; shift 2 ;;
        -e|--seonly) seonly=1 ; shift ;;
        -A|--adapter_seq) adapter_seq=$2 ; shift 2 ;;
        -P|--indexpath) INDEXES=$2 ; shift 2 ;;
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
        -F|--nofilter) nofilter=1 ; shift ;;
        -G|--notrim) notrim=1 ; shift ;;
        -w|--queue) queue=$2 ; shift 2 ;;
        -z|--gzip) gzip=1 ; shift ;;
        --) shift ; break ;;
        *) break ;;
    esac
done
echo "extra arguments: $@"
echo "run type is " $run_type
#exit

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

preprocess_flags="-i $INDEXES -t $procs -Q $min_qual -L $min_length -H $percent_high_quality"

if [[ $gzip -eq 1 ]]
then
    preprocess_flags="$preprocess_flags -z"
fi

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

if [[ $nofilter -ne 0 ]]
then
    preprocess_flags="$preprocess_flags -n"
fi

if [[ $notrim -ne 0 ]]
then
    preprocess_flags="$preprocess_flags -N"
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
                cutadapt=$(which cutadapt)
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

# hisatmcmd, pe_extra_cmd and singles_extra_cmd contain the actual invocation of
# hisat, so they contain arguments and options passed to the program
#
STARcmd="$STARcmd --runThreadN $procs --alignIntronMin $min_intron_length_i --alignIntronMax $max_intron_length_I"

pe_extra_cmd="--alignMatesGapMax $mate_inner_distance_r --genomeDir $INDEXES -1 read_1 -2 read_2 "
singles_extra_cmd="--genomeDir $INDEXES -U read_1.1,read_2.1 "

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
echo "getting ready to run STAR"
if [[ $run_STAR -eq 1 ]]
then

#    $(mkdir -p align)
#    $(cd align)


    export INDEXES=$INDEXES # this is the directory containing the index files created with hisat-build
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

#           STAR options to include:
#           STARcmd="$STARcmd --runThreadN $procs --alignIntronMin $min_intron_length_i --alignIntronMax $max_intron_length_I"
#           pe_extra_cmd="--alignMatesGapMax $mate_inner_distance_r --genomeDir $INDEXES -1 read_1 -2 read_2 "
#           singles_extra_cmd="--genomeDir $INDEXES -U read_1.1,read_2.1 "

            echo "run type will never be transcripts. run type is " $run_type
            pwd=`pwd`
            echo "working directory: "$pwd
            if [[ $seonly -eq 0 ]]
            then
#                echo "creating STAR sbatch file with these options: create_STAR_cmd_sbatch_file.py --index $INDEXES --threads $procs --template STARcmdSE.t \
#                --minIntronLength $min_intron_length_i --maxIntronLength $max_intron_length_I"
#                $(create_STAR_cmd_sbatch_file.py --index $INDEXES --threads $procs --template STARcmdSE.t --minIntronLength $min_intron_length_i \
#                    --maxIntronLength $max_intron_length_I > cmd)
                starcmd="create_STAR_cmd_batch_file.py --index $INDEXES --threads $procs --template STARcmd.t --minIntronLength $min_intron_length_i \
                    --maxIntronLength $max_intron_length_I --alignGapMax $mate_inner_distance_r --queue $queue"

                if [[ $gzip -eq 1 ]]
                then
                    starcmd="$starcmd --gzip"
                fi

#                echo "creating STAR batch file with these options: create_STAR_cmd_batch_file.py --index $INDEXES --threads $procs --template STARcmd.t \
#                --minIntronLength $min_intron_length_i --maxIntronLength $max_intron_length_I --alignGapMax $mate_inner_distance_r --queue $queue"
#                $(create_STAR_cmd_batch_file.py --index $INDEXES --threads $procs --template STARcmd.t --minIntronLength $min_intron_length_i \
#                    --maxIntronLength $max_intron_length_I --alignGapMax $mate_inner_distance_r --queue $queue > cmd)
                echo $starcmd
                $(${starcmd} > cmd)
            else
                echo "creating STAR batch file with these options: create_STAR_cmd_batch_file.py --index $INDEXES --threads $procs --template STARcmdSE.t \
                --minIntronLength $min_intron_length_i --maxIntronLength $max_intron_length_I --queue $queue"
                $(create_STAR_cmd_batch_file.py --index $INDEXES --threads $procs --template STARcmdSE.t --minIntronLength $min_intron_length_i \
                    --maxIntronLength $max_intron_length_I --queue $queue > cmd)
            fi
        fi
    fi
                
fi

function help_messg () {

echo "

        -f|--full) run_type='full' ; shift ;;
        -p|--partial) run_type='partial' ; shift ;;
        -r|--mate_inner_distance) mate_inner_distance_r=$2 ; shift 2 ;;
        -i|--min_intron_length) min_intron_length_i=$2 ; shift 2 ;;
        -I|--max_intron_length) max_intron_length_I=$2 ; shift 2 ;;
        -t|--procs) procs=$2 ; shift 2 ;;
        -e|--seonly) seonly=1 ; shift ;;
        -A|--adapter_seq) adapter_seq=$2 ; shift 2 ;;
        -P|--indexpath) INDEXES=$2 ; shift 2 ;;
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

        "
}
