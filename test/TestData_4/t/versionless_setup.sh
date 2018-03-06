#!/bin/bash
#
#   Copyright Scott A. Givan, University of Missouri, July 6, 2012.
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
module load Python-shared
module load R-3.3.0-sharedlib
module load bowtie2-2.3.2
module load stringtie-1.3.0
module load HISAT2-2.1.0

# Sanity check (Are we where we want to be and have the files been copied?):
sanity_check() {
    echo "------------- Files in `pwd` and subfolders -------------"
    ls -A *
    echo "---------------------------------------------------------"
}

sanity_check

echo "running using hisat2 version '`hisat2 --version | head -n1`' (specifically: `which hisat2`)"

for file in *.gz; do echo "gunzip $file"; gunzip $file; done

sanity_check

echo "creating symbolic link to RNAseq tools directory"
ln -sf ../../bin ./
echo "setting PATH"
wd=`pwd`
export PATH=".:$wd:$wd/bin:$PATH"
echo "making index directory"
mkdir index
ln -sf index hisat_index
echo "moving reference sequence and undesireables into index directory"
cp Chr19.fa index
cp Contaminants.fa index
cp Chr19.gtf transcripts.gtf
echo "creating symbolic links"
cd index
ln -s Chr19.fa refseq.fa
ln -s Contaminants.fa filter.fa

hisat2_extract_exons.py ../transcripts.gtf > exons.txt
hisat2_extract_splice_sites.py ../transcripts.gtf > splice_sites.txt
hisat2-build --threads 4 --exon exons.txt --ss splice_sites.txt refseq.fa refseq.fa

sanity_check

cd ..

sanity_check

#echo "creating sample directories"
mkdir -p s_1 s_2 s_3 s_4
echo "moving sample fastq files into their respective directories"
cp s_1_hits.fastq s_1
cp s_2_hits.fastq s_2
cp s_3_hits.fastq s_3
cp s_4_hits.fastq s_4

sanity_check

echo "creating symbolic links inside of sample directories"
cd s_1
ln -sf s_1_hits.fastq set1.fq
cd ../s_2
ln -sf s_2_hits.fastq set1.fq
cd ../s_3
ln -sf s_3_hits.fastq set1.fq
cd ../s_4
ln -sf s_4_hits.fastq set1.fq
cd ..
#
#echo "making bowtie indices"
#
cd index
#echo "building refseq index"
##bowtie-build chrom3.fa refseq
#hisat2-build --threads 4 refseq.fa refseq.fa
echo "building filter index"
bowtie-build filter.fa filter.fa

sanity_check

cd ..

sanity_check
