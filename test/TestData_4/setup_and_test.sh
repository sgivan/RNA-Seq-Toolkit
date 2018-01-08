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
module load bowtie2-2.3.2
module load stringtie-1.3.0
module load HISAT2-2.1.0
module load Python-shared
module load R-3.3.0-sharedlib

echo "creating symbolic link to RNAseq tools directory"
ln -sf ../../bin ./
echo "setting PATH"
wd=`pwd`
export PATH=".:$wd:$wd/bin:$PATH"
#echo "making index directory"
mkdir index
ln -sf index hisat_index
#echo "moving reference sequence and undesireables into index directory"
mv Chr19.fa index
mv Contaminants.fa index
echo "creating sybolic links"
cd index
ln -s chr19.fa refseq.fa
ln -s Contaminants.fa filter.fa
#cd ..
#echo "creating sample directories"
mkdir -p s_1 s_2 s_3 s_4
echo "moving sample fastq files into their respective directories"
mv s_1.fq s_1
mv s_2.fq s_2
mv s_3.fq s_3
mv s_4.fq s_4
echo "creating symbolic links inside of sample directories"
cd s_1
ln -sf s_1.fq set1.fq
cd ../s_2
ln -sf s_2.fq set1.fq
cd ../s_3
ln -sf s_3.fq set1.fq
cd ../s_4
ln -sf s_4.fq set1.fq
cd ..
#
#echo "making bowtie indices"
#
cd index
#echo "building refseq index"
##bowtie-build chrom3.fa refseq
#hisat2-build refseq.fa refseq.fa
echo "building filter index"
bowtie-build filter.fa filter.fa
cd ..
cp Chr19.gtf transcripts.gtf
echo "running RNAseq_process_data.sh"
bash cmd
echo "finished"
