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

echo "creating symbolic link to RNAseq tools directory"
ln -s ../../bin ./
echo "setting PATH"
wd=`pwd`
export PATH=".:$wd:$wd/bin:$PATH"
echo "making index directory"
mkdir index
ln -sf index hisat_index
echo "moving reference sequence and undesireables into index directory"
mv chrom3.fa index
mv Contaminants.fa index
echo "creating sybolic links"
cd index
ln -s chrom3.fa refseq.fa
ln -s Contaminants.fa filter.fa
cd ..
echo "creating sample directories"
mkdir -p s_5 s_6 s_7 s_8
echo "moving sample fastq files into their respective directories"
mv s_5_* s_5
mv s_6_* s_6
mv s_7_* s_7
mv s_8_* s_8
echo "creating symbolic links inside of sample directories"
cd s_5
ln -sf s_5_1_test.txt set1.fq
cd ../s_6
ln -sf s_6_1_test.txt set1.fq
cd ../s_7
ln -sf s_7_1_test.txt set1.fq
cd ../s_8
ln -sf s_8_1_test.txt set1.fq
cd ..
#
echo "making bowtie indices"
#
cd index
echo "building refseq index"
#bowtie-build chrom3.fa refseq
hisat2-build refseq.fa refseq.fa
echo "building filter index"
bowtie-build filter.fa filter.fa
cd ..
bash cmd
echo "finished"
