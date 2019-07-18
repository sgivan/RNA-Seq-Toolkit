#!/bin/bash

cp Sample_1/gene_cnts.txt joined.txt
for file in `ls -v Sample_{2..5}/gene_cnts.txt`; do echo $file; join --header joined.txt $file > joined2.txt; mv joined2.txt joined.txt; done
for file in `ls -v Sample_{6..10}/gene_cnts.txt`; do echo $file; join --header joined.txt $file > joined2.txt; mv joined2.txt joined.txt; done
mv joined.txt C_v_E.txt

#cp Sample_1/gene_cnts.txt joined.txt
#for file in `ls -v Sample_{2..5}/gene_cnts.txt`; do echo $file; join --header joined.txt $file > joined2.txt; mv joined2.txt joined.txt; done
#for file in `ls -v Sample_{11..15}/gene_cnts.txt`; do echo $file; join --header joined.txt $file > joined2.txt; mv joined2.txt joined.txt; done
#mv joined.txt C_v_R.txt

