#!/bin/bash

for dir
do
    echo $dir
    cd $dir
    rm -fr Sample* cutadapt* finished preprocess* read_1* set1.fq set1_noadapters.fq read_2* set2.fq  set2_noadapters.fq
    ln -sf set1_00.fq set1.fq
    ln -sf set2_00.fq set2.fq
    cd ..
done


