#!/bin/bash

for dir
do
    echo $dir
    cd $dir
    rm -fr Sample* cutadapt* finished preprocess* read_1.1 set1.fq set1_noadapters.fq
    ln -sf set1_00.fq set1.fq
    cd ..
done


