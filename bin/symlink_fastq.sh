#!/bin/sh

dir=$(pwd | sed -E 's/.+\/Sample_//');
cnt=0;
for file in ../original_files/*S${dir}_*;
do
        echo $file;
        ln -fs $file ./;
        nopathname=$(echo $file | sed -r 's/.+\/(.+)$/\1/')
        echo "nopathname :'"$nopathname"'"
        let "cnt++";
        ln -s $nopathname  set${cnt}.fq;
done
