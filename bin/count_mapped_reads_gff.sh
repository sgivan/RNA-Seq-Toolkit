#!/bin/bash
file=$1

if [ $# != 1 ] ; then
    echo "what file?"
    exit
fi

cat $file | cut -f 9 | cut -f 2 -d " " | sort | uniq | wc -l

