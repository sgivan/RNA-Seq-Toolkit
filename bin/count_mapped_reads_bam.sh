#!/bin/bash
samtools view merged.bam | cut -f 1 | psort --parallel=4 | uniq | wc -l > merged.cnt
