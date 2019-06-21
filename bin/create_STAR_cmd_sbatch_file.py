#!/usr/bin/env python

import os, sys
import argparse
from string import Template

parser = argparse.ArgumentParser()
parser.add_argument("--memory", help="amount of ram; ie '25G'", type=str, default='50G')
parser.add_argument("--threads", help="number of threads; ie '5'", type=int, default=5)
parser.add_argument("--template", help="name of template file", type=str, default="STARcmdSE.t")
parser.add_argument("--index", help="name of index directory; ie, 'index'", type=str, default="index")
parser.add_argument("--minIntronLength", help="minimum allowable intron length", type=int, default=21)
parser.add_argument("--maxIntronLength", help="maximum allowable intron length", type=int, default=500000)
parser.add_argument("--alignGapMax", help="maximum allowable gap between pairs of reads", type=int, default=0)
parser.add_argument("--queue", help="cluster queue to submit job", type=str, default='shortQ')
#parser.add_argument("--inputFiles", help="input file string; ie, read_1 read_2", type=str, default="read_1 read_2")


args = parser.parse_args()

wd = os.path.dirname(os.path.realpath(__file__))

template_file = open(wd + "/../templates/" + args.template, 'r').read()

s = Template(template_file)

memory='78G'

new_t = s.substitute(mem=args.memory, threads=args.threads, index=args.index, minIntronLength=args.minIntronLength, maxIntronLength=args.maxIntronLength, alignGapMax=args.alignGapMax, clusterQ=args.queue)

print new_t
