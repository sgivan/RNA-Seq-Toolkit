#!/usr/bin/env python

import os, sys
import argparse
from string import Template

parser = argparse.ArgumentParser()
parser.add_argument("--memory", help="amount of ram [default=50G]", type=str, default='50G')
parser.add_argument("--template", help="name of template file [default=DESeq2.R]", type=str, default="DESeq2.R")
parser.add_argument("--numberOfControls", help="the number of control samples", type=str, default=3)
parser.add_argument("--numberOfExperimentals", help="the number of experimental samples", type=str, default=3)
parser.add_argument("--prefix", help="prefix for output files [default=output_]", type=str, default="output_")
parser.add_argument("--datafile", help="name of data file containing gene count data [default=gene_cnt_matrix.tab]", type=str, default="gene_cnt_matrix.tab")

args = parser.parse_args()

wd = os.path.dirname(os.path.realpath(__file__))

template_file = open(wd + "/../templates/" + args.template, 'r').read()

s = Template(template_file)

new_t = s.substitute(memory=args.memory, cntCont=args.numberOfControls, cntExp=args.numberOfExperimentals, prefix=args.prefix, datafile=args.datafile)

print new_t
