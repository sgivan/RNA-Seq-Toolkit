#!/usr/bin/env python

import os, sys
import argparse
from string import Template

parser = argparse.ArgumentParser()
parser.add_argument("--memory", help="amount of ram; ie '25G'", type=str, default='50G')
parser.add_argument("--threads", help="number of threads; ie '5'", type=int, default=5)
parser.add_argument("--template", help="name of template file", type=str, default="STARcmdSE.t")

args = parser.parse_args()

wd = os.path.dirname(os.path.realpath(__file__))

template_file = open(wd + "/../templates/" + args.template, 'r').read()

s = Template(template_file)

memory='78G'

new_t = s.substitute(mem=args.memory, threads=args.threads)

print new_t
