#!/usr/bin/env python

import os, sys, argparse
import re, yaml

# start with command line options
argparser = argparse.ArgumentParser(description="Parse tab-delimited file")

argparser.add_argument("--infile", type=str,  help="file to parse", default="file")
argparser.add_argument("--verbose", action="store_true", help="verbose messages to terminal")

args = argparser.parse_args()

print 'input file: %(filename)s' % { "filename": args.infile }

f_yaml = open(args.infile, 'r') # open infile

config = yaml.load(f_yaml)

print "dumping config"
print yaml.dump(config)

#print config['input']['control'][0]
print "first experimental file of first experimental set: '%(e1)s'." % { "e1": config['input']['experimental'][0][0] }

elength=len(config['input']['experimental'])
print "experimental data sets: %(length)i" % { "length": elength }
for i in range(0,elength):
    print "experimental file set %(eset)i: %(filenames)s" % { "eset": i, "filenames": config['input']['experimental'][i] }
#for fileset in config['input']['experimental']:
#    print "file set: '%(setname)s'" % { "setname": fileset }

print "OK"
