#!/usr/bin/env python

import os, sys, argparse
import re, yaml

# start with command line options
argparser = argparse.ArgumentParser(description="Parse tab-delimited file")

argparser.add_argument("--infile", type=str,  help="file to parse", default="file")
argparser.add_argument("--verbose", action="store_true", help="verbose messages to terminal")

args = argparser.parse_args()

if args.verbose: print 'input file: %(filename)s' % { "filename": args.infile }

f_yaml = open(args.infile, 'r') # open infile

config = yaml.load(f_yaml)

if args.verbose: print "dumping config"
print yaml.dump(config)

#
# create new working directory
# fail if the directory already exists
#
if not os.access(config['working_datadir'], os.F_OK):
    if args.verbose: print "creating new directory to place renamed files: %(newdir)s." % { "newdir": config['working_datadir'] }
    if os.access(os.path.split(config['working_datadir'])[0], os.W_OK):
        os.mkdir(config['working_datadir'])
    else:
        print "can't create the directory"
        sys.exit(2)
else:
    print "%(newdir)s already exists. Will not overwrite -- please rename or move the diretory." % { "newdir": config['working_datadir'] }
    print "Exiting now."
    sys.exit(3)

#
# end of working directory section
#

elength=len(config['input']['experimental'])
sample_number=0
if args.verbose: print "experimental data sets: %(length)i" % { "length": elength }
for i in range(0,elength):
#   This cut corresponds to the sample replicate. There can be any number of sample replicates. which will have either a single file (non-PE) or a pair of files (Paired End)
    print "sample replicates in set %(eset)i: %(filenames)s" % { "eset": i, "filenames": config['input']['experimental'][i] }
    number_of_reps=len(config['input']['experimental'][i])
    print "number of replicates: %(numseqs)i." % { "numseqs": number_of_reps }
    for j in range(0,number_of_reps):
        number_of_seq_files=len(config['input']['experimental'][i][j])
        print "number of sequences in %(seqfiles)s: %(numseqfiles)i" % { "seqfiles": config['input']['experimental'][i][j], "numseqfiles": number_of_seq_files }
        if number_of_seq_files == 2:
            if args.verbose: print "\tworking with paired-end data"
        else:
            if args.verbose: print "\tworking with non-paired-end data"

        sample_number += 1
        if args.verbose: print "\tThis replicate will be given symbolic name 'Sample_%(sint)s'" % { "sint": sample_number }


print str(sample_number) + ' samples'
print "OK"

