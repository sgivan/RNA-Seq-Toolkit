#!/usr/bin/env python

import os, sys, argparse
import gffutils, string

parser = argparse.ArgumentParser()
parser.add_argument('--datafile', help='full path to data file', default='datafile')
parser.add_argument('--gfffile', help='full path to gtf file', default='gfffile')
parser.add_argument('--verbose', help='verbose output to terminal', action="store_true")
args = parser.parse_args()

outfilename = 'outfile'

if (os.path.exists(args.datafile)):
    datafile = open(args.datafile, 'r')
else:
    print(args.datafile + " doesn't exist")
    sys.exit(2)
    

if (not os.path.exists(args.gfffile)):
    print(args.gfffile + " doesn't exist")
    sys.exit(2)

outfile = open(outfilename, 'w')

if (not os.path.exists('gff.sqlite')): 
    db = gffutils.create_db(args.gfffile,'gff.sqlite')
else:
    db = gffutils.FeatureDB('gff.sqlite')

if (args.verbose):
    print("number of transcripts: %i" % db.count_features_of_type('transcript'))

cnt, cnt2 = 0, 0
gene2name = { }

for gene in db.features_of_type('gene'):

    try:
        gene.attributes['Name'][0] 
    except Exception as e: 
        cnt += 1
        continue

    gene2name[gene.id] = gene.attributes['Name'][0] 

    cnt2 += 1
    if (args.verbose and cnt2 <= 10):
        print(gene.id)

if (args.verbose):
    print("%i genes had no 'Name' attribute" % cnt)

cnt = 0
for dataline in datafile:

    geneid = dataline.split()[1]
    print("geneid: '%s'" % geneid)
    cnt += 1
    if (cnt >= 10):
        break
#    outfile.write(dataline.rstrip() + "\t" + gene2name[geneid] + "\n")

    gene_name = 'n/a'

    try:
        gene_name = gene2name[geneid]
    except:
        continue

    print(dataline.rstrip() + "\t" + gene_name + "\n")
    outfile.write(dataline.rstrip() + "\t" + gene_name + "\n")

#outfile.write("hello world")

outfile.close()

