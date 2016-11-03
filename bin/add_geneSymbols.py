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

if (not os.path.exists('gxf.sqlite')): 
    #db = gffutils.create_db(args.gfffile,'gxf.sqlite', gtf_subfeature='exon', id_spec=None, merge_strategy="merge", disable_infer_transcripts=True, disable_infer_genes=True)
    db = gffutils.create_db(args.gfffile,'gxf.sqlite', keep_order=True)
else:
    db = gffutils.FeatureDB('gxf.sqlite')

if (args.verbose):
    print("number of transcripts: %i" % db.count_features_of_type('transcript'))
    print("number of genes: %i" % db.count_features_of_type('gene'))

cnt, cnt2 = 0, 0
gene2name = { }

#for ft in db.featuretypes():
#    print("feature type: '%s'\n" % ft)

for gene in db.features_of_type('gene'):

#    print("attribute: '%s'" % gene.attributes)

    try:
        gene.attributes['Name'][0] 
    except KeyError: 
        cnt += 1
        continue
    else:
        gene2name[gene.id] = gene.attributes['Name'][0] 


    cnt2 += 1
    if (args.verbose and cnt2 <= 10):
        print(gene.id)

if (args.verbose):
    print("%i genes had no 'Name' attribute" % cnt)

cnt = 0
for dataline in datafile:

    geneid = dataline.split()[1]
    gene_name = 'n/a'

#   In the following exception handler, you want to pass  even
#   if an exception is thrown, because some genes may not have a
#   gene name, but the data is still valid. So, use pass instead
#   of continue after KeyError.
    try:
        gene_name = gene2name[geneid]
    except KeyError:
        pass

    if (args.verbose):
        print(dataline.rstrip() + "\t" + gene_name)

    outfile.write(dataline.rstrip() + "\t" + gene_name + "\n")

outfile.close()

