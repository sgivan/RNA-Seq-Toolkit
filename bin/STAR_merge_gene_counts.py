#!/usr/bin/env python

import os, sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--pe_file", help="name of gene counts file generated with Paired-End data", type=str, default="pe_gene_counts.txt")
parser.add_argument("--se_file", help="name of gene counts file generated with Paired-End data", type=str, default="se_gene_counts.txt")

args = parser.parse_args()

#print "pe gene counts file is '{}'".format(args.pe_file)
#print "se gene counts file is '{}'".format(args.se_file)

#pe_file = open(args.pe_file, 'r').read()
#se_file = open(args.se_file, 'r').read()
out_file = open("merged_ReadsPerGene.out.tab", 'w')

pe = {}
se = {}

with open(args.pe_file, 'r') as pe_file:
    for pe_line in pe_file:
        vals = pe_line.rstrip().split("\t")
        pe[vals[0]] = [vals[1], vals[2], vals[3]]

with open(args.se_file, 'r') as se_file:
    for se_line in se_file:
        vals = se_line.rstrip().split("\t")
        se[vals[0]] = [vals[1], vals[2], vals[3]]


pe_se = {}

for key in pe:
    pe_se[key] = [
            int(pe[key][0]) + int(se[key][0]),
            int(pe[key][1]) + int(se[key][1]),
            int(pe[key][2]) + int(se[key][2])
            ]


#print pe_se

out_file.write( "{}\t{}\n".format( 'N_unmapped', "\t".join(map(str,pe_se['N_unmapped']))))
del pe_se['N_unmapped']
out_file.write( "{}\t{}\n".format( 'N_multimapping', "\t".join(map(str,pe_se['N_multimapping']))))
del pe_se['N_multimapping']
out_file.write( "{}\t{}\n".format( 'N_noFeature', "\t".join(map(str,pe_se['N_noFeature']))))
del pe_se['N_noFeature']
out_file.write( "{}\t{}\n".format( 'N_ambiguous', "\t".join(map(str,pe_se['N_ambiguous']))))
del pe_se['N_ambiguous']

for key in pe_se:
    #print "{}\t{}".format( key, "\t".join(map(str,pe_se[key])) )
    out_file.write( "{}\t{}\n".format( key, "\t".join(map(str,pe_se[key])) ))

#print "OK"
