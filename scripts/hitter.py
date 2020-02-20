#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
if len(sys.argv[1:]) == 0:
    print(
"""
Build a table of annotated sequences and their mapped read abundance.
usage: python3 hitter.py annotation.txt hitfile.txt sample
Miguel Romero github.com/romeromig
"""
)
    sys.exit()

infiles = sys.argv[1:]

print('Loading files.')

# load annotation table
handle = open(infiles[0])
annot = handle.read().split('\n')[:-1]
handle.close()

# each sequence points to its annotation
annot_dict = {} 
for line in annot:
    linesep = line.split('\t')
    annot_dict[linesep[0]] = linesep[-1]

# load read abundance table
handle = open(infiles[1])
hits = handle.read().split('\n')[:-1]
handle.close()

# each sequence points to its abundance
hit_dict = {}
for line in hits:
    linesep = line.split(' ')[-2:]
    hit_dict[linesep[1]] = int(linesep[0])

print('Looking for annotated and mapped sequences.')
chidas = []
for seqid in list(annot_dict.keys()):
    if seqid in list(hit_dict.keys()):
        chidas.append(seqid)

print('Building table.')

# create a row for each reference ID
annotcount_dict = {} 
for seqid in chidas:
    annotcount_dict[annot_dict[seqid]] = 0

print('Adding read counts')

for seqid in chidas:
    annotcount_dict[annot_dict[seqid]] += hit_dict[seqid]

print('Sorting')

rownames = list(annotcount_dict.keys())
rownames.sort()

print('Writing table.')

# saving in a text file
outfile = open(infiles[2]+'.hout','w') 
for annotid in rownames:
    outfile.write(annotid +'\t'+ str(annotcount_dict[annotid]) +'\n')
outfile.close()

print('Done. Final table saved in '+ infiles[2] +'.hout')
