#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
if len(sys.argv[1:]) == 0:
    print(
"""
Merge tables built with the hitter.py script.
usage: python3 hitter_table.py hout_names_list project
Miguel Romero github.com/romeromig
"""
)
    sys.exit()

import os

argvs = sys.argv[1:]

handle = open(argvs[0],'r')
infiles = handle.read().split('\n')[:-1]
handle.close()

os.system('cat '+ ' '.join(infiles) +" | awk '{print $1}' | sort | uniq > "+ argvs[1] +'.tmp1')

placeholder = []
for i in range(len(infiles)):
    placeholder.append('0')

handle = open(argvs[1] +'.tmp1','r')
md5s = handle.read().split('\n')[:-1]
handle.close()

proyect_dict = {} # final table
for md5 in md5s:
    proyect_dict[md5] = placeholder.copy()

print('Building table.')
index = 0
for infile in infiles:
    print('Loading file: '+ infile)
    handle = open(infile,'r')
    sample = handle.read().split('\n')[:-1]
    handle.close()
    for line in sample:
        linesep = line.split('\t')
        proyect_dict[linesep[0]][index] = linesep[1]
    index += 1

outfile = open(argvs[1]+'.tsv','w') # save in a text file
outfile.write('md5\t'+ '\t'.join(infiles) +'\n')
for md5 in proyect_dict.keys():
    outfile.write(md5 +'\t'+ '\t'.join(proyect_dict[md5]) +'\n')
outfile.close()

os.system('rm '+ argvs[1] +'.tmp1')

print('Done. Final table saved in '+ argvs[1] +'.tsv')
